#!/usr/bin/env python3
import argparse
import csv
import re
import gzip
import tempfile
import os
import shutil
from typing import Dict, Tuple, Optional

def decompress_gzip(path: str) -> Tuple[str, bool]:
    """
    If the file is gzip-compressed (.gz), decompress to a temp file and return
    (decompressed_path, True). Otherwise return (original_path, False).
    """
    if not path.endswith(".gz"):
        return path, False
    tmp_fd, tmp_path = tempfile.mkstemp(suffix=".fa")
    os.close(tmp_fd)
    with gzip.open(path, "rb") as f_in, open(tmp_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
        return tmp_path, True

def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]

def parse_args():
    p = argparse.ArgumentParser(
        description="Extract <upstream> b upstream sequences for Ensembl gene IDs using genome FASTA + GTF."
    )
    p.add_argument("-m","--mapping", required=True, help="CSV/TSV with columns: gene_id, ensmbl_id")
    p.add_argument("--ens-col", default="ensmbl_id",
                   help="Mapping column to use as Ensembl ID key (default: ensmbl_id)")
    p.add_argument("--gene-col", default="gene_id",
                   help="Mapping column to use as gene label in FASTA header (default: gene_id)")
    p.add_argument("--gtf", required=True, help="GTF file (galGal6 assembly)")
    p.add_argument("--genome", required=True, help="Genome FASTA (galGal6.fa)")
    p.add_argument("--out", required=True, help="Output FASTA")
    p.add_argument("--upstream", type=int, default=2000, help="Upstream length (default 2000)")
    p.add_argument("--delim", default=None,
                   help="Delimiter for mapping file (default: auto-detect comma vs tab)")
    p.add_argument("--feature", default="gene",
                   help="GTF feature to use for coordinates (default: gene; fallback: transcript)")
    p.add_argument("--padN", action="store_true",
                   help="Pad with Ns to guarantee exactly upstream length even at chromosome ends")
    p.add_argument("--debug", action="store_true",
                   help="Print debug info about mapping/GTF columns, IDs, and match counts")
    return p.parse_args()


def normalize_ensembl_id(ens_id: str) -> str:
    """
    Normalize Ensembl-like IDs for matching by removing optional version suffix,
    e.g. ENSGALG00000012345.1 -> ENSGALG00000012345.
    """
    return ens_id.split(".", 1)[0].strip()


def extract_gtf_id(attrs: str) -> Optional[str]:
    """Extract and normalize a gene-like ID from GTF attributes."""
    re_ens_id_q = re.compile(r'ens(?:embl|mbl)_id "([^"]+)"')
    re_ens_id_eq = re.compile(r'ens(?:embl|mbl)_id=([^;]+)')
    re_gene_id_q = re.compile(r'gene_id "([^"]+)"')
    re_gene_id_eq = re.compile(r'gene_id=([^;]+)')
    re_gene_q = re.compile(r'gene "([^"]+)"')

    gid = None
    m = re_ens_id_q.search(attrs)
    if m:
        gid = m.group(1).strip()
    else:
        m = re_ens_id_eq.search(attrs)
        if m:
            gid = m.group(1).strip()
    if gid is None:
        m = re_gene_id_q.search(attrs)
        if m:
            gid = m.group(1).strip()
        else:
            m = re_gene_id_eq.search(attrs)
            if m:
                gid = m.group(1).strip()
            else:
                m = re_gene_q.search(attrs)
                if m:
                    gid = m.group(1).strip()

    if gid:
        return normalize_ensembl_id(gid)
    return None


def get_gtf_id_debug_info(
    gtf_path: str,
    needed_ens_ids: set,
    feature: str,
    sample_limit: int = 5,
) -> Dict[str, object]:
    scanned_rows = 0
    unique_gtf_ids = set()
    overlap_ids = set()
    with open(gtf_path, "r") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            _, _, feat, _, _, _, _, _, attrs = fields
            if feat != feature:
                continue
            scanned_rows += 1
            gid = extract_gtf_id(attrs)
            if not gid:
                continue
            unique_gtf_ids.add(gid)
            if gid in needed_ens_ids:
                overlap_ids.add(gid)

    return {
        "feature": feature,
        "scanned_rows": scanned_rows,
        "unique_gtf_ids": len(unique_gtf_ids),
        "sample_gtf_ids": sorted(list(unique_gtf_ids))[:sample_limit],
        "overlap_count": len(overlap_ids),
        "sample_overlap_ids": sorted(list(overlap_ids))[:sample_limit],
    }

def detect_delim(path:str) -> str:
    with open(path, "r", newline="") as f:
        sample = f.read(4096)
    if "\t" in sample and sample.count("\t") >= sample.count(","):
        return "\t"
    return ","

def load_mapping(
    path: str,
    delim: Optional[str],
    ens_col_arg: str,
    gene_col_arg: str,
) -> Tuple[Dict[str, str], Dict[str, object]]:
    delim = delim or detect_delim(path)
    mapping = {}
    rows_seen = 0
    rows_with_ens_id = 0
    rows_with_gene_and_ens = 0
    headers = {}
    gene_col = ""
    ens_col = ""
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter=delim)
        if not reader.fieldnames:
            raise ValueError("Mapping file appears to have no header row.")

        headers = {h.lower(): h for h in reader.fieldnames}
        ens_col_name = ens_col_arg.strip().lower()
        gene_col_name = gene_col_arg.strip().lower()
        if gene_col_name not in headers or ens_col_name not in headers:
            raise ValueError(
                "Mapping file is missing requested columns. "
                f"Requested gene column='{gene_col_arg}', ensembl column='{ens_col_arg}'. "
                f"Found: {reader.fieldnames}"
            )

        gene_col = headers[gene_col_name]
        ens_col = headers[ens_col_name]

        for row in reader:
            rows_seen += 1
            gene_id = (row.get(gene_col) or "").strip()
            ens_id_raw = (row.get(ens_col) or "").strip()
            ens_id = normalize_ensembl_id(ens_id_raw)
            if ens_id:
                rows_with_ens_id += 1
                if gene_id:
                    rows_with_gene_and_ens += 1
                mapping[ens_id] = gene_id or ens_id
    debug_info = {
        "delimiter": delim,
        "fieldnames": list(headers.values()),
        "gene_col": gene_col,
        "ens_col": ens_col,
        "rows_seen": rows_seen,
        "rows_with_ens_id": rows_with_ens_id,
        "rows_with_gene_and_ens": rows_with_gene_and_ens,
        "unique_ens_ids": len(mapping),
    }
    return mapping, debug_info


def parse_gtf_coords(
    gtf_path: str,
    needed_ens_ids: set,
    feature: str = "gene",
) -> Dict[str, Tuple[str, int, int, str]]:
    """
    Return dict: ens_gene_id -> (chrom, start_1based, end_1based, strand)

    We try common Ensembl-ish attribute keys used in GTF/GFF:
      - ensembl_id "ENSGALG..." / ensmbl_id "ENSGALG..."
      - ensembl_id=ENSGALG... / ensmbl_id=ENSGALG...
      - gene_id "ENSGALG..."
      - gene_id=ENSGALG...
      - gene "ENSGALG..." (rare)
    """
    out: Dict[str,Tuple[str,int,int,str]] = {}

    with open(gtf_path, "r") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _, feat, start, end, _, strand, _, attrs = fields
            if feat != feature:
                continue

            gid = extract_gtf_id(attrs)

            if not gid or gid not in needed_ens_ids:
                continue

            out[gid] = (chrom, int(start), int(end), strand)
    return out

def normalize_chrom(chrom:str, genome_keys: set) -> Optional[str]:
    if chrom in genome_keys:
        return chrom
    if chrom.startswith("chr"):
        alt = chrom.replace("chr", "", 1)
        if alt in genome_keys:
            return alt
    else:
        alt = "chr" + chrom
        if alt in genome_keys:
            return alt
    return None

def pad_to_length(seq:str, desired:int, pad_left:int,pad_right:int) -> str:
    # pad_left/right: number of missing bases on left/right side in genomic sense (not strand)
    # We just pad with Ns on the appropriate side to make the final length exactly desired.
    if len(seq) >= desired:
        return seq[:desired]
    return ("N" * pad_left) + seq + ("N"*pad_right)

def main():
    args = parse_args()

    try:
        from pyfaidx import Fasta
    except ImportError:
        raise SystemExit("Install dependency: pip install pyfaidx")

    mapping, mapping_dbg = load_mapping(
        args.mapping,
        args.delim,
        args.ens_col,
        args.gene_col,
    )
    needed = set(mapping.keys())

    if args.debug:
        print("[DEBUG] Mapping delimiter:", repr(mapping_dbg["delimiter"]))
        print("[DEBUG] Mapping headers:", mapping_dbg["fieldnames"])
        print("[DEBUG] Mapping gene column:", mapping_dbg["gene_col"])
        print("[DEBUG] Mapping Ensembl column:", mapping_dbg["ens_col"])
        print("[DEBUG] Mapping rows seen:", mapping_dbg["rows_seen"])
        print("[DEBUG] Mapping rows with Ensembl ID:", mapping_dbg["rows_with_ens_id"])
        print("[DEBUG] Mapping rows with both IDs:", mapping_dbg["rows_with_gene_and_ens"])
        print("[DEBUG] Mapping unique normalized Ensembl IDs:", mapping_dbg["unique_ens_ids"])
        print("[DEBUG] Mapping sample Ensembl IDs:", list(needed)[:5])

        gtf_dbg = get_gtf_id_debug_info(args.gtf, needed, feature=args.feature)
        print("[DEBUG] GTF debug feature:", gtf_dbg["feature"])
        print("[DEBUG] GTF rows scanned for feature:", gtf_dbg["scanned_rows"])
        print("[DEBUG] GTF unique extracted IDs:", gtf_dbg["unique_gtf_ids"])
        print("[DEBUG] GTF sample extracted IDs:", gtf_dbg["sample_gtf_ids"])
        print("[DEBUG] GTF overlap with mapping IDs:", gtf_dbg["overlap_count"])
        print("[DEBUG] GTF sample overlap IDs:", gtf_dbg["sample_overlap_ids"])

    coords = parse_gtf_coords(args.gtf, needed, feature=args.feature)
    if len(coords) == 0 and args.feature == "gene":
        coords = parse_gtf_coords(args.gtf, needed, feature="transcript")

    if args.debug:
        print("[DEBUG] GTF matched IDs:", len(coords))
        if coords:
            print("[DEBUG] GTF sample matched IDs:", list(coords.keys())[:5])

    genome_path,  decompressed = decompress_gzip(args.genome)
    try:
        genome = Fasta(genome_path, as_raw=True, sequence_always_upper=True)
        genome_keys = set(genome.keys())

        missing = 0
        written = 0

        with open(args.out, "w") as out_fa:
            for ens_id, gene_id in mapping.items():
                if ens_id not in coords:
                    missing += 1
                    continue

                chrom, start1, end1, strand = coords[ens_id]
                chrom_norm = normalize_chrom(chrom, genome_keys)
                if chrom_norm is None:
                    missing += 1
                    continue
                chrom = chrom_norm

                chrom_len = len(genome[chrom])

                if strand == "+":
                    tss = start1
                    raw_start1 = tss - args.upstream
                    raw_end1 = tss - 1
                    if raw_end1 < 1:
                        if args.padN:
                            seq = "N" * args.upstream
                        else:
                            continue
                    else:
                        up_start1 = max(1, raw_start1)
                        up_end1 = max(1, raw_end1)
                        seq = genome[chrom][up_start1 - 1: up_end1]
                        if args.padN and raw_start1 < 1:
                            pad_left = 1-raw_start1
                            seq = pad_to_length(seq,args.upstream, pad_left,pad_right=0)

                    header = f">{gene_id}|{ens_id}|{chrom}:{max(1, raw_start1)}-{max(1, raw_end1)}({strand})"

                elif strand == "-":
                    tss = end1
                    raw_start1 = tss + 1
                    raw_end1 = tss + args.upstream

                    if raw_start1 > chrom_len:
                        if args.padN:
                            seq = "N" * args.upstream
                        else:
                            continue

                    else:
                        up_start1 = min(chrom_len, raw_start1)
                        up_end1 = min(chrom_len, raw_end1)
                        seq_raw = genome[chrom][up_start1 - 1:up_end1]
                        seq = revcomp(seq_raw)
                        if args.padN and raw_end1 > chrom_len:
                            pad_right = raw_end1 - chrom_len
                            seq = pad_to_length(seq, args.upstream, pad_left=0, pad_right=pad_right)

                    header = f">{gene_id}|{ens_id}|{chrom}:{min(chrom_len, raw_start1)}-{min(chrom_len, raw_end1)}({strand})"

                else:
                    missing += 1
                    continue
                out_fa.write(header + "\n")
                for i in range(0, len(seq), 60):
                    out_fa.write(seq[i:i+60] + "\n")
                written += 1
    finally:
        if decompressed:
            for ext in ("", ".fai"):
                tmp = genome_path + ext
                if os.path.exists(tmp):
                    os.remove(tmp)

    print("Done.")
    print(f"Written: {written}")
    print(f"Missing (not in GTF or chrom mismatch): {missing}")

if __name__ == "__main__":
    main()