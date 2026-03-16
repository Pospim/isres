#!/usr/bin/env python3
import argparse
import csv
import re
import gzip
import tempfile
import os
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
        return tmp_path, True

def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]

def parse_args():
    p = argparse.ArgumentParser(
        description="Extract 1kb upstream sequences for Ensembl gene IDs using genome FASTA + GTF."
    )
    p.add_argument("-m","--mapping", required=True, help="CSV/TSV with columns: gene_id, ensmbl_id")
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
    return p.parse_args()

def detect_delim(path:str) -> str:
    with open(path, "r", newline="") as f:
        sample = f.read(4096)
    if "\t" in sample and sample.count("\t") >= sample.count(","):
        return "\t"
    return ","

def load_mapping(path: str, delim: Optional[str]) -> Dict[str,str]:
    delim = delim or detect_delim(path)
    mapping = {}
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter=delim)
        if not reader.fieldnames:
            raise ValueError("Mapping file appears to have no header row.")

        headers = {h.lower(): h for h in reader.fieldnames}
        if "gene_id" not in headers or "ensmbl_id" not in headers:
            raise ValueError(
                f"Mapping file must have headers 'gene_id' and 'ensmbl_id'. Found: {reader.fieldnames}"
            )

        gene_col = headers["gene_id"]
        ens_col = headers["ensmbl_id"]

        for row in reader:
            gene_id = (row.get(gene_col) or "").strip()
            ens_id = (row.get(ens_col) or "").strip()
            if gene_id and ens_id:
                mapping[ens_id] = gene_id
    return mapping


def parse_gtf_coords(
    gtf_path: str,
    needed_ens_ids: set,
    feature: str = "gene",
) -> Dict[str, Tuple[str, int, int, str]]:
    """
    Return dict: ens_gene_id -> (chrom, start_1based, end_1based, strand)

    We try common attribute keys used in GTF:
      - gene_id "ENSGALG..."
      - gene_id=ENSGALG...
      - gene "ENSGALG..." (rare)
    """
    out: Dict[str,Tuple[str,int,int,str]] = {}

    re_gene_id_q = re.compile(r'gene_id "([^"]+)"')
    re_gene_id_eq = re.compile(r'gene_id=([^;]+)')
    re_gene_q = re.compile(r'gene "([^"]+)"')

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

            gid = None
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

    mapping = load_mapping(args.mapping, args.delim)
    needed = set(mapping.keys())

    coords = parse_gtf_coords(args.gtf, needed, feature=args.feature)
    if len(coords) == 0 and args.feature == "gene":
        coords = parse_gtf_coords(args.gtf, needed, feature="transcript")

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