#!/usr/bin/env python3
"""Validate upstream FASTA sequences against a reference genome.

Supports two FASTA header styles:
1. >GENE|ENSGALG...|chr:start-end(strand)
2. >ENSGALG...

For style 1, coordinates are taken directly from the FASTA header.
For style 2, coordinates are recovered from the GTF using the Ensembl gene ID
and the requested upstream length.
"""

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple


HEADER_WITH_COORDS_RE = re.compile(
    r"^(?P<gene>[^|]+)\|(?P<ens>[^|]+)\|(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)\((?P<strand>[+-])\)$"
)


@dataclass(frozen=True)
class FastaRecord:
    header: str
    sequence: str


def parse_fasta(path: str) -> Iterator[FastaRecord]:
    header = None
    chunks: List[str] = []
    with open(path, "r") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield FastaRecord(header=header, sequence="".join(chunks).upper())
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        yield FastaRecord(header=header, sequence="".join(chunks).upper())


def revcomp(seq: str) -> str:
    trans = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(trans)[::-1]


def normalize_chrom(chrom: str, genome_keys: set) -> Optional[str]:
    if chrom in genome_keys:
        return chrom
    if chrom.startswith("chr"):
        alt = chrom[3:]
        if alt in genome_keys:
            return alt
    else:
        alt = "chr" + chrom
        if alt in genome_keys:
            return alt
    return None


def parse_gtf_coords(gtf_path: str) -> Dict[str, Tuple[str, int, int, str]]:
    out: Dict[str, Tuple[str, int, int, str]] = {}
    re_gene_id_q = re.compile(r'gene_id "([^"]+)"')
    re_gene_id_eq = re.compile(r"gene_id=([^;]+)")
    re_gene_q = re.compile(r'gene "([^"]+)"')

    with open(gtf_path, "r") as handle:
        for raw in handle:
            if not raw or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _, feat, start, end, _, strand, _, attrs = fields
            if feat != "gene":
                continue

            gene_id = None
            for rx in (re_gene_id_q, re_gene_id_eq, re_gene_q):
                match = rx.search(attrs)
                if match:
                    gene_id = match.group(1).strip()
                    break
            if gene_id:
                out[gene_id] = (chrom, int(start), int(end), strand)
    return out


def extract_expected_sequence(
    genome,
    chrom: str,
    start1: int,
    end1: int,
    strand: str,
) -> str:
    seq = genome[chrom][start1 - 1:end1]
    if strand == "-":
        return revcomp(seq)
    return seq


def coords_from_gtf(
    ens_id: str,
    gtf_coords: Dict[str, Tuple[str, int, int, str]],
    upstream: int,
    genome,
    genome_keys: set,
) -> Optional[Tuple[str, int, int, str]]:
    if ens_id not in gtf_coords:
        return None
    chrom, gene_start, gene_end, strand = gtf_coords[ens_id]
    chrom = normalize_chrom(chrom, genome_keys)
    if chrom is None:
        return None
    chrom_len = len(genome[chrom])

    if strand == "+":
        start1 = max(1, gene_start - upstream)
        end1 = max(1, gene_start - 1)
    else:
        start1 = min(chrom_len, gene_end + 1)
        end1 = min(chrom_len, gene_end + upstream)
    return chrom, start1, end1, strand


def infer_issue(observed: str, expected: str) -> str:
    if observed == revcomp(expected):
        return "reverse-complemented"
    if len(observed) != len(expected):
        return f"length-mismatch observed={len(observed)} expected={len(expected)}"
    if observed[1:] == expected[:-1]:
        return "shifted-by-1-left"
    if observed[:-1] == expected[1:]:
        return "shifted-by-1-right"
    return "sequence-mismatch"


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Validate upstream FASTA against reference genome.")
    parser.add_argument("-f", "--fasta", required=True, help="Upstream FASTA to validate")
    parser.add_argument("-g", "--genome", required=True, help="Reference genome FASTA")
    parser.add_argument("--gtf", help="Reference GTF; required when FASTA headers do not include coordinates")
    parser.add_argument("-u", "--upstream", type=int, default=2000, help="Upstream length for coordinate recovery")
    parser.add_argument("--limit", type=int, default=20, help="Max mismatches to print")
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()

    try:
        from pyfaidx import Fasta
    except ImportError:
        raise SystemExit("Install dependency: pip install pyfaidx")

    genome = Fasta(args.genome, as_raw=True, sequence_always_upper=True)
    genome_keys = set(genome.keys())
    gtf_coords = parse_gtf_coords(args.gtf) if args.gtf else {}

    checked = 0
    matches = 0
    mismatches: List[str] = []

    for record in parse_fasta(args.fasta):
        header_match = HEADER_WITH_COORDS_RE.match(record.header)
        if header_match:
            chrom = normalize_chrom(header_match.group("chrom"), genome_keys)
            if chrom is None:
                mismatches.append(f"{record.header}\tchrom-not-found")
                continue
            start1 = int(header_match.group("start"))
            end1 = int(header_match.group("end"))
            strand = header_match.group("strand")
        else:
            ens_id = record.header.split()[0]
            if not args.gtf:
                raise SystemExit("Headers lack coordinates; provide --gtf to recover them.")
            recovered = coords_from_gtf(ens_id, gtf_coords, args.upstream, genome, genome_keys)
            if recovered is None:
                mismatches.append(f"{record.header}\tcoords-not-found")
                continue
            chrom, start1, end1, strand = recovered

        expected = extract_expected_sequence(genome, chrom, start1, end1, strand)
        checked += 1
        if record.sequence == expected:
            matches += 1
        else:
            issue = infer_issue(record.sequence, expected)
            mismatches.append(
                f"{record.header}\t{issue}\tobserved_len={len(record.sequence)}\texpected_len={len(expected)}"
            )

    print(f"Checked: {checked}")
    print(f"Exact matches: {matches}")
    print(f"Mismatches: {len(mismatches)}")
    for line in mismatches[: args.limit]:
        print(line)
    if len(mismatches) > args.limit:
        print(f"... {len(mismatches) - args.limit} more")


if __name__ == "__main__":
    main()
