#!/usr/bin/env python3
"""Find motif-pattern occurrences in nucleotide FASTA files.

Motifs are read from a FASTA-like text file where each record contains a motif
name on the header line and the regex pattern on subsequent line(s), for
example:

- ISRE_core
  YAGTTTC(A/T)YTTTYCC
- ISRE_alt
  AGTTTCNNTTTCNC/T
- ISRE_spaced
  A/GGTTTCN(1-2)TTTCC/T

Each pattern can be raw regex or shorthand using IUPAC symbols and slash
alternatives, for example:
- YAGTTTC(A/T)YTTTYCC
- AGTTTCNNTTTCNC/T
- A/GGTTTCN(1-2)TTTCC/T

By default, searches both forward and reverse-complement strands.
"""

import argparse
from dataclasses import dataclass
from pathlib import Path
import re
from typing import Dict, Iterator, List, Sequence, Set, Tuple


IUPAC_BASES = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "U": "T",
    "R": "AG",
    "Y": "CT",
    "S": "GC",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}


@dataclass
class Match:
    seq_id: str
    strand: str
    start_1based: int
    end_1based: int
    motif_label: str
    matched_seq: str


@dataclass(frozen=True)
class NamedMotif:
    name: str
    pattern: str


def parse_fasta(path: str) -> Iterator[Tuple[str, str]]:
    """Yield (header_id, sequence) records from FASTA."""
    seq_id = None
    chunks: List[str] = []

    with open(path, "r") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_id is not None:
                    yield seq_id, "".join(chunks).upper().replace("U", "T")
                seq_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)

    if seq_id is not None:
        yield seq_id, "".join(chunks).upper().replace("U", "T")


def revcomp(seq: str) -> str:
    trans = str.maketrans("ACGTRYMKBDHVN", "TGCAYRKMVHDBN")
    return seq.translate(trans)[::-1]


def _bases_to_regex_atom(bases: str) -> str:
    ordered = "".join(base for base in "ACGT" if base in set(bases))
    if len(ordered) == 1:
        return ordered
    return f"[{ordered}]"


def _symbols_to_regex_atom(symbols: Sequence[str]) -> str:
    bases = "".join(IUPAC_BASES.get(symbol, symbol) for symbol in symbols)
    return _bases_to_regex_atom(bases)


def motif_to_regex(motif: str) -> str:
    """Convert shorthand motifs into regex while allowing raw regex input."""
    text = motif.strip().upper().replace(" ", "")
    if not text:
        return text

    out: List[str] = []
    i = 0
    while i < len(text):
        ch = text[i]

        if ch.isalpha():
            symbols = [ch]
            j = i
            while j + 2 < len(text) and text[j + 1] == "/" and text[j + 2].isalpha():
                symbols.append(text[j + 2])
                j += 2

            atom = _symbols_to_regex_atom(symbols)
            i = j + 1

            if i < len(text):
                m = None
                if text[i] == "(":
                    m = re.match(r"\((\d+)-(\d+)\)", text[i:])
                else:
                    m = re.match(r"(\d+)-(\d+)", text[i:])
                if m:
                    atom = f"{atom}{{{m.group(1)},{m.group(2)}}}"
                    i += len(m.group(0))

            out.append(atom)
            continue

        if ch == "(":
            close = text.find(")", i + 1)
            if close == -1:
                raise ValueError(f"Unclosed parenthesis in motif: {motif}")

            inside = text[i + 1 : close]
            if re.fullmatch(r"[A-Z](?:/[A-Z])+", inside):
                out.append(_symbols_to_regex_atom(inside.split("/")))
            else:
                out.append(f"({inside})")
            i = close + 1
            continue

        out.append(ch)
        i += 1

    return "".join(out)


def load_motifs(path: str) -> List[NamedMotif]:
    """Load FASTA-like motif records as (header, pattern) pairs."""
    motifs: List[NamedMotif] = []
    name = None
    pattern_chunks: List[str] = []

    with open(path, "r") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    pattern = "".join(pattern_chunks).strip()
                    if not pattern:
                        raise ValueError(f"Motif '{name}' has no pattern in {path}")
                    motifs.append(NamedMotif(name=name, pattern=pattern))
                name = line[1:].strip()
                if not name:
                    raise ValueError(f"Encountered empty motif header in {path}")
                pattern_chunks = []
                continue
            if name is None:
                raise ValueError(f"Expected motif header starting with '>' before pattern in {path}")
            pattern_chunks.append(line)

    if name is not None:
        pattern = "".join(pattern_chunks).strip()
        if not pattern:
            raise ValueError(f"Motif '{name}' has no pattern in {path}")
        motifs.append(NamedMotif(name=name, pattern=pattern))

    if not motifs:
        raise ValueError(f"No motifs found in motif file: {path}")
    return motifs


def find_overlapping(sequence: str, motif: str) -> Iterator[Tuple[int, str]]:
    """Yield 0-based start and matched sequence for overlapping motif hits."""
    regex = re.compile(f"(?=({motif_to_regex(motif)}))")
    for m in regex.finditer(sequence):
        yield m.start(), m.group(1)


def safe_output_name(name: str) -> str:
    """Convert a motif header into a safe basename for output files."""
    sanitized = name.strip().replace("/", "_").replace("\\", "_")
    if not sanitized:
        raise ValueError("Motif header cannot be empty")
    return sanitized


def find_matches(
    seq_id: str,
    sequence: str,
    forward_motifs: Sequence[NamedMotif],
    search_reverse: bool,
) -> List[Match]:
    seq_len = len(sequence)
    out: List[Match] = []

    for motif in forward_motifs:
        for start0, matched in find_overlapping(sequence, motif.pattern):
            out.append(
                Match(
                    seq_id=seq_id,
                    strand="+",
                    start_1based=start0 + 1,
                    end_1based=start0 + len(matched),
                    motif_label=motif.name,
                    matched_seq=matched,
                )
            )

    if search_reverse:
        rc_sequence = revcomp(sequence)
        for motif in forward_motifs:
            for rc_start0, rc_matched in find_overlapping(rc_sequence, motif.pattern):
                match_len = len(rc_matched)
                start_1based = seq_len - (rc_start0 + match_len) + 1
                end_1based = seq_len - rc_start0
                out.append(
                    Match(
                        seq_id=seq_id,
                        strand="-",
                        start_1based=start_1based,
                        end_1based=end_1based,
                        motif_label=motif.name,
                        matched_seq=revcomp(rc_matched),
                    )
                )

    out.sort(key=lambda x: (x.seq_id, x.start_1based, x.end_1based, x.strand, x.motif_label))
    return out


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Find motif-pattern occurrences in FASTA sequences. "
            "Outputs TSV with coordinates relative to each sequence."
        )
    )
    p.add_argument(
        "-f",
        "--fasta",
        help="Input nucleotide FASTA file")
    p.add_argument(
        "-o",
        "--output",
        default=".",
        help="Output directory for per-motif TSV files (default: current directory)",
    )
    p.add_argument(
        "-m",
        "--motifs-file",
        required=True,
        help=(
            "FASTA-like file of motif records. Each entry must be '>name' "
            "followed by the regex pattern on the next line(s). Supports raw "
            "regex or shorthand with IUPAC symbols (e.g. N, Y), slash "
            "alternatives (A/G or (A/T)), and ranges like N(1-2)."
        ),
    )
    p.add_argument(
        "--forward-only",
        action="store_true",
        help="Search only forward motifs (do not scan reverse complement motifs)",
    )
    return p


def main() -> None:
    args = build_arg_parser().parse_args()
    hits: List[Match] = []
    motifs = load_motifs(args.motifs_file)
    fasta_records = list(parse_fasta(args.fasta))
    total_genes = len(fasta_records)

    for seq_id, sequence in fasta_records:
        hits.extend(
            find_matches(
                seq_id=seq_id,
                sequence=sequence,
                forward_motifs=motifs,
                search_reverse=not args.forward_only,
            )
        )

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    input_prefix = Path(args.fasta).stem
    header = "seq_id\tstrand\tstart_1based\tend_1based\tmotif\tmatched_seq"
    hits_by_motif: Dict[str, List[Match]] = {motif.name: [] for motif in motifs}
    for hit in hits:
        hits_by_motif[hit.motif_label].append(hit)

    for motif in motifs:
        motif_hits = hits_by_motif[motif.name]
        unique_genes_found = len({m.seq_id for m in motif_hits})
        lines = [header] + [
            f"{m.seq_id}\t{m.strand}\t{m.start_1based}\t{m.end_1based}\t{m.motif_label}\t{m.matched_seq}"
            for m in motif_hits
        ]
        lines.append(f"#unique_genes_found/total_genes\t{unique_genes_found}/{total_genes}")
        output_path = output_dir / f"{input_prefix}_{safe_output_name(motif.name)}_hits.tsv"
        with open(output_path, "w") as handle:
            handle.write("\n".join(lines) + "\n")

    unique_hits: List[Match] = []
    seen_hits: Set[Tuple[str, str, int, int, str, str]] = set()
    for hit in hits:
        key = (
            hit.seq_id,
            hit.strand,
            hit.start_1based,
            hit.end_1based,
            hit.motif_label,
            hit.matched_seq,
        )
        if key in seen_hits:
            continue
        seen_hits.add(key)
        unique_hits.append(hit)

    unique_hits.sort(key=lambda x: (x.seq_id, x.start_1based, x.end_1based, x.strand, x.motif_label))
    unique_genes_found = len({m.seq_id for m in unique_hits})
    all_lines = [header] + [
        f"{m.seq_id}\t{m.strand}\t{m.start_1based}\t{m.end_1based}\t{m.motif_label}\t{m.matched_seq}"
        for m in unique_hits
    ]
    all_lines.append(f"#unique_genes_found/total_genes\t{unique_genes_found}/{total_genes}")
    all_hits_path = output_dir / f"{input_prefix}_all_unique_hits.tsv"
    with open(all_hits_path, "w") as handle:
        handle.write("\n".join(all_lines) + "\n")


if __name__ == "__main__":
    main()
