#!/usr/bin/env python3
import sys, re, os
import argparse

def parse_fasta(path):
    seqs = {}
    with open(path, encoding="utf-8", errors="ignore") as f:
        sid, buf = None, []
        for line in f:
            if line.startswith(">"):
                if sid is not None:
                    seqs[sid] = "".join(buf).upper()
                sid = line[1:].rstrip("\n")
                buf = []
            else:
                buf.append(line.strip())
        if sid is not None:
            seqs[sid] = "".join(buf).upper()
    return seqs

def revcomp(s):
    comp = str.maketrans("ACGTRYMKBDHVNacgtrymkbdhvn", "TGCAYRKMVHDBNtgcayrkmvhdbn")
    return s.translate(comp)[::-1]

def extract_gene_info(fasta_id):
    """
    Extract gene_name and ensembl_id from FASTA ID.
    Examples:
        "ZNFX1 (ENSG00000124201) upstream 2000bp" -> ("ZNFX1", "ENSG00000124201")
        "ENSG00000124201 upstream 2000bp" -> (None, "ENSG00000124201")
        "ZNFX1 upstream 2000bp" -> ("ZNFX1", None)
    """
    # Remove leading '>' if present
    fasta_id = fasta_id.lstrip('>')

    # Try to extract pattern: GENE_NAME (ENSEMBL_ID)
    match = re.match(r'([^\s(]+)\s*\((\w+)\)', fasta_id)
    if match:
        gene_name = match.group(1)
        ensembl_id = match.group(2)
        return gene_name, ensembl_id

    # Try to extract just ENSEMBL ID
    match = re.match(r'(ENS\w+)', fasta_id)
    if match:
        return None, match.group(1)

    # Otherwise, assume it's a gene name
    gene_name = fasta_id.split()[0]
    return gene_name, None

def main():
    parser = argparse.ArgumentParser(
        description="Convert HOMER motif hits to coordinate format with matched sequences from FASTA"
    )
    parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file with sequences")
    parser.add_argument("-i", "--hits", required=True, help="HOMER motif hits file (offset format)")
    parser.add_argument("-o", "--output-dir", dest="output_dir", required=True, help="Output directory for results")
    parser.add_argument("--output-file", dest="output_file", default="homerMotifs.all_coords.tsv",
                        help="Output filename (default: homerMotifs.all_coords.tsv)")

    args = parser.parse_args()

    fa = args.fasta
    hits_txt = args.hits
    outdir = args.output_dir
    output_filename = args.output_file

    os.makedirs(outdir, exist_ok=True)

    seqs = parse_fasta(fa)
    if not seqs:
        print(f"[ERROR] No sequences found in FASTA: {fa}", file=sys.stderr)
        sys.exit(1)

    print(f"[INFO] Loaded {len(seqs)} sequences from {fa}")
    out_path = os.path.join(outdir, output_filename)

    out = open(out_path, "w", encoding="utf-8")

    out.write("\t".join(["motif", "gene_name", "ensembl_id", "start", "end", "strand", "score", "matched_sequence", "reverse"]) + "\n")

    matches_found = 0
    warnings = 0

    # Store all hits for deduplication: key = (motif, gene, start, end, strand, matched_seq), value = max score
    all_hits = {}

    with open(hits_txt, encoding="utf-8", errors="ignore") as f:
        _=f.readline()

        for line in f:
            if not line.strip():
                continue

            parts = re.split(r"\t+|\s{2,}", line.rstrip("\n"))
            if len(parts) < 6:
                parts = line.rstrip("\n").split("\t")

                if len(parts)<6:
                    continue

            fid, offset_str, motif_seq, motif, strand, score = parts[:6]

            if fid not in seqs:
                cand = None

                for k in seqs.keys():
                    if k.startswith(fid):
                        cand = k; break

                if cand is None:
                    warnings += 1
                    if warnings <= 10:  # Limit warning spam
                        print(f"[WARN] FASTA ID not found in FASTA: {fid}", file=sys.stderr)
                    continue
                fid = cand
            seq = seqs[fid]
            L = len(seq)
            try:
                offset = int(offset_str)
            except:
                continue

            center0 = L // 2
            expected0 = center0 + offset
            motif_len = len(motif_seq)

            expected0 = max(0, min(expected0, L - motif_len))

            window_lo = max(0, expected0 - 200)
            window_hi = min(L, expected0 + 200 + motif_len)
            window = seq[window_lo:window_hi]
            # The "matched_seq" string is the genomic text in the FASTA orientation.
            # On '-' strand, matched_seq should still appear as-is in the FASTA.
            # So we search for it directly.
            pos = window.find(motif_seq.upper())
            if pos >= 0:
                start0 = window_lo + pos
            else:
                # As a fallback, if strand '+' and not found, try anyway;
                # If strand '-' and not found, we're still expecting the same text in FASTA.
                # Last resort: global search (can be slower but safe).
                rc = revcomp(motif_seq.upper())
                pos3 = seq.find(rc)
                if pos3 >= 0:
                    start0 = pos3
                else:
                    warnings += 1
                    if warnings <= 10:  # Limit warning spam
                        print(f"[WARN] Could not locate match in FASTA for {fid} offset {offset} motif {motif}", file=sys.stderr)
                    continue
            start1 = start0 + 1
            end1 = start1 + motif_len - 1

            # Extract gene name and ENSEMBL ID
            gene_name, ensembl_id = extract_gene_info(fid)
            gene_name_str = gene_name if gene_name else ""
            ensembl_id_str = ensembl_id if ensembl_id else ""

            # Get the motif base sequence (after the '-' in motif name like "1-TGTAATCCCA")
            motif_base = motif.split('-', 1)[1] if '-' in motif else motif

            # Calculate reverse complement of the motif
            reverse_seq = revcomp(motif_base)

            # Write to all hits file
            out.write("\t".join([motif, gene_name_str, ensembl_id_str, str(start1), str(end1), strand, score, motif_seq, reverse_seq]) + "\n")
            matches_found += 1

            # Store for unique hits (keep highest score for duplicate entries)
            key = (motif, gene_name_str, ensembl_id_str, start1, end1, strand, motif_seq, reverse_seq)
            score_float = float(score)
            if key not in all_hits or score_float > all_hits[key]:
                all_hits[key] = score_float

        out.close()

    print(f"[INFO] Wrote {matches_found} motif matches to {out_path}")

    # Write unique hits file (deduplicated by keeping highest score)
    unique_path = out_path.replace('.tsv', '_unique.tsv').replace('.txt', '_unique.txt')
    if unique_path == out_path:  # fallback if no extension
        unique_path = out_path + '_unique'

    with open(unique_path, "w", encoding="utf-8") as unique_out:
        unique_out.write("\t".join(["motif", "gene_name", "ensembl_id", "start", "end", "strand", "score", "matched_sequence", "reverse"]) + "\n")

        # Sort by motif, gene_name, start position for consistent output
        sorted_hits = sorted(all_hits.items(), key=lambda x: (x[0][0], x[0][1], x[0][2]))

        for (motif, gene_name, ensembl_id, start, end, strand, matched_seq, reverse_seq), score in sorted_hits:
            unique_out.write("\t".join([motif, gene_name, ensembl_id, str(start), str(end), strand, str(score), matched_seq, reverse_seq]) + "\n")

    print(f"[INFO] Wrote {len(all_hits)} unique motif matches to {unique_path}")

    if warnings > 10:
        print(f"[INFO] Total warnings: {warnings} (only first 10 shown)", file=sys.stderr)

if __name__ == "__main__":
    main()