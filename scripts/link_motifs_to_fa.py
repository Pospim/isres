#!/usr/bin/env python3
import sys,re,html, os
from collections import defaultdict

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

def main():
    if len(sys.argv) != 4:
        print("Usage: python convert_homer_offset_hits.py <fasta.fa> <homer_offset_hits.txt> <out_hits_dir>", file=sys.stderr)
        sys.exit(1)

    fa, hits_txt, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
    os.makedirs(outdir, exist_ok=True)

    seqs = parse_fasta(fa)
    if not seqs:
        sys.exit(1)
    out_path = os.path.join(outdir, "homerMotifs.all.txt")

    out = open(out_path, "w", encoding="utf-8")

    with open(hits_txt, encoding="utf-8", errors="ignore") as f:
        header=f.readline()

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
                    print(f"[WARN] Could not locate match in FASTA for {fid} offset {offset} motif {motif}", file=sys.stderr)
                    continue
            start1 = start0 + 1
            end1 = start1 + motif_len - 1
            out.write("\t".join([motif, f">{fid}", str(start1), str(end1), strand, score, motif_seq]) + "\n")

        out.close()

if __name__ == "__main__":
    main()