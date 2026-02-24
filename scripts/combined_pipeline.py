#!/usr/bin/env python3
"""
Pipeline script that runs the full motif analysis:
  1. find_hits.sh            — find motif hits in FASTA using HOMER
  2. link_motifs_to_fasta.py — convert HOMER hits to coordinate format
  3. map_motifs_to_html.py   — generate interactive HTML visualization
"""

import argparse
import subprocess
import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIND_HITS_SCRIPT = "find_hits.sh"
LINK_MOTIFS_SCRIPT = "link_motifs_to_fasta.py"
MAP_TO_HTML_SCRIPT = "map_motifs_to_html.py"

def run_cmd(cmd, description):
    """Run a shell command and exit on failure."""
    print(f"\n{'='*60}")
    print(f"STEP: {description}")
    print(f"CMD:  {' '.join(cmd)}")
    print(f"{'='*60}\n")

    result = subprocess.run(cmd, cwd=SCRIPT_DIR)
    if result.returncode != 0:
        print(f"\n[ERROR] Step failed: {description}", file=sys.stderr)
        sys.exit(result.returncode)

def main():
    parser = argparse.ArgumentParser(
        description="Run the full motif analysis pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Pipeline steps:
  1. find_hits.sh
     Runs HOMER findMotifs.pl to scan FASTA for motif hits.
     Input:  --fasta, --motifs
     Output: <outdir>/motif_hits/homerMotifs.all.motifs_hits.txt

  2. link_motifs_to_fasta.py
     Converts HOMER offset-based hits into coordinate TSV.
     Input:  --fasta, HOMER hits from step 1
     Output: <outdir>/homerMotifs.all_coords.tsv
             <outdir>/homerMotifs.all_coords_unique.tsv

  3. map_motifs_to_html.py
     Generates interactive HTML + heatmap from coordinate TSV.
     Input:  --fasta, coordinate TSV from step 2
     Output: <outdir>/motif_map.html
             <outdir>/motif_map_heatmap_topN.html

Example usage:
  python combined_pipeline.py \\
      --fasta data/sequences.fa \\
      --motifs data/homerMotifs.all.motifs \\
      --outdir results/ \\
      --format gene_ensembl \\
      --top-motifs 10
        """,
    )

    parser.add_argument(
        "-f", "--fasta", required=True, help="FASTA file with upstream sequences"
    )
    parser.add_argument(
        "-m", "--motifs", required=True, help="HOMER motifs file to search for"
    )
    parser.add_argument(
        "-o", "--outdir", required=True, help="Output directory for all results"
    )
    parser.add_argument(
        "--format",
        choices=["ensembl", "gene_ensembl", "gene"],
        default="gene_ensembl",
        help="FASTA header format (default: gene_ensembl)",
    )
    parser.add_argument(
        "-n",
        "--top-motifs",
        type=int,
        default=5,
        help="Number of top motifs for heatmap (default: 5)",
    )
    parser.add_argument(
        "--html-name",
        default="motif_map.html",
        help="Output HTML filename (default: motif_map.html)",
    )
    parser.add_argument(
        "--use-unique",
        action="store_true",
        default=False,
        help="Use deduplicated hits for visualization (default: False)",
    )
    parser.add_argument(
        "--skip-homer",
        action="store_true",
        default=False,
        help="Skip step 1 (find_hits.sh) if HOMER hits already exist",
    )

    args = parser.parse_args()
    # Resolve all paths to absolute (relative to user's CWD, not SCRIPT_DIR)
    args.fasta = os.path.abspath(args.fasta)
    args.motifs = os.path.abspath(args.motifs)
    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    homer_hits = os.path.join(outdir, "motif_hits", "homerMotifs.all.motifs_hits.txt")

    # ── Step 1: find_hits.sh ─────────────────────────────────────
    # Input:  --fasta FILE, --motifs FILE, --work-dir DIR
    # Output: <work-dir>/motif_hits/homerMotifs.all.motifs_hits.txt
    #   Uses HOMER findMotifs.pl internally with -find flag
    if args.skip_homer and os.path.isfile(homer_hits):
        print(f"[INFO] Skipping HOMER (--skip-homer), using existing: {homer_hits}")
    else:
        find_hits = os.path.join(SCRIPT_DIR, FIND_HITS_SCRIPT)
        run_cmd(
            [
                "bash",
                find_hits,
                "--fasta", args.fasta,
                "--motifs", args.motifs,
                "--work-dir", outdir,
            ],
            "Find motif hits using HOMER (findMotifs.pl)"
        )
    if not os.path.isfile(homer_hits):
        print(f"[ERROR] Expected HOMER hits not found: {homer_hits}", file=sys.stderr)
        sys.exit(1)
    print(f"[OK] HOMER hits: {homer_hits}")

    # ── Step 2: link_motifs_to_fasta.py ──────────────────────────
    # Input:  -f/--fasta FILE, -i/--hits FILE, -o/--output-dir DIR, --output-file NAME
    #   --hits: the HOMER hits file from step 1
    # Output: <output-dir>/homerMotifs.all_coords.tsv       (all hits)
    #         <output-dir>/homerMotifs.all_coords_unique.tsv (deduplicated, highest score kept)
    #   TSV columns: motif | gene_name | ensembl_id | start | end | strand | score | matched_sequence | reverse
    coords_file = "homerMotifs.all_coords.tsv"
    link_script = os.path.join(SCRIPT_DIR, LINK_MOTIFS_SCRIPT)
    run_cmd(
        [
            sys.executable,
            link_script,
            "--fasta", args.fasta,
            "--hits", homer_hits,
            "--output-dir", outdir,
            "--output-file", coords_file,
        ],
        "Convert HOMER hits to coordinate formate."
    )
    coords_tsv_all = os.path.join(outdir, coords_file)
    coords_tsv_unique = os.path.join(outdir, coords_file.replace(".tsv", "_unique.tsv"))

    if args.use_unique and os.path.isfile(coords_tsv_unique):
        coords_tsv = coords_tsv_unique
        print(f"[OK] Using deduplicated hits: {coords_tsv}")
    elif os.path.isfile(coords_tsv_all):
        coords_tsv = coords_tsv_all
        print(f"[OK] Using all hits: {coords_tsv}")
    else:
        print(
            f"[ERROR] No coordinate TSV found. "
            f"Checked:\n  {coords_tsv_unique}\n  {coords_tsv_all}",
            file=sys.stderr,
        )
        sys.exit(1)

    # ── Step 3: map_motifs_to_html.py ────────────────────────────
    # Input:  -t/--tsv FILE, -f/--fasta FILE, -o/--outfile FILE, --format STR, -n/--top_motifs INT
    #   --tsv: coordinate TSV from step 2
    #   --fasta: same FASTA as step 1 & 2
    #   --format: ensembl | gene_ensembl | gene
    #   --top_motifs: number of top motifs for heatmap
    # Output: <outfile>                              (interactive HTML with sequence view)
    #         <outfile_stem>_heatmap_top<N>.html      (plotly heatmap in same directory)
    html_out = os.path.join(outdir, args.html_name)
    map_script = os.path.join(SCRIPT_DIR, MAP_TO_HTML_SCRIPT)
    run_cmd(
        [
            sys.executable,
            map_script,
            "--tsv", coords_tsv,
            "--fasta", args.fasta,
            "--outfile", html_out,
            "--format", args.format,
            "--top_motifs", str(args.top_motifs),
        ],
        "Generate interactive HTML visualization (map_motifs_to_html.py)"
    )

    # ── Summary ──────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print("Pipeline finished successfully!")
    print(f"  HOMER hits       : {homer_hits}")
    print(f"  Coords (all)     : {coords_tsv_all}")
    print(f"  Coords (unique)  : {coords_tsv_unique}")
    print(f"  Used for HTML    : {coords_tsv}")
    print(f"  Interactive HTML : {html_out}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
