#!/bin/bash

# Script to find motif hits in FASTA sequences using HOMER
# Usage: find_hits.sh --fasta FILE --motifs FILE [OPTIONS]

set -e

# Default values
FASTA_FILE=""
MOTIF_FILE=""
OUTPUT_FILE="motif_hits.txt"
WORK_DIR="."

# Function to show usage
usage() {
    cat << EOF
Usage: $0 --fasta FILE --motifs FILE [OPTIONS]

Required arguments:
  -f, --fasta FILE          Input FASTA file with sequences
  -m, --motifs FILE         HOMER motifs file to search for

Optional arguments:
  -w, --work-dir DIR        Working directory to cd into before running (default: current)
  -h, --help            Show this help message

Example:
  $0 --fasta isg_2k_upstream.fa --motifs homerMotifs.all.motifs --work-dir /path/to/dir
EOF
    exit 1
}

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -f|--fasta)
            FASTA_FILE="$2"
            shift 2
            ;;
        -m|--motifs)
            MOTIF_FILE="$2"
            shift 2
            ;;
        -w|--work-dir)
            WORK_DIR="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Error: Unknown argument: $1"
            usage
            ;;
    esac
done

# Check required arguments
if [[ -z "$FASTA_FILE" ]]; then
    echo "Error: --fasta is required"
    usage
fi

if [[ -z "$MOTIF_FILE" ]]; then
    echo "Error: --motifs is required"
    usage
fi

OUTPUT_FILE="homerMotifs.all.motifs_hits.txt"

# Convert all paths to absolute
FASTA_FILE=$(realpath "$FASTA_FILE")
MOTIF_FILE=$(realpath "$MOTIF_FILE")
WORK_DIR=$(realpath "$WORK_DIR")

# Set output directory to motif_hits within working directory
OUTPUT_DIR="$WORK_DIR/motif_hits"

echo "[INFO] Working directory: $WORK_DIR"
echo "[INFO] Input FASTA: $FASTA_FILE"
echo "[INFO] Motif file: $MOTIF_FILE"
echo "[INFO] Output directory: $OUTPUT_DIR"
echo "[INFO] Output file: $OUTPUT_FILE"

# Check if input files exist
if [[ ! -f "$FASTA_FILE" ]]; then
    echo "Error: FASTA file not found: $FASTA_FILE"
    exit 1
fi

if [[ ! -f "$MOTIF_FILE" ]]; then
    echo "Error: Motif file not found: $MOTIF_FILE"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run HOMER findMotifs.pl
echo "[INFO] Running HOMER findMotifs.pl..."
findMotifs.pl "$FASTA_FILE" fasta "$OUTPUT_DIR" \
  -find "$MOTIF_FILE" \
  > "$OUTPUT_DIR/$OUTPUT_FILE"

echo "[INFO] Results written to: $OUTPUT_DIR/$OUTPUT_FILE"
