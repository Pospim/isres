#!/bin/bash

# Script to remove duplicate rows from a TSV file while preserving the header

if [ $# -ne 2 ]; then
    echo "Usage: $0 <input.tsv> <output.tsv>"
    echo "Example: $0 input.tsv output_unique.tsv"
    exit 1
fi

input_file="$1"
output_file="$2"

if [ ! -f "$input_file" ]; then
    echo "Error: Input file '$input_file' not found"
    exit 1
fi

# Extract header (first line)
head -n 1 "$input_file" > "$output_file"

# Remove duplicates from the rest of the file and append to output
# sort -u removes duplicate lines, keeping only unique rows
tail -n +2 "$input_file" | sort -u >> "$output_file"

echo "Duplicates removed. Output saved to: $output_file"
echo "Original lines: $(wc -l < "$input_file")"
echo "Unique lines: $(wc -l < "$output_file")"
