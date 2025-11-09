#!/usr/bin/env bash
# make_negatives.sh
# Build GC-matched negative promoters for HOMER from a list of human genes.
# Requires: HOMER (annotatePeaks.pl, pos2bed.pl), coreutils (awk, sort, shuf), optional: bedtools
set -euo pipefail

usage() {
  cat <<'USAGE'
USAGE:
  make_negatives.sh -g hg38 -p -2000,0 -i genes_pos.txt -o out_dir [-t 8] [--no-bedtools] [--skip-find]

Inputs
  -g   HOMER genome key (e.g., hg38)
  -p   Promoter window (e.g., -2000,0 for 2 kb upstream)
  -i   genes_pos.txt (one gene per line; symbols or Ensembl IDs)
  -o   Output directory (created if missing)

Options
  -t   Threads for HOMER (default: 8)
  --no-bedtools   Do not subtract overlaps with bedtools
  --skip-find     Build negatives only; skip findMotifsGenome.pl

Outputs (in -o)
  all_promoters.bed
  pos_promoters.bed
  pool_neg_candidates.bed
  pool_neg_nooverlap.bed (if bedtools)
  pos.gc.tsv / pool.gc.tsv
  neg_promoters.bed
  out_denovo/ (if HOMER run was not skipped)

Example
  ./make_negatives.sh -g hg38 -p -2000,0 -i genes_pos.txt -o human/out_2kb_upstream -t 12
USAGE
}

# ---- parse args ----
GENOME=""; PROM=""; POSLIST=""; OUTDIR=""
THREADS=8; USE_BEDTOOLS=1; DO_FIND=1
while [[ $# -gt 0 ]]; do
  case "$1" in
    -g) GENOME="$2"; shift 2 ;;
    -p) PROM="$2"; shift 2 ;;
    -i) POSLIST="$2"; shift 2 ;;
    -o) OUTDIR="$2"; shift 2 ;;
    -t) THREADS="$2"; shift 2 ;;
    --no-bedtools) USE_BEDTOOLS=0; shift ;;
    --skip-find) DO_FIND=0; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1"; usage; exit 1 ;;
  esac
done

[[ -z "$GENOME" || -z "$PROM" || -z "$POSLIST" || -z "$OUTDIR" ]] && { usage; exit 1; }

# ---- preflight checks ----
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found in PATH."; exit 1; }; }
need annotatePeaks.pl
need pos2bed.pl
need awk
need sort
need shuf

if [[ $USE_BEDTOOLS -eq 1 ]]; then
  if ! command -v bedtools >/dev/null 2>&1; then
    echo "WARN: bedtools not found; proceeding without overlap subtraction (--no-bedtools)."
    USE_BEDTOOLS=0
  fi
fi

# ---- robust abspath (works on Linux/WSL) ----
abspath() {
  if command -v realpath >/dev/null 2>&1; then
    realpath -m "$1"
  elif command -v python3 >/dev/null 2>&1; then
    python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$1"
  else
    # POSIX fallback
    echo "$(cd "$(dirname -- "$1")" && pwd -P)/$(basename -- "$1")"
  fi
}

# ---- resolve absolute paths & prepare outdir ----
POSLIST_ABS="$(abspath "$POSLIST")"
OUTDIR_ABS="$(abspath "$OUTDIR")"

mkdir -p "$OUTDIR_ABS"
cd "$OUTDIR_ABS"

# ---- step 1: build promoter universe once, then positives by Gene Name ----
echo "[1/7] Building promoter universe for ${GENOME} with window ${PROM} ..."
annotatePeaks.pl tss "$GENOME" -size "$PROM" > all_tss.annot

# Convert full annotation to BED = the universe
pos2bed.pl all_tss.annot > all_promoters.bed
N_ALL=$(wc -l < all_promoters.bed | tr -d ' ')
echo "    total promoters (universe): ${N_ALL}"

echo "[1b/7] Selecting positives by Gene Name from ${POSLIST_ABS} ..."
# Keep only rows whose "Gene Name" matches your positives (symbols)
awk -v FS='\t' -v OFS='\t' '
  NR==FNR { pos[$1]=1; next }          # load genes_pos.txt (symbols)
  NR==1  {                              # find "Gene Name" column
    for (i=1;i<=NF;i++)
      if ($i=="Gene Name" || $i=="GeneName" || $i=="Gene_Name") gn=i
    if (!gn) { print "ERROR: Gene Name column not found in annotatePeaks output." > "/dev/stderr"; exit 1 }
    print; next
  }
  ($gn in pos)                          # keep rows where Gene Name ∈ positives
' "$POSLIST_ABS" all_tss.annot > pos_tss.annot

# Convert that subset to BED = your positive promoters
pos2bed.pl pos_tss.annot > pos_promoters.bed
N_POS=$(wc -l < pos_promoters.bed | tr -d ' ')
echo "    positives: ${N_POS} promoters"
if [[ "$N_POS" -eq 0 ]]; then
  echo "ERROR: No promoters matched your gene list (by Gene Name). Check symbols/build."; exit 1
fi

# ---- step 2: exclude positives by Gene Name in the annotation, then to BED ----
echo "[2/7] Excluding positives (by Gene Name) to form candidate pool ..."
# Build a set of positive symbols
awk '{print $1}' "$POSLIST_ABS" | sort -u > pos_symbols.set

# Keep rows whose Gene Name is NOT in positives
awk -v FS='\t' -v OFS='\t' '
  NR==FNR { pos[$1]=1; next }                 # load pos gene symbols
  NR==1  {
    for (i=1;i<=NF;i++)
      if ($i=="Gene Name" || $i=="GeneName" || $i=="Gene_Name") gn=i
    if (!gn) { print "ERROR: Gene Name column not found in annotatePeaks output." > "/dev/stderr"; exit 1 }
    print; next
  }
  !($gn in pos)                               # keep rows where Gene Name ∉ positives
' pos_symbols.set all_tss.annot > pool_neg_candidates.annot

# Convert pool annotation to BED
pos2bed.pl pool_neg_candidates.annot > pool_neg_candidates.bed

N_POOL=$(wc -l < pool_neg_candidates.bed | tr -d ' ')
echo "    candidate negatives (by gene): ${N_POOL}"

# ---- step 3: optional overlap subtraction ----
NEGPOOL="pool_neg_candidates.bed"
if [[ $USE_BEDTOOLS -eq 1 ]]; then
  echo "[3/7] Removing overlaps with positives using bedtools ..."
  BEFORE="$N_POOL"
  bedtools subtract -a pool_neg_candidates.bed -b pos_promoters.bed > pool_neg_nooverlap.bed
  NEGPOOL="pool_neg_nooverlap.bed"
  N_POOL2=$(wc -l < "$NEGPOOL" | tr -d ' ')
  REMOVED=$(( BEFORE - N_POOL2 ))
  echo "    non-overlapping pool: ${N_POOL2}  (removed overlaps: ${REMOVED})"
else
  echo "[3/7] Skipping overlap subtraction (no bedtools)."
fi

# ---- step 4: GC% annotate both sets ----
echo "[4/7] Computing GC% ..."
annotatePeaks.pl pos_promoters.bed "$GENOME" -size given -noann -nogene -GC > pos.gc.tsv
annotatePeaks.pl "$NEGPOOL"       "$GENOME" -size given -noann -nogene -GC > pool.gc.tsv

# ---- step 5: GC-binned sampling to match positives 1:1 ----
echo "[5/7] GC-binned sampling of negatives to match positives ..."
N="$N_POS"

# extract GC values
awk 'NR>1{gc=$(NF); print gc}' pos.gc.tsv | sort -n > pos.gc.values
awk 'NR>1{gc=$(NF); print gc"\t"$0}' pool.gc.tsv > pool.gc.withline

# choose number of bins
K=5
# quantile cutpoints from positives
awk -v k="$K" '{
  a[NR]=$1
} END {
  if (NR==0) exit 0;
  for (i=1;i<=k-1;i++){
    p=i/k; idx=int(p*NR); if(idx<1) idx=1; if(idx>NR) idx=NR;
    print a[idx];
  }
}' pos.gc.values > gc.cuts

# count positives per bin
awk 'NR==FNR{cut[++m]=$1; next}
     {gc=$1; b=1; for(i=1;i<=m;i++) if(gc>cut[i]) b=i+1; cnt[b]++}
     END{for(i=1;i<=m+1;i++) printf "%d\t%d\n", i, (cnt[i]?cnt[i]:0)}' gc.cuts pos.gc.values > pos.bin.counts

# tag pool lines with bin id
awk 'NR==FNR{cut[++m]=$1; next}
     {gc=$1; b=1; for(i=1;i<=m;i++) if(gc>cut[i]) b=i+1; print b"\t"$0}' gc.cuts pool.gc.withline \
| sort -k1,1n > pool.binned

# sample per bin
: > neg_promoters.bed
while read -r BIN COUNT; do
  [[ "$COUNT" -eq 0 ]] && continue
  awk -v bin="$BIN" '$1==bin{print $0}' pool.binned \
    | shuf | head -n "$COUNT" \
    | cut -f3- \
    | awk -v OFS="\t" '{print $2,$3,$4,$5,$6,$7}' >> neg_promoters.bed
done < pos.bin.counts

# top-up if short (rare)
LINES=$(wc -l < neg_promoters.bed | tr -d ' ')
if (( LINES < N )); then
  echo "    Top-up negatives to reach ${N} (short by $((N-LINES)))."
  shuf "$NEGPOOL" | head -n $((N-LINES)) >> neg_promoters.bed
fi

# ---- step 6: run HOMER (optional) ----
if [[ $DO_FIND -eq 1 ]]; then
  echo "[6/7] Running HOMER de novo motif discovery with custom background ..."
  findMotifsGenome.pl pos_promoters.bed "$GENOME" out_denovo \
    -bg neg_promoters.bed -size given -mask -len 8,10,12 -p "$THREADS"
  echo "DONE. Results in: $(pwd)/out_denovo"
else
  echo "[6/7] Skipped HOMER run (--skip-find). Negatives ready."
fi

# ---- summary ----
echo "Summary:"
echo "  positives:        ${N_POS}"
echo "  negatives (final): $(wc -l < neg_promoters.bed | tr -d ' ')"
echo "  outputs in:       $(pwd)"
