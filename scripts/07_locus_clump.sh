#!/usr/bin/env bash
set -euo pipefail
# Usage: 07_locus_clump.sh FDR_RESULTS.tsv OUTPUT_PREFIX SIGNIFICANCE
# SIGNIFICANCE is 0.05 for conjFDR or 0.01 for condFDR.
input=${1:?input file}; output=${2:?output prefix}; threshold=${3:?threshold}
: "${FDR_CLUMP_PY:?}" "${FDR_REF_BFILE_CHR:?}" "${PLINK_BIN:?}"
python "$FDR_CLUMP_PY" clump --sumstats "$input" --out "$output" \
  --clump-field FDR --clump-snp-field SNP --chr CHR \
  --clump-p1 "$threshold" --candidate-fdr 0.1 \
  --indep-r2 0.6 --lead-r2 0.1 --ld-window-kb 10000 \
  --loci-merge-kb 250 --bfile-chr "$FDR_REF_BFILE_CHR" \
  --plink "$PLINK_BIN" --force
