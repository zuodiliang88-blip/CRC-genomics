#!/usr/bin/env bash
set -euo pipefail
# Usage: 03_twas.sh TRAIT
trait=${1:?trait}
: "${CRC_DATA_DIR:?}" "${FUSION_R:?}" "${FUSION_WEIGHTS:?}"
: "${FUSION_WEIGHTS_DIR:?}" "${FUSION_LD:?}" "${FUSION_PSEUDOVAL:?}"
mkdir -p "$CRC_DATA_DIR/results/twas/$trait"
for chr in $(seq 1 22); do
  Rscript "$FUSION_R" \
    --sumstats "$CRC_DATA_DIR/twas_sumstats/$trait.sumstats.gz" \
    --weights "$FUSION_WEIGHTS" --weights_dir "$FUSION_WEIGHTS_DIR" \
    --ref_ld_chr "$FUSION_LD" --pseudovalidation "$FUSION_PSEUDOVAL" \
    --chr "$chr" --out "$CRC_DATA_DIR/results/twas/$trait/$trait.chr$chr.dat"
done
