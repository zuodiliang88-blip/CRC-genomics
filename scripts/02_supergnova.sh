#!/usr/bin/env bash
set -euo pipefail
# Usage: 02_supergnova.sh TRAIT N_TRAIT
trait=${1:?trait}; n=${2:?sample size}
: "${CRC_DATA_DIR:?}" "${SUPERGNOVA_PY:?}" "${SUPERGNOVA_BFILE:?}" "${SUPERGNOVA_PARTITION:?}"
mkdir -p "$CRC_DATA_DIR/results/supergnova"
python "$SUPERGNOVA_PY" \
  "$CRC_DATA_DIR/ldsc/$trait.sumstats.gz" "$CRC_DATA_DIR/ldsc/CRC.sumstats.gz" \
  --N1 "$n" --N2 185616 \
  --bfile "$SUPERGNOVA_BFILE" --partition "$SUPERGNOVA_PARTITION" \
  --out "$CRC_DATA_DIR/results/supergnova/$trait.CRC.2353.txt" --thread 32
