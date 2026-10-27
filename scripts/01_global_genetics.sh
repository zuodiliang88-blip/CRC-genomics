#!/usr/bin/env bash
set -euo pipefail
# Usage: 01_global_genetics.sh munge|h2|rg|gnova TRAIT [N]
# Inputs are harmonized GWAS files under $CRC_DATA_DIR/gwas and LDSC sumstats under $CRC_DATA_DIR/ldsc.
mode=${1:?mode}; trait=${2:?trait}; pcol=${3:-pval}
: "${CRC_DATA_DIR:?}" "${LDSC_DIR:?}" "${LDSC_LD:?}"
mkdir -p "$CRC_DATA_DIR/ldsc" "$CRC_DATA_DIR/results/global"
case "$mode" in
  munge)
    : "${HM3_SNPLIST:?}"
    python "$LDSC_DIR/munge_sumstats.py" \
      --sumstats "$CRC_DATA_DIR/gwas/$trait.txt.gz" \
      --snp SNP --a1 A1 --a2 A2 --signed-sumstats beta,0 \
      --p "$pcol" --N-col n --merge-alleles "$HM3_SNPLIST" \
      --chunksize 100000 --out "$CRC_DATA_DIR/ldsc/$trait"
    ;;
  h2)
    python "$LDSC_DIR/ldsc.py" --h2 "$CRC_DATA_DIR/ldsc/$trait.sumstats.gz" \
      --ref-ld-chr "$LDSC_LD" --w-ld-chr "$LDSC_LD" \
      --out "$CRC_DATA_DIR/results/global/$trait.h2"
    ;;
  rg)
    python "$LDSC_DIR/ldsc.py" \
      --rg "$CRC_DATA_DIR/ldsc/$trait.sumstats.gz,$CRC_DATA_DIR/ldsc/CRC.sumstats.gz" \
      --ref-ld-chr "$LDSC_LD" --w-ld-chr "$LDSC_LD" \
      --out "$CRC_DATA_DIR/results/global/$trait.CRC.rg"
    ;;
  gnova)
    : "${GNOVA_PY:?}"
    python "$GNOVA_PY" \
      "$CRC_DATA_DIR/gnova/$trait.clean.sumstats.gz" \
      "$CRC_DATA_DIR/gnova/CRC.clean.sumstats.gz" \
      --bfile "${GNOVA_BFILE:?}" \
      --out "$CRC_DATA_DIR/results/global/$trait.CRC.gnova.txt"
    ;;
  *) echo "Unknown mode: $mode" >&2; exit 2;;
esac
