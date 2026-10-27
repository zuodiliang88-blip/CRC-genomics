#!/usr/bin/env bash
set -euo pipefail

# Main MiXeR analysis: 20 fit/test replicates for CRC and each of 13 factors.
: "${MIXER_PYTHON:?}" "${MIXER_SCRIPT:?}" "${MIXER_LIB:?}"
: "${MIXER_REF_PREFIX:?}" "${MIXER_SUMSTATS_DIR:?}" "${MIXER_OUT_DIR:?}"

crc=CRC_pmid36539618
traits=(
  BMI_GCST009004
  Body_fat_percentage_GCST90468160
  WHR_GCST008996
  Smoking_initiation_GSCAN_SmkInit_2022_EUR
  Smoking_cessation_GSCAN_SmkCes_2022_EUR
  Cigarettes_per_day_GSCAN_CigDay_2022_EUR
  Educational_attainment_pmid35361970
  ISEI_occupational_score_GCST90446160
  Household_income_GCST90566699_pmid40360725
  MVPA_GCST90104341
  Age_at_smoking_initiation_GSCAN_AgeSmk_2022_EUR
  Red_wine_GCST90132987
  White_wine_GCST90132986
)
threads=${MIXER_THREADS:-40}
ref=$MIXER_REF_PREFIX
common=(--bim-file "${ref}.@.bim" --ld-file "${ref}.@.run4.ld"
        --lib "$MIXER_LIB" --threads "$threads")

for trait in "$crc" "${traits[@]}"; do
  input="$MIXER_SUMSTATS_DIR/${trait}_noMHC.sumstats.gz"
  out="$MIXER_OUT_DIR/univariate/$trait"
  mkdir -p "$out"
  for i in $(seq 1 20); do
    fit="$out/${trait}_noMHC.fit.rep$i"
    test="$out/${trait}_noMHC.test.rep$i"
    "$MIXER_PYTHON" "$MIXER_SCRIPT" fit1 \
      --trait1-file "$input" --out "$fit" \
      --extract "${ref}.prune_maf0p05_rand2M_r2p8.rep$i.snps" "${common[@]}"
    "$MIXER_PYTHON" "$MIXER_SCRIPT" test1 \
      --trait1-file "$input" --load-params-file "$fit.json" \
      --out "$test" "${common[@]}"
  done
done

crc_input="$MIXER_SUMSTATS_DIR/${crc}_noMHC.sumstats.gz"
for trait in "${traits[@]}"; do
  input="$MIXER_SUMSTATS_DIR/${trait}_noMHC.sumstats.gz"
  out="$MIXER_OUT_DIR/bivariate/${crc}__${trait}"
  mkdir -p "$out"
  for i in $(seq 1 20); do
    fit="$out/${crc}__${trait}.fit.rep$i"
    test="$out/${crc}__${trait}.test.rep$i"
    "$MIXER_PYTHON" "$MIXER_SCRIPT" fit2 \
      --trait1-file "$crc_input" --trait2-file "$input" \
      --trait1-params-file "$MIXER_OUT_DIR/univariate/$crc/${crc}_noMHC.fit.rep$i.json" \
      --trait2-params-file "$MIXER_OUT_DIR/univariate/$trait/${trait}_noMHC.fit.rep$i.json" \
      --out "$fit" --extract "${ref}.prune_maf0p05_rand2M_r2p8.rep$i.snps" \
      "${common[@]}"
    "$MIXER_PYTHON" "$MIXER_SCRIPT" test2 \
      --trait1-file "$crc_input" --trait2-file "$input" \
      --load-params-file "$fit.json" --out "$test" "${common[@]}"
  done
done
