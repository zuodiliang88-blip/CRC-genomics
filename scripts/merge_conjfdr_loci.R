library(data.table)
source("config/config.R")
fdr_dir <- file.path(data_dir, "conjfdr")
input_file <- file.path(fdr_dir, "All_13_phenotypes_CRC_conjfdr_0.05_fuma.loci.csv")
dir.create(file.path(out_root, "conjfdr"), recursive = TRUE, showWarnings = FALSE)

# 1. 读取 13 个表型与 CRC 的 conjFDR loci 汇总表。
#    显著性阈值为 conjFDR < 0.05。
loci <- fread(input_file)

# 2. 检查必要列，并确认确实包含 13 个表型。
stopifnot(
  all(c("PHENOTYPE", "locusnum", "CHR", "LEAD_SNP", "LEAD_BP", "MinBP", "MaxBP", "FDR") %in% names(loci)),
  uniqueN(loci$PHENOTYPE) == 13L
)

# 3. 保留有染色体和区间坐标的 loci，并把坐标和 FDR 转为数值。
loci <- loci[!is.na(CHR) & !is.na(MinBP) & !is.na(MaxBP)]
loci[, CHR := as.integer(CHR)]
loci[, MinBP := as.integer(MinBP)]
loci[, MaxBP := as.integer(MaxBP)]
loci[, LEAD_BP := as.integer(LEAD_BP)]
loci[, FDR := as.numeric(FDR)]

# 4. 按染色体和基因组位置排序。
setorder(loci, CHR, MinBP, MaxBP, FDR)

# 5. 给每个原始 locus 分配跨表型合并后的 region 编号。
#    同一染色体上，如果下一个 locus 的起点距离当前 region 终点
#    不超过 250 kb，就把它们归入同一个 region。
loci[, region_number := NA_integer_]
region_number <- 0L

for (chromosome in sort(unique(loci$CHR))) {
  row_index <- which(loci$CHR == chromosome)
  merged_region_end <- -Inf

  for (i in row_index) {
    starts_new_region <- is.infinite(merged_region_end) ||
      loci$MinBP[i] > merged_region_end + 250000

    if (starts_new_region) {
      region_number <- region_number + 1L
      merged_region_end <- loci$MaxBP[i]
    } else {
      merged_region_end <- max(merged_region_end, loci$MaxBP[i])
    }

    loci$region_number[i] <- region_number
  }
}

# 6. 每个 region 选择 FDR 最小的原始 locus 作为代表 locus。
representative_loci <- loci[order(FDR), .SD[1], by = region_number]

# 7. 汇总每个 region 的边界、原始 loci 数量、表型和原始 lead SNP。
merged_regions <- loci[, .(
  CHR = CHR[1],
  REGION_START = min(MinBP),
  REGION_END = max(MaxBP),
  N_ORIGINAL_LOCI = .N,
  N_PHENOTYPES = uniqueN(PHENOTYPE),
  PHENOTYPES = paste(sort(unique(PHENOTYPE)), collapse = ";"),
  ORIGINAL_LEAD_SNPS = paste(unique(LEAD_SNP), collapse = ";"),
  ORIGINAL_LOCUS_IDS = paste(paste(PHENOTYPE, locusnum, sep = ":"), collapse = ";")
), by = region_number]

# 8. 加入代表 locus 的表型、locus 编号、lead SNP、位置和 FDR。
merged_regions <- merge(
  merged_regions,
  representative_loci[, .(
    region_number,
    REPRESENTATIVE_PHENOTYPE = PHENOTYPE,
    REPRESENTATIVE_LOCUSNUM = locusnum,
    REPRESENTATIVE_LEAD_SNP = LEAD_SNP,
    REPRESENTATIVE_LEAD_BP = LEAD_BP,
    REPRESENTATIVE_FDR = FDR
  )],
  by = "region_number"
)

# 9. 计算 region 宽度，并区分表型特异与表型共享 region。
#    只出现在 1 个表型中：phenotype_specific。
#    出现在至少 2 个表型中：phenotype_shared。
merged_regions[, REGION_WIDTH := REGION_END - REGION_START + 1L]
merged_regions[, LOCUS_CLASS := ifelse(
  N_PHENOTYPES == 1,
  "phenotype_specific",
  "phenotype_shared"
)]

# 10. 按基因组位置排序并生成正式 region ID。
setorder(merged_regions, CHR, REGION_START, REGION_END)
merged_regions[, merged_region_id := sprintf(
  "conjFDR_region_%03d",
  seq_len(.N)
)]

# 11. 把正式 region ID 回填到每一个原始 locus。
original_loci_with_region <- merge(
  loci,
  merged_regions[, .(
    region_number,
    merged_region_id,
    REGION_START,
    REGION_END
  )],
  by = "region_number",
  all.x = TRUE
)

setcolorder(
  original_loci_with_region,
  c(
    "merged_region_id",
    "REGION_START",
    "REGION_END",
    setdiff(
      names(original_loci_with_region),
      c("merged_region_id", "REGION_START", "REGION_END")
    )
  )
)
setorder(original_loci_with_region, CHR, REGION_START, MinBP, FDR)

# 12. 最终 region 表不保留内部使用的数字编号。
merged_regions[, region_number := NULL]
setcolorder(
  merged_regions,
  c("merged_region_id", setdiff(names(merged_regions), "merged_region_id"))
)

# 13. 输出最终非冗余 region 表。
fwrite(
  merged_regions,
  file.path(out_root, "conjfdr", "merged_regions.csv")
)

# 14. 输出每个原始 locus 与最终 region 的对应表。
fwrite(
  original_loci_with_region,
  file.path(out_root, "conjfdr", "original_loci_with_region.csv")
)

# 15. 打印检查结果。
cat("原始 loci 数量：", nrow(original_loci_with_region), "\n")
cat("合并后 regions 数量：", nrow(merged_regions), "\n")
cat("表型特异 regions：", merged_regions[LOCUS_CLASS == "phenotype_specific", .N], "\n")
cat("表型共享 regions：", merged_regions[LOCUS_CLASS == "phenotype_shared", .N], "\n")
print(merged_regions[, .N, by = N_PHENOTYPES][order(N_PHENOTYPES)])
