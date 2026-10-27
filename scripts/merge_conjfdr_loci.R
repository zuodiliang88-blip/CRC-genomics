library(data.table)
source("config/config.R")
fdr_dir <- file.path(data_dir, "conjfdr")
input_file <- file.path(fdr_dir, "All_13_phenotypes_CRC_conjfdr_0.05_fuma.loci.csv")
dir.create(file.path(out_root, "conjfdr"), recursive = TRUE, showWarnings = FALSE)

loci <- fread(input_file)

stopifnot(
  all(c("PHENOTYPE", "locusnum", "CHR", "LEAD_SNP", "LEAD_BP", "MinBP", "MaxBP", "FDR") %in% names(loci)),
  uniqueN(loci$PHENOTYPE) == 13L
)

loci <- loci[!is.na(CHR) & !is.na(MinBP) & !is.na(MaxBP)]
loci[, CHR := as.integer(CHR)]
loci[, MinBP := as.integer(MinBP)]
loci[, MaxBP := as.integer(MaxBP)]
loci[, LEAD_BP := as.integer(LEAD_BP)]
loci[, FDR := as.numeric(FDR)]

setorder(loci, CHR, MinBP, MaxBP, FDR)

# Merge loci on the same chromosome when separated by at most 250 kb.
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

representative_loci <- loci[order(FDR), .SD[1], by = region_number]

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

merged_regions[, REGION_WIDTH := REGION_END - REGION_START + 1L]
merged_regions[, LOCUS_CLASS := ifelse(
  N_PHENOTYPES == 1,
  "phenotype_specific",
  "phenotype_shared"
)]

setorder(merged_regions, CHR, REGION_START, REGION_END)
merged_regions[, merged_region_id := sprintf(
  "conjFDR_region_%03d",
  seq_len(.N)
)]

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

merged_regions[, region_number := NULL]
setcolorder(
  merged_regions,
  c("merged_region_id", setdiff(names(merged_regions), "merged_region_id"))
)

fwrite(
  merged_regions,
  file.path(out_root, "conjfdr", "merged_regions.csv")
)

fwrite(
  original_loci_with_region,
  file.path(out_root, "conjfdr", "original_loci_with_region.csv")
)

cat("Original loci:", nrow(original_loci_with_region), "\n")
cat("Merged regions:", nrow(merged_regions), "\n")
cat("Trait-specific regions:", merged_regions[LOCUS_CLASS == "phenotype_specific", .N], "\n")
cat("Shared regions:", merged_regions[LOCUS_CLASS == "phenotype_shared", .N], "\n")
print(merged_regions[, .N, by = N_PHENOTYPES][order(N_PHENOTYPES)])
