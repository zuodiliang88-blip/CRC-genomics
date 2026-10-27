source("config/config.R")
plink_bin <- required_path("PLINK_BIN")
eur_bfile <- required_path("EUR_BFILE")
library(MRlap)
library(data.table)

#========================
# 路径和固定参数
#========================
exp_name <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(exp_name) || !nzchar(exp_name)) stop("Usage: Rscript scripts/mrlap.R TRAIT_ID")
result_dir <- file.path(out_root, "MRlap", exp_name)
tmp_dir <- file.path(out_root, "MRlap_tmp", exp_name)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
setwd(tmp_dir)

exposure_file <- file.path(factor_dir, paste0(exp_name, "_noMHC.txt.gz"))
outcome_file <- crc_file

ld_path <- required_path("LDSC_LD")
hm3_path <- required_path("HM3_SNPLIST")

bfile_ref <- eur_bfile

#========================
# 读取 outcome
#========================
outcome <- fread(outcome_file, data.table = FALSE)
outcome_name <- "colorectalcancer_pmid36539618"

setnames(
  outcome,
  old = c("SNP", "A1", "A2", "beta", "se", "n", "CHR", "BP"),
  new = c("snp", "a1", "a2", "beta", "se", "N", "chr", "pos")
)

outcome$snp <- as.character(outcome$snp)
outcome$a1 <- toupper(as.character(outcome$a1))
outcome$a2 <- toupper(as.character(outcome$a2))
outcome$beta <- as.numeric(outcome$beta)
outcome$se <- as.numeric(outcome$se)
outcome$N <- as.numeric(outcome$N)
outcome$chr <- as.integer(outcome$chr)
outcome$pos <- as.integer(outcome$pos)

outcome <- outcome[
  grepl("^rs[0-9]+$", outcome$snp) &
    outcome$a1 %in% c("A", "C", "G", "T") &
    outcome$a2 %in% c("A", "C", "G", "T") &
    outcome$a1 != outcome$a2 &
    is.finite(outcome$beta) & is.finite(outcome$se) & outcome$se > 0 &
    is.finite(outcome$N) & outcome$N > 0 &
    outcome$chr %in% 1:22 & is.finite(outcome$pos) & outcome$pos > 0,
  c("snp", "a1", "a2", "beta", "se", "N", "chr", "pos")
]
outcome <- outcome[!duplicated(outcome$snp), ]

#========================
# 读取 exposure
#========================
exposure <- fread(exposure_file, data.table = FALSE)

setnames(
  exposure,
  old = c("SNP", "A1", "A2", "beta", "se", "n", "chr", "bp"),
  new = c("snp", "a1", "a2", "beta", "se", "N", "chr", "pos")
)

exposure$snp <- as.character(exposure$snp)
exposure$a1 <- toupper(as.character(exposure$a1))
exposure$a2 <- toupper(as.character(exposure$a2))
exposure$beta <- as.numeric(exposure$beta)
exposure$se <- as.numeric(exposure$se)
exposure$N <- as.numeric(exposure$N)
exposure$chr <- as.integer(exposure$chr)
exposure$pos <- as.integer(exposure$pos)

exposure <- exposure[
  grepl("^rs[0-9]+$", exposure$snp) &
    exposure$a1 %in% c("A", "C", "G", "T") &
    exposure$a2 %in% c("A", "C", "G", "T") &
    exposure$a1 != exposure$a2 &
    is.finite(exposure$beta) & is.finite(exposure$se) & exposure$se > 0 &
    is.finite(exposure$N) & exposure$N > 0 &
    exposure$chr %in% 1:22 & is.finite(exposure$pos) & exposure$pos > 0,
  c("snp", "a1", "a2", "beta", "se", "N", "chr", "pos")
]
exposure <- exposure[!duplicated(exposure$snp), ]

#========================
# MRlap
#========================
res <- MRlap(
  exposure = exposure,
  exposure_name = exp_name,
  outcome = outcome,
  outcome_name = outcome_name,

  ld = ld_path,
  hm3 = hm3_path,

  do_pruning = TRUE,
  MR_threshold = 5e-8,
  MR_pruning_dist = 10000,
  MR_pruning_LD = 0.001,

  MR_plink = plink_bin,
  MR_bfile = bfile_ref,

  verbose = FALSE
)

#========================
# 整理并保存结果
#========================
MRcorr <- res$MRcorrection
MRcorr <- MRcorr[!names(MRcorr) %in% "IVs"]

LDSC <- res$LDSC
GA <- res$GeneticArchitecture

res_df <- data.frame(
  exposure = exp_name,
  outcome = outcome_name,
  t(unlist(MRcorr)),
  t(unlist(LDSC)),
  t(unlist(GA)),
  row.names = NULL,
  check.names = FALSE
)

write.table(
  res_df,
  file = file.path(result_dir, paste0(exp_name, "_MRlap_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

res_df
