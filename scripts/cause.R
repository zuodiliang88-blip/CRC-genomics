source("config/config.R")
plink_bin <- required_path("PLINK_BIN")
eur_bfile <- required_path("EUR_BFILE")
library(readr)
library(dplyr)
library(cause)
library(data.table)
library(ieugwasr)

exposure_files <- sort(list.files(factor_dir, full.names = TRUE))
outcome_file <- crc_file
out_root <- file.path(out_root, "CAUSE")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

for (exposure_file in exposure_files) {
  exposure_name <- sub("_noMHC.*$", "", basename(exposure_file))
  outcome_name <- tools::file_path_sans_ext(
    tools::file_path_sans_ext(basename(outcome_file))
  )
  analysis_name <- paste0(exposure_name, "_to_", outcome_name)
  analysis_dir <- file.path(out_root, analysis_name)
  dir.create(analysis_dir, showWarnings = FALSE, recursive = TRUE)

  tryCatch({
    # Step 1: Format Data for CAUSE ####
    X1 <- fread(exposure_file)
    X2 <- fread(outcome_file)
    X <- gwas_merge(
      X1, X2,
      snp_name_cols = c("SNP", "SNP"),
      beta_hat_cols = c("beta", "beta"),
      se_cols = c("se", "se"),
      A1_cols = c("A1", "A1"),
      A2_cols = c("A2", "A2")
    )

    # 添加 p 列；gwas_merge 函数不具有 pval 参数
    X1_snp_p <- X1[, .(SNP, pval)]
    X2_snp_p <- X2[, .(SNP, p)]
    colnames(X1_snp_p)[2] <- "p1"
    colnames(X2_snp_p)[2] <- "p2"
    colnames(X1_snp_p)[1] <- "snp"
    colnames(X2_snp_p)[1] <- "snp"
    # Newer CAUSE versions may add p1/p2 during gwas_merge. Remove those
    # before reproducing the senior script's joins from the source GWAS.
    X$p1 <- NULL
    X$p2 <- NULL
    a1 <- left_join(X, X1_snp_p, by = c("snp"))
    a2 <- left_join(a1, X2_snp_p, by = c("snp"))
    X <- a2

    # Step 2: Calculate nuisance parameters ####
    set.seed(100)
    varlist <- with(X, sample(snp, size = 1000000, replace = FALSE))
    params <- est_cause_params(X, varlist)

    # Step 3: LD Pruning ####
    dat <- X[, c("snp", "p1")]
    colnames(dat)[1:2] <- c("rsid", "pval")
    dat <- mutate(dat, pval = as.numeric(dat$pval))
    expo_clump <- ld_clump(
      dat = dat,
      clump_r2 = 0.01,
      clump_p = 1e-3,
      plink_bin = plink_bin,
      bfile = eur_bfile
    )
    X_clump <- X[X$snp %in% expo_clump$rsid, ]
    top_vars <- X_clump$snp

    # Step 4: Fit CAUSE ####
    res <- cause(X = X, variants = top_vars, param_ests = params, force = TRUE)

    # Step 5: Look at and save Results ####
    saveRDS(res, file = file.path(analysis_dir, paste0(analysis_name, "_cause_result.rds")))
    write.table(
      res$elpd,
      file = file.path(analysis_dir, paste0(analysis_name, "_elpd.txt")),
      sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
    )
    capture.output(
      summary(res, ci_size = 0.95),
      file = file.path(analysis_dir, paste0(analysis_name, "_summary.txt"))
    )
    capture.output(
      loo::pareto_k_table(res$loos[[2]]),
      file = file.path(analysis_dir, paste0(analysis_name, "_sharing_pareto_k.txt"))
    )
    capture.output(
      loo::pareto_k_table(res$loos[[3]]),
      file = file.path(analysis_dir, paste0(analysis_name, "_causal_pareto_k.txt"))
    )

    message("DONE: ", analysis_name)
  }, error = function(e) {
    message("ERROR in exposure file: ", exposure_file, " | ", e$message)
  })
}
