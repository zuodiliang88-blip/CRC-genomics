source("config/config.R")
plink_bin <- required_path("PLINK_BIN")
eur_bfile <- required_path("EUR_BFILE")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(TwoSampleMR)
  library(ieugwasr)
  library(plinkbinr)
  library(RadialMR)
})

out_root <- file.path(out_root, "MR")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
crc_name <- "colorectalcancer_pmid36539618"

short_trait <- function(x) {
  x <- sub("_noMHC.*$", "", x)
  x <- sub("_addchrBP_GRCh37$", "", x)
  x <- sub("_GRCh37$", "", x)
  x <- sub("_pmid[0-9]+$", "", x)
  x <- sub("_GCST[0-9]+$", "", x)
  x <- sub("_GSCAN_SmkInit_2022_EUR$", "", x)
  x <- sub("_GSCAN_SmkCes_2022_EUR$", "", x)
  x <- sub("_GSCAN_CigDay_2022_EUR$", "", x)
  x <- sub("_GSCAN_AgeSmk_2022_EUR$", "", x)
  x
}

safe_num <- function(x) {
  if (length(x) == 0) return(NA_real_)
  as.numeric(x[1])
}

ivw_summary <- function(dat) {
  if (is.null(dat) || nrow(dat) < 2) {
    return(data.table(ivw_b = NA_real_, ivw_se = NA_real_, ivw_p = NA_real_,
                      median_b = NA_real_, median_se = NA_real_, median_p = NA_real_))
  }
  res <- tryCatch(
    mr(dat, method_list = c("mr_ivw", "mr_weighted_median")),
    error = function(e) NULL
  )
  if (is.null(res) || nrow(res) == 0) {
    return(data.table(ivw_b = NA_real_, ivw_se = NA_real_, ivw_p = NA_real_,
                      median_b = NA_real_, median_se = NA_real_, median_p = NA_real_))
  }
  ivw <- res[res$method == "Inverse variance weighted", ]
  median <- res[res$method == "Weighted median", ]
  data.table(ivw_b = safe_num(ivw$b), ivw_se = safe_num(ivw$se), ivw_p = safe_num(ivw$pval),
             median_b = safe_num(median$b), median_se = safe_num(median$se),
             median_p = safe_num(median$pval))
}

het_summary <- function(dat) {
  if (is.null(dat) || nrow(dat) < 2) return(NA_real_)
  het <- tryCatch(mr_heterogeneity(dat), error = function(e) NULL)
  if (is.null(het) || nrow(het) == 0) return(NA_real_)
  safe_num(het$Q_pval[het$method == "Inverse variance weighted"])
}

egger_summary <- function(dat) {
  if (is.null(dat) || nrow(dat) < 3) return(NA_real_)
  ple <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)
  if (is.null(ple) || nrow(ple) == 0) return(NA_real_)
  safe_num(ple$pval)
}

remove_radial_outliers <- function(dat) {
  if (is.null(dat) || nrow(dat) < 3) {
    return(list(dat = dat, removed = 0L, status = "not_run_lt3_snps"))
  }
  radial <- tryCatch(
    ivw_radial(tsmr_to_rmr_format(dat), alpha = 0.05),
    error = function(e) e
  )
  if (inherits(radial, "error")) {
    return(list(dat = dat, removed = NA_integer_, status = paste0("radial_error: ", radial$message)))
  }
  outliers <- character(0)
  if (!is.null(radial$outliers) && nrow(radial$outliers) > 0 && "SNP" %in% names(radial$outliers)) {
    p_col <- intersect(c("p.value", "P.value", "pval", "P"), names(radial$outliers))[1]
    if (!is.na(p_col)) {
      outliers <- radial$outliers$SNP[radial$outliers[[p_col]] < 0.05]
    } else {
      outliers <- radial$outliers$SNP
    }
  }
  outliers <- unique(outliers[!is.na(outliers)])
  list(dat = dat[!dat$SNP %in% outliers, ], removed = length(outliers), status = "ok")
}

apply_steiger <- function(dat) {
  if (is.null(dat) || nrow(dat) == 0) {
    return(list(dat = dat, removed = 0L, status = "not_run_0_snps"))
  }
  st <- tryCatch(steiger_filtering(dat), error = function(e) e)
  if (inherits(st, "error")) {
    return(list(dat = dat, removed = NA_integer_, status = paste0("steiger_error: ", st$message)))
  }
  keep <- st$steiger_dir == TRUE | st$steiger_dir == "TRUE"
  keep[is.na(keep)] <- FALSE
  list(dat = st[keep, ], removed = sum(!keep), status = "ok")
}

summarise_pair <- function(dat, trait, direction, radial_status = NA_character_, steiger_status = NA_character_) {
  before <- dat$before
  after <- dat$after
  ivw_before <- ivw_summary(before)
  ivw_after <- ivw_summary(after)
  data.table(
    Trait = trait,
    Direction = direction,
    SNPs_before_filtering = ifelse(is.null(before), 0L, nrow(before)),
    Removed_by_radial = dat$removed_radial,
    Removed_by_Steiger = dat$removed_steiger,
    SNPs_after_filtering = ifelse(is.null(after), 0L, nrow(after)),
    IVW_beta_before_filtering = ivw_before$ivw_b,
    IVW_SE_before_filtering = ivw_before$ivw_se,
    IVW_P_before_filtering = ivw_before$ivw_p,
    IVW_beta_after_filtering = ivw_after$ivw_b,
    IVW_SE_after_filtering = ivw_after$ivw_se,
    IVW_P_after_filtering = ivw_after$ivw_p,
    Weighted_median_beta_after_filtering = ivw_after$median_b,
    Weighted_median_SE_after_filtering = ivw_after$median_se,
    Weighted_median_P_after_filtering = ivw_after$median_p,
    Cochran_Q_P_before_filtering = het_summary(before),
    Cochran_Q_P_after_filtering = het_summary(after),
    Egger_intercept_P_before_filtering = egger_summary(before),
    Egger_intercept_P_after_filtering = egger_summary(after),
    radial_status = radial_status,
    steiger_status = steiger_status
  )
}

prepare_exposure <- function(file, p_col, name) {
  exposure <- read_exposure_data(
    file,
    sep = "\t",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "freq",
    pval_col = p_col,
    samplesize_col = "n"
  )
  exposure <- exposure[exposure$pval.exposure < 5e-8, ]
  clump_input <- exposure[, c("SNP", "pval.exposure")]
  colnames(clump_input) <- c("rsid", "pval")
  clump_input <- mutate(clump_input, pval = as.numeric(pval))
  clumped <- ld_clump(
    dat = clump_input,
    clump_kb = 10000,
    clump_r2 = 0.001,
    plink_bin = plink_bin,
    bfile = eur_bfile
  )
  exposure <- exposure[exposure$SNP %in% clumped$rsid, ]
  exposure$r2 <- get_r_from_pn(exposure$pval.exposure, exposure$samplesize.exposure)
  exposure <- mutate(exposure, F = beta.exposure^2 / se.exposure^2)
  exposure <- subset(exposure, F >= 10)
  exposure$exposure <- name
  exposure
}

harmonised_sets <- function(exposure, outcome) {
  dat <- harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
  if ("mr_keep" %in% names(dat)) dat <- dat[dat$mr_keep, ]
  dat
}

filter_sets <- function(dat) {
  before <- dat
  radial <- remove_radial_outliers(before)
  steiger <- apply_steiger(radial$dat)
  list(
    before = before,
    after = steiger$dat,
    removed_radial = radial$removed,
    removed_steiger = steiger$removed,
    radial_status = radial$status,
    steiger_status = steiger$status
  )
}

factor_files <- sort(list.files(factor_dir, full.names = TRUE))

crc_outcome <- read_outcome_data(
  filename = crc_file,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "freq",
  pval_col = "p",
  samplesize_col = "n"
)
crc_outcome$outcome <- crc_name

message("Preparing CRC instruments for reverse-direction MR...")
crc_exposure <- prepare_exposure(crc_file, p_col = "p", name = crc_name)

results <- list()

for (file in factor_files) {
  factor_raw <- sub("_noMHC.*$", "", basename(file))
  factor_short <- short_trait(factor_raw)
  message("Processing factor -> CRC: ", factor_raw)

  tryCatch({
    factor_exposure <- prepare_exposure(file, p_col = "pval", name = factor_raw)
    crc_subset <- crc_outcome[crc_outcome$SNP %in% factor_exposure$SNP, ]
    dat <- harmonised_sets(factor_exposure, crc_subset)
    filtered <- filter_sets(dat)
    results[[paste0(factor_raw, "__factor_to_crc")]] <- summarise_pair(
      filtered,
      trait = factor_short,
      direction = "factor_to_CRC",
      radial_status = filtered$radial_status,
      steiger_status = filtered$steiger_status
    )
  }, error = function(e) {
    results[[paste0(factor_raw, "__factor_to_crc_error")]] <<- data.table(
      Trait = factor_short,
      Direction = "factor_to_CRC",
      error = e$message
    )
  })

  message("Processing CRC -> factor: ", factor_raw)
  tryCatch({
    factor_outcome <- read_outcome_data(
      filename = file,
      sep = "\t",
      snp_col = "SNP",
      beta_col = "beta",
      se_col = "se",
      effect_allele_col = "A1",
      other_allele_col = "A2",
      eaf_col = "freq",
      pval_col = "pval",
      samplesize_col = "n"
    )
    factor_outcome$outcome <- factor_raw
    factor_subset <- factor_outcome[factor_outcome$SNP %in% crc_exposure$SNP, ]
    dat <- harmonised_sets(crc_exposure, factor_subset)
    filtered <- filter_sets(dat)
    results[[paste0(factor_raw, "__crc_to_factor")]] <- summarise_pair(
      filtered,
      trait = factor_short,
      direction = "CRC_to_factor",
      radial_status = filtered$radial_status,
      steiger_status = filtered$steiger_status
    )
  }, error = function(e) {
    results[[paste0(factor_raw, "__crc_to_factor_error")]] <<- data.table(
      Trait = factor_short,
      Direction = "CRC_to_factor",
      error = e$message
    )
  })
}

summary_table <- rbindlist(results, fill = TRUE)
setorder(summary_table, Trait, Direction)

tsv_file <- file.path(out_root, "MR_instrument_filtering_before_after_summary.tsv")
fwrite(summary_table, tsv_file, sep = "\t")
message("Wrote: ", tsv_file)
message("Rows: ", nrow(summary_table))
