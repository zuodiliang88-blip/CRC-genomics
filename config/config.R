required_path <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) stop("Set ", name, " before running this script")
  value
}

data_dir <- required_path("CRC_DATA_DIR")
crc_file <- required_path("CRC_GWAS")
factor_dir <- file.path(data_dir, "factors_noMHC")
plink_bin <- Sys.getenv("PLINK_BIN", unset = "")
eur_bfile <- Sys.getenv("EUR_BFILE", unset = "")
out_root <- file.path(data_dir, "results")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
