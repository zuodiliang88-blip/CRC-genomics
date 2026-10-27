library(data.table)

source("config/config.R")
in_dir <- file.path(data_dir, "factors_GRCh37")
out_dir <- file.path(data_dir, "twas_sumstats")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
files <- list.files(in_dir, full.names = TRUE)
files <- c(files, crc_file)

for (f in files) {
  d <- fread(f)
  names(d) <- toupper(names(d))
  d <- d[is.finite(as.numeric(BETA)) & is.finite(as.numeric(SE)) & as.numeric(SE) > 0]
  d[, Z := as.numeric(BETA) / as.numeric(SE)]
  d <- d[, .(SNP, A1, A2, Z, N)]
  d <- na.omit(d)
  d <- d[!duplicated(SNP)]
  trait <- if (normalizePath(f) == normalizePath(crc_file)) "CRC" else sub("\\.(txt|tsv)\\.gz$", "", basename(f))
  out <- file.path(out_dir, paste0(trait, ".sumstats.gz"))
  fwrite(d, out, sep = "\t")
}
