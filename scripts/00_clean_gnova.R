library(data.table)
source("config/config.R")
input_dir <- file.path(data_dir, "ldsc")
output_dir <- file.path(data_dir, "gnova")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
files <- list.files(input_dir, pattern = "\\.sumstats\\.gz$", full.names = TRUE)
if (!length(files)) stop("No LDSC sumstats in ", input_dir)
qc <- lapply(files, function(file) {
  x <- fread(cmd = paste("gzip -dc", shQuote(file)), na.strings = c("NA", "NaN", "Inf", "-Inf", ""))
  before <- nrow(x)
  x <- x[, .(SNP, A1, A2, Z = as.numeric(Z), N = as.numeric(N))]
  x <- x[!is.na(SNP) & !is.na(A1) & !is.na(A2) & is.finite(Z) &
           abs(Z) < 80 & is.finite(N) & N > 0]
  trait <- sub("\\.sumstats\\.gz$", "", basename(file))
  fwrite(x, file.path(output_dir, paste0(trait, ".clean.sumstats.gz")), sep = "\t")
  data.table(trait = trait, before = before, after = nrow(x))
})
fwrite(rbindlist(qc), file.path(out_root, "gnova_qc.tsv"), sep = "\t")
