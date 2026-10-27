library(data.table)
source("config/config.R")
files <- list.files(file.path(out_root, "supergnova"), pattern = "\\.2353\\.txt$", full.names = TRUE)
stopifnot(length(files) > 0)
results <- lapply(files, function(file) {
  x <- fread(file)
  stopifnot(all(c("chr", "start", "end", "rho", "p") %in% names(x)))
  x[, trait := sub("\\.CRC\\.2353\\.txt$", "", basename(file))]
  x[, FDR := p.adjust(p, method = "BH")]
  x[, direction := fifelse(rho > 0, "positive", "negative")]
  x
})
fwrite(rbindlist(results, fill = TRUE), file.path(out_root, "supergnova_fdr.tsv"), sep = "\t")
