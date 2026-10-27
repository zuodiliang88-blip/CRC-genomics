library(data.table)
source("config/config.R")
direction <- fread("config/rg_direction.tsv")
twas_dir <- file.path(out_root, "twas")
read_trait <- function(trait) {
  files <- list.files(file.path(twas_dir, trait), pattern = "\\.dat$", full.names = TRUE)
  if (!length(files)) stop("Missing FUSION output for ", trait)
  x <- rbindlist(lapply(files, fread), fill = TRUE)
  stopifnot(all(c("ID", "TWAS.Z", "TWAS.P") %in% names(x)))
  x <- x[!is.na(ID) & is.finite(`TWAS.P`) & `TWAS.P` >= 0 & `TWAS.P` <= 1]
  x[, FDR := p.adjust(`TWAS.P`, "BH")]
  x[, .(ID, TWAS_Z = `TWAS.Z`, TWAS_P = `TWAS.P`, FDR)]
}
crc <- read_trait("CRC")
setnames(crc, c("TWAS_Z", "TWAS_P", "FDR"), c("CRC_Z", "CRC_P", "CRC_FDR"))
shared <- rbindlist(lapply(seq_len(nrow(direction)), function(i) {
  trait <- direction$Trait[i]
  x <- merge(read_trait(trait), crc, by = "ID")
  x[, `:=`(Trait = trait, rg_sign = direction$rg_sign[i])]
  x[, direction_concordant := sign(TWAS_Z * CRC_Z) == rg_sign]
  x
}), fill = TRUE)
fwrite(shared, file.path(out_root, "all_shared_twas.tsv"), sep = "\t")
eligible <- shared[FDR < 0.05 & CRC_FDR < 0.05 & direction_concordant == TRUE]
fwrite(eligible[, .(Trait, Gene_ID = ID)], file.path(out_root, "coloc_candidates.tsv"), sep = "\t")
