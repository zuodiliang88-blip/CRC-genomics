library(data.table)
library(coloc)
source("config/config.R")

trait <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(trait) || !nzchar(trait)) stop("Usage: Rscript scripts/coloc.R TRAIT_ID")
types <- fread("config/coloc_trait_types.tsv", na.strings = c("", "NA"))
row <- types[Trait == trait]
if (nrow(row) != 1L) stop("Trait missing or duplicated in coloc_trait_types.tsv")
if (row$type == "cc" && (!is.finite(row$case_fraction) || row$case_fraction <= 0 || row$case_fraction >= 1)) {
  stop("Case-control traits require a valid case_fraction in coloc_trait_types.tsv")
}

read_gwas <- function(file) {
  x <- fread(file)
  setnames(x, toupper(names(x)))
  pcol <- intersect(c("PVAL", "P"), names(x))[1]
  fcol <- intersect(c("FREQ", "FRQ"), names(x))[1]
  if (is.na(pcol) || is.na(fcol)) stop("Missing P or allele frequency in ", file)
  x <- x[, .(SNP = as.character(SNP), A1 = toupper(as.character(A1)),
             A2 = toupper(as.character(A2)), p = as.numeric(get(pcol)),
             n = as.numeric(N), freq = as.numeric(get(fcol)))]
  x <- x[!duplicated(SNP) & nzchar(SNP) & is.finite(p) & p > 0 & p <= 1 &
           is.finite(n) & n > 0 & is.finite(freq) & freq > 0 & freq < 1]
  x
}

trait_gwas <- read_gwas(file.path(data_dir, "factors_GRCh37", paste0(trait, ".txt.gz")))
crc_gwas <- read_gwas(crc_file)
maf <- fread(required_path("EQTLGEN_MAF"))[, .(SNP = as.character(SNP), MAF = as.numeric(MAF))]
maf <- maf[!duplicated(SNP) & is.finite(MAF) & MAF > 0 & MAF < 0.5]
setkey(maf, SNP)
genes <- unique(fread(file.path(out_root, "coloc_candidates.tsv"))[Trait == trait, Gene_ID])
if (!length(genes)) stop("No directionally concordant shared TWAS genes for ", trait)

make_pair <- function(eqtl, gwas) {
  x <- merge(eqtl, gwas, by = "SNP", suffixes = c("_eqtl", "_gwas"))
  x <- x[(A1_eqtl == A1_gwas & A2_eqtl == A2_gwas) |
           (A1_eqtl == A2_gwas & A2_eqtl == A1_gwas)]
  x <- x[!paste0(A1_eqtl, A2_eqtl) %chin% c("AT", "TA", "CG", "GC")]
  x[, MAF := maf[J(SNP), MAF]]
  x[is.finite(MAF) & MAF > 0 & MAF < 0.5]
}

run_coloc <- function(x, type, case_fraction = NA_real_) {
  if (nrow(x) < 2L) return(NULL)
  dataset2 <- list(snp = x$SNP, pvalues = x$p, N = x$n, type = type)
  if (type == "cc") dataset2$s <- case_fraction
  coloc.abf(
    dataset1 = list(snp = x$SNP, pvalues = x$Pvalue, N = x$NrSamples, type = "quant"),
    dataset2 = dataset2, MAF = x$MAF
  )$summary
}

results <- list()
for (gene_id in genes) {
  eqtl <- fread(file.path(required_path("EQTLGEN_DIR"), paste0(gene_id, ".txt")))
  eqtl <- eqtl[, .(SNP = as.character(SNP), CHR = as.integer(SNPChr),
                   BP = as.numeric(SNPPos), A1 = toupper(as.character(AssessedAllele)),
                   A2 = toupper(as.character(OtherAllele)), Pvalue = as.numeric(Pvalue),
                   NrSamples = as.numeric(NrSamples), GeneSymbol = as.character(GeneSymbol))]
  eqtl <- eqtl[!duplicated(SNP) & nzchar(SNP) & is.finite(Pvalue) &
                 Pvalue > 0 & Pvalue <= 1 & is.finite(NrSamples) & NrSamples > 0]
  if (!nrow(eqtl)) next
  lead <- eqtl[which.min(Pvalue)]
  region <- eqtl[CHR == lead$CHR & BP > lead$BP - 250000 & BP < lead$BP + 250000]
  for (comparison in c("Trait", "CRC")) {
    is_crc <- comparison == "CRC"
    x <- make_pair(region, if (is_crc) crc_gwas else trait_gwas)
    summary <- run_coloc(x, if (is_crc) "cc" else row$type,
                         if (is_crc) 78473 / 185616 else row$case_fraction)
    if (is.null(summary)) next
    results[[length(results) + 1L]] <- data.table(
      Trait = trait, Gene_ID = gene_id, GeneSymbol = lead$GeneSymbol,
      comparison = comparison, nsnps = as.integer(summary["nsnps"]),
      PP.H0 = as.numeric(summary["PP.H0.abf"]),
      PP.H1 = as.numeric(summary["PP.H1.abf"]),
      PP.H2 = as.numeric(summary["PP.H2.abf"]),
      PP.H3 = as.numeric(summary["PP.H3.abf"]),
      PP.H4 = as.numeric(summary["PP.H4.abf"])
    )
  }
}
dir.create(file.path(out_root, "coloc"), recursive = TRUE, showWarnings = FALSE)
fwrite(rbindlist(results, fill = TRUE), file.path(out_root, "coloc", paste0(trait, ".tsv")), sep = "\t")
