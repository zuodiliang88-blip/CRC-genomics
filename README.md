# CRC genomics analysis code

Analysis scripts for the study of shared genetic architecture between colorectal cancer (CRC) and 37 potentially modifiable or socially patterned factors. This repository is a compact, path-configurable presentation of the analysis code used on the study server. It does not contain GWAS summary statistics, LD reference panels, eQTL files, single-cell data, generated results, or figure scripts.

## Scope

| Analysis | Script | Original analysis |
| --- | --- | --- |
| GWAS preparation and global correlation | `01_global_genetics.sh`, `00_clean_gnova.R` | LDSC h2/rg and GNOVA, 37 CRC-factor pairs |
| Local genetic covariance | `02_supergnova.sh`, `04_supergnova_fdr.R` | SUPERGNOVA, 2,353-region EUR partition, BH within trait pair |
| Polygenic overlap | `MiXeR/README.md` | Placeholder; code to be supplied later |
| Cross-trait enrichment | `06_pleiofdr.sh`, `config/pleiofdr.example.txt` | pleioFDR condFDR/conjFDR; external tool required |
| Locus definition | `07_locus_clump.sh`, `merge_conjfdr_loci.R` | FUMA-like clumping and 250-kb cross-trait merging |
| Gene-level analysis | `twas_prepare.R`, `03_twas.sh`, `05_shared_twas.R`, `coloc.R` | FUSION TWAS and eQTLGen colocalization |
| Directionality analyses | `mr_bidirectional.R`, `cause.R`, `mrlap.R` | Bidirectional MR, CAUSE and MRlap |

Single-cell analysis and all plotting code are intentionally excluded.

## Setup

1. Run from the repository root. Copy `config/paths.example.sh` to a private file, set absolute paths, and source it. Do not commit the private file.
2. Provide harmonized GRCh37 GWAS files and the external tools/reference panels listed below. The table of GWAS sources and accessions is Supplementary Table S1 of the manuscript.
3. Install the R packages used in each script (`data.table`, `dplyr`, `TwoSampleMR`, `ieugwasr`, `RadialMR`, `cause`, `MRlap`, and `coloc`) and their system dependencies. FUSION, LDSC, GNOVA, SUPERGNOVA, pleioFDR, MATLAB, PLINK and 1000 Genomes EUR reference data are external dependencies, not vendored here.

`CRC_DATA_DIR` should contain `gwas/` (LDSC input), `factors_GRCh37/` (TWAS/coloc input), and `factors_noMHC/` (MR/CAUSE/MRlap input). GWAS tables should provide SNP, alleles, beta, standard error, P value, allele frequency and sample size; see the input-column handling in each script. Keep results under `CRC_DATA_DIR/results/`.

Examples:

```bash
source config/paths.local.sh
bash scripts/01_global_genetics.sh munge BMI_GCST009004
bash scripts/01_global_genetics.sh h2 BMI_GCST009004
bash scripts/01_global_genetics.sh rg BMI_GCST009004
Rscript scripts/00_clean_gnova.R
bash scripts/01_global_genetics.sh gnova BMI_GCST009004
bash scripts/02_supergnova.sh BMI_GCST009004 697734
Rscript scripts/04_supergnova_fdr.R
Rscript scripts/twas_prepare.R
bash scripts/03_twas.sh BMI_GCST009004_addchrBP_GRCh37
Rscript scripts/05_shared_twas.R
Rscript scripts/mr_bidirectional.R
Rscript scripts/cause.R
Rscript scripts/mrlap.R BMI_GCST009004
```

For the CRC LDSC input, run the `munge` mode with the CRC file and its P-value column (`p`) after placing it at `gwas/CRC.txt.gz`. Sample sizes for SUPERGNOVA should be taken from the source GWAS rather than inferred from the filename.

## Method notes

- LDSC uses HapMap3 SNPs and EUR LD scores with an unconstrained bivariate intercept. GNOVA preprocessing follows the original exclusion of nonfinite values and `|Z| >= 80`.
- SUPERGNOVA is shown for the 2,353-region partition used for inference. FDR adjustment is performed separately for each CRC-factor pair on the covariance P values.
- pleioFDR uses 500 random LD-pruning iterations, MAF >= 0.005, and excludes the MHC and 8p23.1 inversion regions from model fitting. The study-server runner additionally set `options.correct_for_sample_overlap = true` in the pleioFDR MATLAB options. The external pleioFDR code and its locally modified runner are not redistributed here; the wrapper requires that same setting. The separate Python conversion and clumping utilities distributed with pleioFDR are also required.
- FUSION used eQTLGen blood prediction models. Shared TWAS genes were selected at BH FDR < 0.05 in each phenotype, then classified by the sign of the CRC and factor TWAS Z scores relative to genome-wide rg. Only directionally concordant genes were sent to coloc.
- `mr_bidirectional.R` follows the instrument-filtering audit requested for the final analysis: genome-wide significant instruments, EUR LD clumping at r2 < 0.001 within 10 Mb, F >= 10, harmonization, radial outlier removal at P < 0.05, then Steiger filtering. It reports pre/post IVW estimates and post-filter weighted median, heterogeneity and Egger tests. This differs from an older forward-MR server script where radial removal was commented out.
- CAUSE uses 1,000,000 sampled variants for nuisance-parameter estimation and LD pruning at r2 = 0.01, P = 1e-3. MRlap uses P < 5e-8 and its 10-Mb/r2 = 0.001 pruning settings. The manuscript reports the seven significant forward associations; the scripts can process all available factors.

## Coloc correction pending rerun

The original study-server coloc scripts used `type = "quant"` for every factor GWAS, including binary outcomes. The unified `coloc.R` here assigns `type = "cc"` to smoking initiation, smoking cessation and moderate-to-vigorous physical activity, and requires the corresponding case fractions (`s`) in `config/coloc_trait_types.tsv`. Those fractions are **not** provided in Supplementary Table S1 and have deliberately been left as `NA` rather than guessed. The corrected script cannot run for these three traits until the original case/control counts are supplied. Its PP.H4 estimates must be rerun and reconciled with the manuscript before claiming exact numerical reproduction. CRC uses `type = "cc"` with 78,473 cases among 185,616 individuals. The other factors remain `quant`, matching their study measurement scales.

## Data and provenance

GWAS sources and accessions are in Supplementary Table S1. [eQTLGen TWAS models](https://zenodo.org/records/7068381) and [eQTLGen cis-eQTL data](https://eqtlgen.org/cis-eqtls.html) are external inputs. No primary or participant-level data are redistributed. `SOURCE_MAP.md` identifies the original server scripts used to assemble each module. The scripts have been syntax-checked, but full reruns require the original GWAS and reference data and have not been performed as part of repository preparation.
