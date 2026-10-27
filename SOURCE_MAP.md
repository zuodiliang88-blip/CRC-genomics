# Original analysis sources

The following directories under the original study workspace were inspected to prepare this compact repository. They are provenance references, not required locations for a new run.

| Repository module | Original server source |
| --- | --- |
| LDSC | `LDSC/code/`, `LDSC/rg_code/`, `LDSC/h2_single_trait_20260816/code/` |
| GNOVA | `GNOVA/code/01_clean_gnova_inputs_38.R`, `GNOVA/code/*_CRC_gnova.sh` |
| SUPERGNOVA | `SUPERGNOVA/code/crc_13_traits/` |
| MiXeR | Server 2: `MiXeR/unimixer/code/*_unimixer.sh`, `MiXeR/bimixer/code/*_bimixer.sh` |
| pleioFDR and loci | `conjfdr2/*_crc/code/`, `conjfdr2/All_result/merge_conjFDR_loci_across_13_phenotypes.R` |
| FUSION TWAS | `TWAS/code/` |
| Colocalization | `coloc/FDR/code/` |
| Bidirectional MR | `MR/mr_*_modifiable_factors*.R`, `MR/filtering_audit/generate_MR_instrument_filtering_summary.R` |
| CAUSE | `cause/cause_13_modifiable_factors_to_crc.R` |
| MRlap | `MR/MRlap/code/` |

Only the main analysis flow is retained from the trait-specific MiXeR scripts; restart/queue and figure code were omitted. The original `MiXeR/01_make_noMHC_gwas_12_traits.R` is noted as upstream preparation but is not reproduced because it did not cover BMI and WHR. Single-cell and plotting scripts were excluded at the author's request.
