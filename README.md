# ICBO

Code for **Efficient Optimization of Genotype Pairs for Intercropping using Genomic Prediction and Bayesian Optimization**, by Sei Kinoshita and Hiroyoshi Iwata.

## Workflow

Run the scripts from the repository root in the following order. Inputs are read from `data/`; analysis outputs are saved under `midstream/`.

| Step | Scripts in `scripts/` |
|---|---|
| Genotype preparation and phenotype simulation | `0.1.modifyGeno.r` → `0.2.sampleTrueMrkEff.r` → `0.3.generatePheno.r` |
| Prediction model evaluation | `0.4.selTrainDataGP.r` → `1.0.GPFitting.r` → `1.1.GPPred.r` → `1.2.GPAccuracy.r` |
| Genotype-pair search | `2.0.BO.R` → `2.1.BOEval.R` |
| Tables and figures | `3.0.makeTables.R` (Tables S1–S3) → `3.2.makeFigures.R` (Figures 2–5, S1–S3, S5–S10) |

Figure sections are labelled in `3.2.makeFigures.R`. Figures 1 and S4 are manually drawn schematics. Functions shared by the scripts are in `scripts/function_scripts/`.

## Genotype data

Sorghum and soybean genotype data were obtained from [Takanashi et al. (2022)](https://doi.org/10.1093/pcp/pcac057) and [Kajiya-Kanegae et al. (2021)](https://doi.org/10.1093/dnares/dsaa032), respectively. Data and results are not included in this repository.

Prepare the MAF/LD-filtered genotype matrices as described in the manuscript and place them in `data/genotype/`:

- `0.1.Sorghum237_geno_Maf_0.1_LD_0.6.rds`
- `0.1.Soybean198_geno_Maf_0.1_LD_0.6.rds`

Matrices use 0/1/2 coding, with accession IDs as row names and SNP IDs as column names. The selected 79/66 accession IDs are read from `data/SelectedIndList_kMedoids_nIndSb_79_nIndGm_66.rds` (a list with `Sb` and `Gm` vectors).

## Environment

R 4.5.1 and BGLR 1.1.4; package versions are listed in `package_versions.tsv`. Launch parallel analyses with `Rscript` from a terminal; worker counts are set in each script.
