# ICBO: genomic prediction and Bayesian optimization for intercropping

Analysis code corresponding to the submitted revised manuscript, its Supplementary File, and the illustrative MCMC checks in the reviewer response. The release retains Models 1–5 for genomic prediction and Model 3 with Mu, PI, EI and UCB for optimization. Earlier exploratory figures, full-panel fitting diagnostics, the unused one-shot comparison, and Monte Carlo random-search runs are excluded from the public workflow.

This is a **code-only update**. No additional RDS inputs or results, manuscripts, reviewer correspondence, figures or MCMC samples are distributed. The three genotype/marker RDS files already public in the earlier repository remain unchanged at the repository root. A fresh clone alone cannot reproduce the submitted numerical results: the selected-accession IDs and historical result archives listed below are not distributed.

## Redraw figures and tables from local results

Run commands from the repository root in separate R processes. The default evaluation scripts require locally supplied historical archives (`data/publicationResults.rds`, and additionally `data/2.1.publicationSelection.rds` for BO); **they do not fit GP models or run BO**. These archives are not included in this repository. If they are absent, the scripts stop with an explanation. After a new full analysis, use `ICBO.rebuildSummaries=TRUE` as shown below instead.

```sh
Rscript scripts/1.2.GPAccuracy.r
Rscript scripts/2.1.BOEval.R
Rscript scripts/3.0.makeTables.R
```

For Figures S1–S2, first supply the selected-accession ID file described below, then prepare the genotype panels and evaluate the subset:

```sh
Rscript scripts/0.1.modifyGeno.r
Rscript scripts/0.1.evaluateSubset.r
```

Install the dependencies if needed:

```r
install.packages(c("BGLR", "gaston", "rrBLUP", "pbmcapply", "ggplot2", "scales",
                   "dplyr", "tidyr", "multcompView", "coda", "gtable"))
```

The tested environment is R 4.5.1 on Ubuntu 22.04.2, including BGLR 1.1.4. See [package_versions.tsv](package_versions.tsv) for direct/transitive dependency versions and [sessionInfo.txt](sessionInfo.txt) for the **2026-11-06 validation session**, not a retrospectively reconstructed original worker session. Installing today's package releases need not reproduce this environment. PDF generation requires an R build with Cairo support.

## Script-to-figure correspondence

Figure numbers refer to the submitted revision. Each PDF contains only the figure indicated in its filename. Multipage PDFs separate panels for subsequent manuscript layout; they do not reproduce the Word document's page layout or annotations added during artwork assembly.

| Manuscript item | Script | Output under `midstream/` | Panel/page interpretation |
|---|---|---|---|
| Figure 1 | No analysis script | Manual schematic in the manuscript | Conceptual BO workflow; the implementation is `2.0.BO.R` |
| Figure 2 | `1.2.GPAccuracy.r` | `1.2.GP_Accuracy/1.2.Figure2_PredictiveAbility.pdf` | Pages 1–2: h² = 0.3, 0.7; pooled sorghum and soybean |
| Figure 3 | `2.1.BOEval.R` | `2.1.BO_Evaluation/2.1.Figure3_CycleAdvantage.pdf` | Pages 1–2: h² = 0.3, 0.7; nine scenarios each |
| Figure 4 | `2.1.BOEval.R` | `2.1.BO_Evaluation/2.1.Figure4_CommonGPSelection.pdf` | Pages 1–2: ρ = −0.7, +0.7; h² = 0.7, γ = 0.2 |
| Figure 5 | `2.1.BOEval.R` | `2.1.BO_Evaluation/2.1.Figure5_DiscoveryRanks_SelectionDifference.pdf` | Page 1: panel a, discovery ranks; page 2: panel b, selected pairs differing from Mu |
| Figure S1 | `0.1.evaluateSubset.r` | `0.1.FigureS1_PCA.pdf` | Full-panel GRM PCA; sorghum and soybean, selected accessions highlighted |
| Figure S2 | `0.1.evaluateSubset.r` | `0.1.FigureS2_SubsetComparison.pdf` | MAF spectra and off-diagonal GRM distributions; 3 × 2 layout |
| Figure S3 | `1.2.GPAccuracy.r` | `1.2.GP_Accuracy/1.2.FigureS3_RealizedCorrelation.pdf` | Pages 1–2: sorghum, soybean; target versus realized correlation |
| Figure S4 | No analysis script | Manual schematic in the supplement | Design-D/Design-B illustration; actual training sets are generated/loaded by `0.4.selTrainDataGP.r` |
| Figure S5 | `1.2.GPAccuracy.r` | `1.2.GP_Accuracy/1.2.FigureS5_CorrelationRecovery.pdf` | Design-D; pages 1–2: sorghum, soybean |
| Figure S6 | `1.2.GPAccuracy.r` | `1.2.GP_Accuracy/1.2.FigureS6_CorrelationRecovery.pdf` | Design-B; pages 1–2: sorghum, soybean |
| Figure S7 | `1.2.GPAccuracy.r` | `1.2.GP_Accuracy/1.2.FigureS7_PredictiveAbility.pdf` | Sorghum; pages 1–2: h² = 0.3, 0.7 |
| Figure S8 | `1.2.GPAccuracy.r` | `1.2.GP_Accuracy/1.2.FigureS8_PredictiveAbility.pdf` | Soybean; pages 1–2: h² = 0.3, 0.7 |
| Figure S9 | `2.1.BOEval.R` | `2.1.BO_Evaluation/2.1.FigureS9_Discovery_CycleSavings.pdf` | Page 1: discovery curves and analytical Random; page 2: cycles saved relative to Mu |
| Figure S10 | `2.1.BOEval.R` | `2.1.BO_Evaluation/2.1.FigureS10_OptimumOutsideMeanTop56.pdf` | Fraction of runs selecting an optimum outside the posterior-mean top 56 |
| Table S1 | `3.0.makeTables.R` | `TableS1_PredictiveAbility.csv` | Mean and SE for each scenario, design and model |
| Table S2 | `3.0.makeTables.R` | `TableS2_BOCycles.csv` | Mean and SE of discovery cycles for each scenario and method |
| Table S3 | `3.0.makeTables.R` | `TableS3_GPTime.csv`, `TableS3_BOTime.csv` | Seconds per unique GP fit/prediction and per optimization run |
| Response Figures 1–2, optional | `3.1.checkMCMC.R` | `ResponseFigure1_GP_RunningMeans.pdf`, `ResponseFigure2_BO_RunningMeans.pdf` | Requires saved fits/traces from a full rerun; no new fits are performed by this check |

Figure S1 is reconstructed from the documented full-panel GRM eigendecomposition: the original plotting script was not available. PC signs and coordinate scaling are conventional; this is not a claim of pixel-identical reproduction of the original artwork. Figures 1 and S4 are explicitly manual illustrations, not missing numerical analyses.

## Data requirements and availability

[INPUT_SHA256.tsv](INPUT_SHA256.tsv) records checksums for the three **previously public** files below. They remain at their original root paths; script 0.1 also accepts local copies under `data/genotype/`. Genotypes use 0/1/2 coding and these filtered panels contain no heterozygous marker calls. `Sb` denotes sorghum and `Gm` soybean.

| Existing public file | Contents |
|---|---|
| `0.1.Sorghum237_geno_Maf_0.1_LD_0.6.rds` | Filtered full sorghum panel |
| `0.1.Soybean198_geno_Maf_0.1_LD_0.6.rds` | Filtered full soybean panel |
| `0.1.Sorghum237_markerInfo_Maf_0.1_LD_0.6.rds` | Sorghum marker positions; also regenerated by script 0.1 |

The following seven local files are **not distributed**. Their paths document the inputs accepted by the scripts; no accession IDs, seeds or numerical archives are embedded in the source code as substitutes.

| Local file | Role and behavior if absent |
|---|---|
| `data/SelectedIndList_kMedoids_nIndSb_79_nIndGm_66.rds` | Required for script 0.1: named list with character vectors `Sb` (79 unique IDs) and `Gm` (66 unique IDs), matching full-panel matrix row names. The exact study IDs are needed to reproduce its subset. |
| `data/0.2.genotypeEffects_used.rds` | Optional historical genotype-level DGE/IGE products. Without it, script 0.2 samples new QTN effects; this is a new simulation, not exact reproduction of the study. |
| `data/phenotype/0.3.seed.rds` | Optional historical scenario/replicate/crop seeds. Script 0.3 creates a new seed file if absent. |
| `data/selTrainData_GP/0.4.trainData_selCombName_allComb.rds` | Optional historical ten Design-D sets (26 × 20 = 520 unique pairs). Script 0.4 generates sets if absent. |
| `data/selTrainData_GP/0.4.trainData_selCombName_useOnce.rds` | Optional historical ten Design-B sets (520 unique pairs, each crop's marginal counts differing by at most one). Script 0.4 generates sets if absent. |
| `data/publicationResults.rds` | Required for default archive-based plotting/tables: parameters, crop-level predictive ability/timings, covariance recovery, BO discovery counts/timings and padded best-value trajectories. Unnecessary when rebuilding summaries from completed analysis outputs. |
| `data/2.1.publicationSelection.rds` | Required for default archive-based BO plotting: first-fit selections, candidate mean/SD and discovery-rank records. Unnecessary when rebuilding from completed BO fits. |

The historical genotype-effect archive, when provided locally, substitutes exact genotype-level products for large marker-effect matrices. The original `phenoNo` indexing is retained, including reuse of some genetic-effect realizations across scenario/replicate combinations. The local historical result archive preserves correlations rounded to three decimals by the evaluation script; it does not substitute new fits or higher-precision correlations. Individual MCMC draws and full GP fits are not distributed.

### Genotype preprocessing provenance

The original preprocessing used `gaston::select.snps(..., maf >= 0.1)`, followed by `gaston::LD.thin(..., threshold = 0.6)` on the **full** 237/198 panels, then removal of markers containing heterozygotes. No window or thinning-order override was passed. After accession selection, the MAF filter was applied again. Script 0.1 now starts from the already public filtered matrices and regenerates selected panels and marker metadata, avoiding a dependency on unavailable raw VCF files. Physical positions are converted using the original constant rates of 2.2 cM/Mb (sorghum) and 2.5 cM/Mb (soybean).

The LD threshold was chosen pragmatically from the retained SNP count, not by a documented optimization or sensitivity study. The original k-medoids call, distance definition and initialization seed were not retained in the available scripts; the actual accession IDs are retained only locally. This code-only release does not reconstruct the original medoid selection. Providing different 79/66 IDs defines a different experiment, not reproduction of the submitted results.

## Full analysis after supplying the required inputs

First supply `data/SelectedIndList_kMedoids_nIndSb_79_nIndGm_66.rds` in the format above. Use a **fresh checkout/output directory** for a new analysis. Existing scripts reuse completed files, and simply changing a seed in a directory containing old results will not refit those files. The full analysis is computationally intensive; it is separate from the default figure-redrawing workflow.

```sh
Rscript scripts/0.1.modifyGeno.r
Rscript scripts/0.2.sampleTrueMrkEff.r
Rscript scripts/0.3.generatePheno.r
Rscript scripts/0.4.selTrainDataGP.r
Rscript scripts/1.0.GPFitting.r
Rscript scripts/1.1.GPPred.r
Rscript -e 'options(ICBO.rebuildSummaries=TRUE); source("scripts/1.2.GPAccuracy.r")'
Rscript scripts/2.0.BO.R
Rscript -e 'options(ICBO.rebuildSummaries=TRUE); source("scripts/2.1.BOEval.R")'
Rscript -e 'options(ICBO.rebuildSummaries=TRUE); source("scripts/3.0.makeTables.R")'
# Optional checks after fits and their parameter traces exist:
Rscript scripts/3.1.checkMCMC.R
```

If the optional genotype-effect archive is present locally, script 0.2 reuses it without resampling QTNs. Otherwise it samples new effects with the documented fallback seeds. New effects, seeds or training sets can change the results from those in the manuscript. Script 1.0 retains its historical filename but only prepares kernels; the former full-panel diagnostic fits are not part of the manuscript analysis. Script 1.1 uses 40 PSOCK workers. Script 2.0 uses 40 fork workers and should be launched from a Linux/macOS terminal with `Rscript`, not sourced in a Positron R session that prohibits forking. Worker counts are explicit script parameters; their original values have not been changed. Memory requirements increase with concurrent workers and dense 5,214 × 5,214 kernels.

Published timing measurements came from an AMD Ryzen Threadripper 3990X, 64 physical cores/128 threads, about 251.5 GiB RAM, Ubuntu 22.04.2. They are per-task elapsed times under parallel execution, not exclusive CPU costs or a guaranteed wall-clock runtime on another machine. Table S3 counts 28,800 unique GP fits (joint multivariate fits counted once) and 1,800 BO runs per method.

## Model, prior and RNG details

All fits use BGLR 1.1.4, one chain, 10,000 iterations, burn-in 2,000 and thinning 20 (400 retained draws). `model="RKHS"` with genomic relationship kernels is the **GBLUP implementation**; no nonlinear Gaussian kernel is introduced.

| Model | Components and implementation |
|---|---|
| 1 | Single-crop DGE; `BGLR()` with one GRM kernel |
| 2 | Single-crop DGE + IGE; `BGLR()` with the two crop GRM kernels |
| 3 | Single-crop DGE + IGE + SMA; `BGLR()` with two crop kernels and their Kronecker interaction kernel; used in BO |
| 4 | `Multitrait()` with two crop kernels; unrestricted 2 × 2 genetic covariance matrices |
| 5 | Model 4 plus the Kronecker interaction kernel with diagonal across-crop SMA covariance |

Models 1, 2 and 4 keep the full 5,214-pair response structure with unobserved responses masked; Models 3 and 5 fit the 520 observed pairs and predict the other 4,694. Each BO fit uses all pairs observed up to that cycle. These implementations and priors were retained, rather than refactored into a different fitting algorithm.

Prior settings follow the installed BGLR 1.1.4 source. For univariate fits, global `R2=0.5`, variance-component `df0=5`, and automatic `S0` are used. Residual `S0=var(y,na.rm=TRUE)*(1-R2)*(df0+2)`. For an RKHS term, `R2_term=R2/nLT` and `S0=var(y,na.rm=TRUE)*R2_term/mean(d)*(df0+2)`, where `d` contains the retained kernel eigenvalues. Thus prior scales depend on the observed data and kernel, rather than being one fixed numerical constant across fits.

For multivariate fits, `R2=0.5`; genetic `UN` covariance priors use `df0=traits+1=3` and `S0=mo*(df0+traits+1)`, with `mo=(R2/nLT)*Sy/MSx`. Here `Sy=cov(y,use="pairwise.complete.obs")` and `MSx` is BGLR's centered kernel-design sum of squares per row. Model 5's `DIAG` SMA term uses `df0=c(3,3)` and `S0=(R2/nLT)*diag(Sy)/MSx*(df0+2)`. Residual covariance is explicitly `resCov=list(S0=diag(2),type="UN")`, with BGLR's default `df0=5`. The fitted objects retain the resulting hyperparameters.

Simulation heritability is `var(DGE+IGE+SMA)/(var(DGE+IGE+SMA)+var(error))`; `gamma=var(SMA)/var(DGE+IGE)`. The residual variance uses `var(GMA+SMA)`, including realized covariance. Original scaling of phenotypes and objective values is retained.

The BO objective is the sum of the two standardized latent genomic values. EI/PI use the maximum **current latent posterior mean over observed pairs**, computed from the same fit as the candidates. Candidate SD is a conditional latent genetic SD obtained by plugging in estimated variance components; it excludes future-observation residual variance and does not integrate over hyperparameter uncertainty. Across-crop variances are added because Model 3 fits the crops independently. One common incumbent is used for the entire batch. The true optimum is used only for initialization exclusion, simulation stopping and retrospective evaluation, not as the EI/PI incumbent.

| Random process | Reproducibility record |
|---|---|
| Historical QTN effects | Original RNG states unavailable; exact used genotype effects are retained locally, not distributed. Script 0.2's fallback seed base `20261106` defines **new** QTN realizations when the local archive is absent. |
| SMA and residuals | Historical per-scenario/crop seeds are retained locally in `0.3.seed.rds`, not distributed; `20261106` is the fallback for a new seed file. |
| Design-D | Original ten sets retained locally, not distributed; fallback seed `20261106` applies to newly generated sets. |
| Design-B | Seed `20260928 + trainingRep`; historical balanced sets retained locally, not distributed. |
| GP MCMC refits | Historical PSOCK-worker RNG states unavailable. New deterministic per-fit seeds: `20261106 + model*100000 + design*10000 + trainingRep*1000 + parNo*2 + cropNo`; joint models omit `cropNo`. These are **new refit seeds**, not historical seeds. |
| BO | Initial-set seed vector generated from `set.seed(20260928)`; task seed `seed[boRep]+1000*parNo`; checkpoint RNG states retained for resume. |
| Figure 3 intervals | 10,000 crossed-factor bootstrap draws, seed `20261101`; pairing preserved across methods and cycles. |

The locally retained numerical archives reproduce the submitted statistical summaries, but are not included in this code-only release. A fresh stochastic refit is not claimed to recover the historical GP MCMC results bit for bit, because those worker RNG states were not saved. Numerical libraries can also affect eigendecompositions and stochastic trajectories. Running means in the optional response checks illustrate stability in selected fits; they do not establish convergence of every chain.

## Statistical summaries and evaluation conventions

- Predictive ability is Pearson correlation with pseudo-phenotype. Means/tests use Fisher z; crop and training-set results are averaged within each phenotype replicate. GP comparisons use ten phenotype blocks, Friedman tests and exact paired Wilcoxon signed-rank tests, Holm-adjusted over ten model contrasts **within** each scenario/design/panel. No across-scenario adjustment is implied by the letters.
- Table S1 reports the back-transformed mean z and delta-method SE across ten phenotype means. Table S2 averages ten initial sets within phenotype, then reports mean/SE over ten phenotype blocks.
- Figure 3 uses all 100 paired runs per scenario and a crossed phenotype × initial-set bootstrap. It carries the attained optimum forward after stopping. Intervals are pointwise, not simultaneous; cycles represent equal field-evaluation budgets, not equal CPU time.
- Figure S9b uses pointwise t intervals (9 df) for paired cycle savings after averaging initial sets within phenotype. Figure 5b and S10 use SE across the ten phenotype summaries.
- The initial 56 pairs are cycle 1 and exclude the unique optimum. For Random, its rank among the remaining 5,158 pairs is uniform: `P(T<=c)=min(56*(c-1)/5158,1)` for integer `c>=1`. The final batch has six pairs; the expected discovery cycle is `1+mean(ceiling((1:5158)/56)) = 47.554091`, rounded to **47.554**. Random is analytical in this release, not a new Monte Carlo method mixed into BO.

## Release validation

Validation used the authors' **local inputs and saved results**, which are not all distributed. It is not evidence that a fresh public checkout contains all required data. The cleanup was checked without rerunning the expensive GP/BO simulations: all R files parsed; eight regenerated genotype/marker objects matched the existing inputs; GP/BO figure scripts and table scripts completed; all 180 Table S1 and 72 Table S2 mean/SE entries matched the previous numerical tables to numerical precision; Figure 3 means and bootstrap intervals were unchanged. Saved-fit reconstructions reproduced the selected first batches for Figures 4/5b and the discovery selections for all 7,200 runs. Optional running-mean plots were regenerated from existing traces. The submitted Word documents and original fitting results were not edited.

Code-only packaging was also checked: missing local archives produce explicit instructions, existing root-level genotype files are accepted by script 0.1, and regenerated table summaries can be used without the historical publication archive.

## Citation

Kinoshita and Iwata (2026). See the manuscript/preprint associated with this repository: [doi:10.64898/2026.05.15.725387](https://doi.org/10.64898/2026.05.15.725387). Figure numbering in this release follows the submitted revised version, which can differ from earlier preprint versions.
