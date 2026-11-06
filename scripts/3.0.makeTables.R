###### Supplementary Tables S1-S3 ######
rm(list = ls())
library(dplyr)
if (isTRUE(getOption("ICBO.rebuildSummaries", FALSE))) {
  gp <- readRDS("midstream/1.2.GP_Accuracy/1.2.GP_PublicationSummary.rds")
  bo <- readRDS("midstream/2.1.BO_Evaluation/2.1.BO_PublicationSummary.rds")
  parameters <- readRDS("data/phenotype/0.3.parmeterSet.rds")
} else {
  if (!file.exists("data/publicationResults.rds")) {
    stop("Publication archives are not distributed with this code-only release. ",
      "Supply the local archive, or rebuild the GP/BO summaries and set ",
      "options(ICBO.rebuildSummaries = TRUE); see README.md.", call. = FALSE)
  }
  publication <- readRDS("data/publicationResults.rds")
  gp <- publication$gpAccuracy
  bo <- publication$boCounts
  parameters <- publication$parameters
}
dir.create("midstream", showWarnings = FALSE)
designLabels <- c(All = "Design-D", Once = "Design-B")

### Table S1: average Fisher z within phenotype replicate, then across replicates.
gpBlocks <- gp %>% mutate(z = atanh(pmax(pmin(Accuracy, 1 - 1e-7), -1 + 1e-7))) %>%
  group_by(genCor, h2, gamma, CVMethod, Model, rep) %>%
  summarise(z = mean(z), nTrainingRep = n_distinct(CV), nCrops = n_distinct(Crop), .groups = "drop")
tableS1 <- gpBlocks %>% group_by(genCor, h2, gamma, CVMethod, Model) %>%
  summarise(meanZ = mean(z), seZ = sd(z) / sqrt(n()), nPhenotypeRep = n(),
    nTrainingRep = first(nTrainingRep), nCrops = first(nCrops), .groups = "drop") %>%
  mutate(Design = unname(designLabels[CVMethod]), Mean = tanh(meanZ),
    SE = (1 - Mean^2) * seZ) %>%
  select(genCor, h2, gamma, Design, Model, Mean, SE, nPhenotypeRep, nTrainingRep, nCrops)
write.csv(tableS1, "midstream/TableS1_PredictiveAbility.csv", row.names = FALSE)

### Table S2: average the ten initial sets within each phenotype replicate.
boBlocks <- bo %>% mutate(phenotypeRep = parameters$rep[ParNo]) %>%
  group_by(GenCor, h2, Gamma, Index, phenotypeRep) %>%
  summarise(cycles = mean(nIter), nInitialSetRep = n(), .groups = "drop")
tableS2 <- boBlocks %>% group_by(GenCor, h2, Gamma, Index) %>%
  summarise(Mean = mean(cycles), SE = sd(cycles) / sqrt(n()),
    nPhenotypeRep = n(), nInitialSetRep = first(nInitialSetRep), .groups = "drop")
write.csv(tableS2, "midstream/TableS2_BOCycles.csv", row.names = FALSE)

### Table S3: joint multivariate fits are counted once, not once per crop.
uniqueFits <- gp %>% distinct(genCor, h2, gamma, rep, CVMethod, CV, Model, FitCrop, .keep_all = TRUE)
tableS3GP <- uniqueFits %>% group_by(Model, CVMethod) %>%
  summarise(MeanSeconds = mean(fitPredictionTime), MedianSeconds = median(fitPredictionTime),
    SDSeconds = sd(fitPredictionTime), Q1Seconds = quantile(fitPredictionTime, .25),
    Q3Seconds = quantile(fitPredictionTime, .75), nFits = n(), .groups = "drop") %>%
  mutate(Design = unname(designLabels[CVMethod])) %>% select(-CVMethod)
tableS3BO <- attr(bo, "runTiming") %>% filter(method %in% c("Mu", "PI", "EI", "UCB")) %>%
  group_by(method) %>% summarise(MeanSeconds = mean(totalOptimizationTime),
    MedianSeconds = median(totalOptimizationTime), SDSeconds = sd(totalOptimizationTime),
    Q1Seconds = quantile(totalOptimizationTime, .25), Q3Seconds = quantile(totalOptimizationTime, .75),
    nRuns = n(), .groups = "drop")
stopifnot(nrow(tableS1) == 180L, nrow(tableS2) == 72L, sum(tableS3GP$nFits) == 28800L,
  all(tableS1$nPhenotypeRep == 10L), all(tableS2$nPhenotypeRep == 10L),
  all(tableS3BO$nRuns == 1800L))
write.csv(tableS3GP, "midstream/TableS3_GPTime.csv", row.names = FALSE)
write.csv(tableS3BO, "midstream/TableS3_BOTime.csv", row.names = FALSE)
