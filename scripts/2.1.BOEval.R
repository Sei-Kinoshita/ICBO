###### 1. Settings ######
rm(list = ls())
require(dplyr)
require(tidyr)
source("scripts/function_scripts/2.0.function.BO.R")
source("scripts/function_scripts/2.1.function.BOEval.R")
scriptID <- "2.1."
saveDir <- "midstream/2.1.BO_Evaluation/"
dir.create(saveDir, recursive = TRUE, showWarnings = FALSE)
indexName <- c("Mu", "PI", "EI", "UCB")
nIndex <- length(indexName)
nRep <- nRepBo <- 10L
nPair <- 5214L
kappa <- round(qnorm(.95), 2)
###### 2. Load optimization results ######
useArchivedResults <- !isTRUE(getOption("ICBO.rebuildSummaries", TRUE))
if (useArchivedResults) {
  archiveFiles <- c("data/publicationResults.rds", "data/2.1.publicationSelection.rds")
  if (!all(file.exists(archiveFiles))) {
    stop("Missing local BO archives: ", paste(archiveFiles[!file.exists(archiveFiles)], collapse = ", "),
      call. = FALSE)
  }
  publication <- readRDS("data/publicationResults.rds")
  parSetDf <- publication$parameters
  boResLst <- publication$boCurves
  boIterDf <- publication$boCounts
  selection <- readRDS("data/2.1.publicationSelection.rds")
} else {
  parSetDf <- readRDS("data/phenotype/0.3.parmeterSet.rds")
  boDir <- "midstream/2.0.BO/Model3_GMA_SMA/"
  rebuilt <- ReadPublicationBoResults(parSetDf, boDir, nRepBo)
  boResLst <- rebuilt$curves
  boIterDf <- rebuilt$counts
  selection <- ReadPublicationSelection(parSetDf, boDir, "data/phenotype/", nRep, nRepBo, kappa)
}
nPar <- nrow(parSetDf)
nParUniq <- nrow(unique(parSetDf[c("genCor", "h2", "gamma")]))

###### 3. Cycle-wise advantage and selection summaries ######
attr(boIterDf, "cycleAdvantage") <- EvaluateCycleAdvantage(boResLst, parSetDf, nRepBo)
attr(boIterDf, "selectionMechanism") <- selection

###### 4. Save summaries ######
saveRDS(boIterDf, paste0(saveDir, scriptID, "BO_PublicationSummary.rds"))
