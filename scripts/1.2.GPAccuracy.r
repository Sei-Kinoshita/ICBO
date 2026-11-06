###### 1.2. Manuscript Figure 2 and Figures S3, S5-S8 ######
rm(list = ls())
require(ggplot2)
require(dplyr)
require(multcompView)
source("scripts/function_scripts/2.1.function.BOEval.R")
source("scripts/function_scripts/1.2.function.GPAccuracy.R")
scriptID <- "1.2."
saveDir <- "midstream/1.2.GP_Accuracy/"
dir.create(saveDir, recursive = TRUE, showWarnings = FALSE)
fileNameAccuDf <- paste0(saveDir, scriptID, "GP_PublicationSummary.rds")
### Default: redraw from locally supplied archives without fitting models.
useArchivedResults <- !isTRUE(getOption("ICBO.rebuildSummaries", FALSE))
if (useArchivedResults) {
  archiveFiles <- "data/publicationResults.rds"
  if (!all(file.exists(archiveFiles))) {
    stop("Publication archives are not distributed with this code-only release. ",
      "Supply the local archives, or rerun the analysis and set ",
      "options(ICBO.rebuildSummaries = TRUE); see README.md.", call. = FALSE)
  }
  publication <- readRDS("data/publicationResults.rds")
  resDf <- publication$gpAccuracy
  correlationRecovery <- publication$correlationRecovery
} else {
  ### Rebuild after running 0.1-1.1; new MCMC fits need not equal historical fits.
  nRep <- nCV <- 10L
  cropName <- c("Sb", "Gm")
  modelName <- c("Model1_DGE", "Model2_DGE_IGE", "Model3_GMA_SMA",
    "Model4_Multi_DGE_IGE", "Model5_Multi_GMA_SMA")
  dirPheno <- "data/phenotype/"
  dirGpPred <- "midstream/1.1.GP_CV/"
  parSetDf <- readRDS(paste0(dirPheno, "0.3.parmeterSet.rds"))
  parNoVec <- seq_len(nrow(parSetDf))
  pairNameDf <- readRDS(paste0(dirPheno, "0.3.pairNameDf_5214Pair.rds"))
  cvCombLst <- list(All = readRDS("data/selTrainData_GP/0.4.trainData_selCombName_allComb.rds"),
    Once = readRDS("data/selTrainData_GP/0.4.trainData_selCombName_useOnce.rds"))
  cvMethod <- names(cvCombLst)
  ### plot and calculate accuracy
  resDf <- do.call(what = rbind, 
                   args = lapply(X = parNoVec, FUN = function(parNo) {
                     
                     ### parmeter set
                     genCor <- parSetDf$genCor[parNo]
                     h2 <- parSetDf$h2[parNo]
                     gamma <- parSetDf$gamma[parNo]
                     repNo <- parSetDf$rep[parNo]
                     ### True phenotype
                     fileNamePheno <- paste0(dirPheno, "0.3.scaled_phenotype_genCor_", genCor,
                                             "_h2_", h2, "_gamma_", gamma, "_nRep_", nRep, ".rds")
                     phenoLst <- readRDS(fileNamePheno)
                     
                     resMat <- do.call(what = rbind, 
                                       args = lapply(X = cvMethod, FUN = function(cvMethodNow) {
                                         
                                         do.call(what = rbind, 
                                                 args = lapply(X = modelName, FUN = function(modelNow) {
                                                   
                                                   do.call(what = rbind, 
                                                           args = lapply(X = 1:nCV, FUN = function(cvNo) {
                                                             
                                                             do.call(what = rbind, 
                                                                     args = lapply(X = cropName, FUN = function(cropNow) {
                                                                       
                                                                       cropNo <- which(cropName %in% cropNow)
                                                                       ### all individuals' phenotype
                                                                       phenoNow <- phenoLst[[repNo]][[cropNow]]
                                                                       
                                                                       ### GP CV result
                                                                       if (modelNow == "Model5_Multi_GMA_SMA") {
                                                                         ### for Model 5
                                                                         fileNameCvRes <- paste0(dirGpPred, modelNow, "/", cvMethodNow, "/CV", cvNo, 
                                                                                                 "/allRes/1.1.", "Multitrait_predicted_yNew_parNo_", parNo, 
                                                                                                 "_", cvMethodNow, "_CV_", cvNo, ".rds")
                                                                         ### predicted phenotype
                                                                         gpResNow <- readRDS(file = fileNameCvRes)
                                                                         yPred <- gpResNow[, cropNow]
                                                                         ### trained data name
                                                                         cvCombLstNow <- cvCombLst[[cvMethodNow]][[cvNo]]
                                                                         trainPos <- match(paste(cvCombLstNow$Sb, cvCombLstNow$Gm),
                                                                                           paste(pairNameDf$Sb, pairNameDf$Gm))
                                                                         ### true phenotype
                                                                         yTrue <- phenoNow[-trainPos]
                                                                         
                                                                       } else if (modelNow == "Model4_Multi_DGE_IGE") {
                                                                         ### for Model 4
                                                                         fileNameCvRes <- paste0(dirGpPred, modelNow, "/", cvMethodNow, "/CV", cvNo, 
                                                                                                 "/allRes/1.1.", "Multitrait_fitting_allRes_parNo_", parNo, 
                                                                                                 "_", cvMethodNow, "_CV_", cvNo, ".rds")
                                                                         
                                                                         ### predicted phenotype
                                                                         gpResNow <- readRDS(file = fileNameCvRes)
                                                                         ### trained data name
                                                                         cvCombLstNow <- cvCombLst[[cvMethodNow]][[cvNo]]
                                                                         trainPos <- match(paste(cvCombLstNow$Sb, cvCombLstNow$Gm),
                                                                                           paste(pairNameDf$Sb, pairNameDf$Gm))
                                                                         yPred <- gpResNow$ETAHat[-trainPos, cropNo]
                                                                         ### true phenotype
                                                                         yTrue <- phenoNow[-trainPos]
                                                                         
                                                                       } else if (modelNow == "Model3_GMA_SMA") {
                                                                         ### for Model 3
                                                                         fileNameCvRes <- paste0(dirGpPred, modelNow, "/", cvMethodNow, "/CV", cvNo, 
                                                                                                 "/allRes/1.1.", cropNow, "_predicted_yNew_parNo_", parNo, 
                                                                                                 "_", cvMethodNow, "_CV_", cvNo, ".rds")
                                                                         ### predicted phenotype
                                                                         yPred <- readRDS(file = fileNameCvRes)
                                                                         ### trained data name
                                                                         cvCombLstNow <- cvCombLst[[cvMethodNow]][[cvNo]]
                                                                         trainPos <- match(paste(cvCombLstNow$Sb, cvCombLstNow$Gm),
                                                                                           paste(pairNameDf$Sb, pairNameDf$Gm))
                                                                         ### true phenotype
                                                                         yTrue <- phenoNow[-trainPos]
                                                                         
                                                                       } else {
                                                                         ### for Model 1 & 2
                                                                         fileNameCvRes <- paste0(dirGpPred, modelNow, "/", cvMethodNow, "/CV", cvNo, 
                                                                                                 "/allRes/1.1.", cropNow, "_fitting_allRes_parNo_", parNo, 
                                                                                                 "_", cvMethodNow, "_CV_", cvNo, ".rds")
                                                                         gpResNow <- readRDS(file = fileNameCvRes)
                                                                         ### predicted phenotype
                                                                         yPred <- gpResNow$yHat[is.na(gpResNow$y)]
                                                                         ### true phenotype
                                                                         yTrue <- phenoNow[is.na(gpResNow$y)]
                                                                       }
                                                                       ### accuracy
                                                                       phenoCor <- round(cor(yTrue, yPred), digits = 3)
                                                                       ### result (correlation)
                                                                       resVec <- data.frame(genCor = genCor, h2 = h2, gamma = gamma, rep = repNo, 
                                                                                            Crop = cropNow, Model = modelNow, CVMethod = cvMethodNow, 
                                                                                            CV = cvNo, Accuracy = phenoCor)
                                                                       
                                                                       timingNow <- if (modelNow == "Model3_GMA_SMA") {
                                                                         attr(yPred, "timing")
                                                                       } else if (modelNow == "Model5_Multi_GMA_SMA") {
                                                                         attr(gpResNow, "timing")
                                                                       } else gpResNow$timing
                                                                       resVec$fitPredictionTime <- timingNow$fitPredictionTime
                                                                       resVec$nTrainingPair <- timingNow$nTrainingPair
                                                                       resVec$nPredictionCandidate <- timingNow$nPredictionCandidate
                                                                       resVec$FitCrop <- timingNow$Crop
                                                                       return(resVec)
                                                                     }))
                                                           }))
                                                 }))
                                       }))
                     
                     return(resMat)
                     
                   }))
  
  resDf$genCor <- factor(resDf$genCor)
  resDf$h2 <- factor(resDf$h2)
  resDf$gamma <- factor(resDf$gamma)
  correlationRecovery <- EvaluateGeneticCorrelationRecovery(parSetDf, pairNameDf,
    dirPheno, dirGpPred, saveDir, scriptID, cvCombLst,
    nIter = 10000L, burnIn = 2000L, thin = 20L)
}

###### Published predictive-ability panels and their paired comparisons ######
figureGroups <- list(Figure2 = c("Sb", "Gm"), FigureS7 = "Sb", FigureS8 = "Gm")
comparison <- list()
for (figure in names(figureGroups)) {
  pdf(paste0(saveDir, scriptID, figure, "_PredictiveAbility.pdf"), width = 12, height = 9)
  for (h in c(.3, .7)) {
    dat <- resDf[resDf$Crop %in% figureGroups[[figure]] & resDf$h2 == h, ]
    p <- MakePlot(dat, paste0("Predictive ability / h2 = ", h))
    print(p)
    comparison[[length(comparison) + 1L]] <- attr(p, "modelComparison")
  }
  dev.off()
}
PlotGeneticCorrelationRecovery(correlationRecovery, saveDir, scriptID)
attr(resDf, "geneticCorrelationRecovery") <- correlationRecovery
attr(resDf, "modelComparison") <- bind_rows(comparison)
saveRDS(resDf, fileNameAccuDf)
