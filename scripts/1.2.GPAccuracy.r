###### 1. Settings ######
rm(list = ls())
require(dplyr)
require(multcompView)
source("scripts/function_scripts/2.1.function.BOEval.R")
source("scripts/function_scripts/1.2.function.GPAccuracy.R")
scriptID <- "1.2."
saveDir <- "midstream/1.2.GP_Accuracy/"
dir.create(saveDir, recursive = TRUE, showWarnings = FALSE)
fileNameAccuDf <- paste0(saveDir, scriptID, "GP_PublicationSummary.rds")
###### 2. Load prediction results ######
### Set ICBO.rebuildSummaries = FALSE to use a local archive.
useArchivedResults <- !isTRUE(getOption("ICBO.rebuildSummaries", TRUE))
if (useArchivedResults) {
  archiveFiles <- "data/publicationResults.rds"
  if (!all(file.exists(archiveFiles))) {
    stop("Missing local archive: data/publicationResults.rds", call. = FALSE)
  }
  publication <- readRDS("data/publicationResults.rds")
  resDf <- publication$gpAccuracy
  correlationRecovery <- publication$correlationRecovery
} else {
  ### Read the fits produced by 1.1.GPPred.r.
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
  ### Calculate predictive ability
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

###### 3. Paired model comparisons ######
figureGroups <- list(Figure2 = c("Sb", "Gm"), FigureS7 = "Sb", FigureS8 = "Gm")
comparison <- list()
for (figure in names(figureGroups)) {
  for (h in c(.3, .7)) {
    dat <- resDf[resDf$Crop %in% figureGroups[[figure]] & resDf$h2 == h, ]
    dat$CVMethod <- factor(dat$CVMethod, levels = c("All", "Once"), labels = c("CV1", "CV2"))
    dat$Model <- factor(dat$Model)
    dat$genCor <- factor(dat$genCor)
    dat$gamma <- factor(dat$gamma)
    lettersDf <- MakeLettersDf(dat)
    comparison[[length(comparison) + 1L]] <- attr(lettersDf, "modelComparison")
  }
}
###### 4. Save summaries ######
attr(resDf, "geneticCorrelationRecovery") <- correlationRecovery
attr(resDf, "modelComparison") <- bind_rows(comparison)
saveRDS(resDf, fileNameAccuDf)
