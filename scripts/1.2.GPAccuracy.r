###### 0. Settings ######
##### 0.1. Remove objects #####
rm(list = ls())
gc(reset = TRUE);gc(reset = TRUE)

##### 0.2. Load Packages and Functions #####
require(ggplot2)
require(scales)
require(dplyr)
require(multcompView)

##### 0.3. Define Parameters #####
scriptID <- "1.2."

cropName <- c("Sb", "Gm")

nRep <- 10 ### repetition for phenotype
nCV <- 10
modelName <- c("Model1_DGE", "Model2_DGE_IGE", "Model3_GMA_SMA", 
               "Model4_Multi_DGE_IGE", "Model5_Multi_GMA_SMA")

##### 0.4. Define Directories #####
dirData <- "data/"
dirPheno <- paste0(dirData, "phenotype/")
dirComb <- paste0(dirData, "selTrainData_GP/")
dirGpPred <- "midstream/1.1.GP_CV/"

saveDir <- paste0("midstream/", scriptID, "GP_Accuracy/")
if (!dir.exists(saveDir)) {
  dir.create(saveDir, recursive = TRUE)
}


###### 1. Load Data ######
##### 1.1. Parameter Set #####
fileNameParSet <- paste0(dirPheno, "0.3.parmeterSet.rds")
parSetDf <- readRDS(file = fileNameParSet)
nPar <- nrow(parSetDf)

##### 1.2. Pair Name (fixed order) #####
fileNamePairName <- paste0(dirPheno, "0.3.pairNameDf_5214Pair.rds")
pairNameDf <- readRDS(file = fileNamePairName)

##### 1.3. Selected combination for CV #####
fileNameSelCombLst1 <- paste0(dirComb, "0.4.trainData_selCombName_random.rds")
fileNameSelCombLst2 <- paste0(dirComb, "0.4.trainData_selCombName_allComb.rds")
fileNameSelCombLst3 <- paste0(dirComb, "0.4.trainData_selCombName_useOnce.rds")
selCombRandomLst <- readRDS(file = fileNameSelCombLst1)
selCombAllLst <- readRDS(file = fileNameSelCombLst2)
selCombOnceLst <- readRDS(file = fileNameSelCombLst3)

cvCombLst <- list(Random = selCombRandomLst, All = selCombAllLst, Once = selCombOnceLst)
cvMethod <- names(cvCombLst)
rm(selCombRandomLst);rm(selCombAllLst);rm(selCombOnceLst)




###### 2. Accuracy Plot ######
##### 2.1. Scatter Plot #####
### for saving accuracy result
# addDf <- expand.grid(Crop = cropName, Model = modelName, 
#                      CVMethod = cvMethod, CV = 1:nCV)
# resDf <- cbind(parSetDf[rep(1:nPar, times = nrow(addDf)), ], 
#                addDf[rep(1:nrow(addDf), each = nPar), ])
# resDf$Accuracy <- rep(NA, nrow(resDf))
fileNameAccuDf <- paste0(saveDir, scriptID, "GP_Accuracy_df.rds")
if (!file.exists(fileNameAccuDf)) {
  ### plot and calculate accuracy
  resDf <- do.call(what = rbind, 
                   args = pbmcapply::pbmclapply(X = 1:nPar, FUN = function(parNo) {
                     
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
                                                                                                 "/allRes/1.1.Multitrait_predicted_yNew_parNo_", parNo, 
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
                                                                                                 "/allRes/1.1.Multitrait_fitting_allRes_parNo_", parNo, 
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
                                                                       ### accuracy plot
                                                                       plotDir <- paste0(saveDir, "GP_Accuracy_ScatterPlot/")
                                                                       if (!dir.exists(plotDir)) {
                                                                         dir.create(path = plotDir, recursive = TRUE)
                                                                       }
                                                                       fileNameScatPlot <- paste0(plotDir, "genCor_", genCor, "_h2_", h2, 
                                                                                                  "_gamma_", gamma, "_repNo_", repNo, "_", 
                                                                                                  cvMethodNow, "_", modelNow, "_", 
                                                                                                  cropNow, "_cvNo_", cvNo, ".png")
                                                                       if (!file.exists(fileNameScatPlot)) {
                                                                         png(filename = fileNameScatPlot, width = 1000, height = 1000)
                                                                         plotTitle <- paste0(modelNow, "/", cropNow, "/", cvMethodNow, 
                                                                                             "/CV:", cvNo, "/genCor:", genCor, 
                                                                                             "/h2:", h2, "/gamma:", gamma, "/repNo:", repNo)
                                                                         plot(yTrue, yPred, main = plotTitle, 
                                                                              xlab = "Observed", ylab = "Predicted")
                                                                         legend("bottomright", legend = paste0("r = ", phenoCor))
                                                                         dev.off()
                                                                       }
                                                                       
                                                                       ### result (correlation)
                                                                       resVec <- data.frame(genCor = genCor, h2 = h2, gamma = gamma, rep = repNo, 
                                                                                            Crop = cropNow, Model = modelNow, CVMethod = cvMethodNow, 
                                                                                            CV = cvNo, Accuracy = phenoCor)
                                                                       
                                                                       return(resVec)
                                                                     }))
                                                           }))
                                                 }))
                                       }))
                     
                     return(resMat)
                     
                   }, mc.cores = 40))
  
  resDf$genCor <- factor(resDf$genCor)
  resDf$h2 <- factor(resDf$h2)
  resDf$gamma <- factor(resDf$gamma)
  ### save
  saveRDS(object = resDf, file = fileNameAccuDf)
} else {
  resDf <- readRDS(file = fileNameAccuDf)
}




##### 2.2. Box Plot #####
#### 2.2.1. for each scenario ####
fileNameAccuBoxPlot <- paste0(saveDir, scriptID, "GP_Accuracy_BoxPlot.pdf")
pdf(file = fileNameAccuBoxPlot, width = 12, height = 10)
lapply(X = c("genCor", "h2", "gamma", "CVMethod"), FUN = function(xNow) {
  ggplot(resDf, 
         aes(x = .data[[xNow]], y = .data[["Accuracy"]], fill = .data[["Model"]])) + 
    geom_boxplot() + 
    facet_wrap(~Crop) + 
    # scale_fill_discrete(guide = "none") + 
    theme(text = element_text(size = 17))
})
dev.off()


#### 2.2.2. 9 cells (1 plot for h2=0.3/0.7, SB, GM) ####
### remove Random CV
resFilDf <- resDf[(resDf$CVMethod != "Random"), ]

resLowSb <- resFilDf[(resFilDf$h2 == 0.3) & (resFilDf$Crop == "Sb"), ]
resLowGm <- resFilDf[(resFilDf$h2 == 0.3) & (resFilDf$Crop == "Gm"), ]
resHighSb <- resFilDf[(resFilDf$h2 == 0.7) & (resFilDf$Crop == "Sb"), ]
resHighGm <- resFilDf[(resFilDf$h2 == 0.7) & (resFilDf$Crop == "Gm"), ]
resLow <- resFilDf[(resFilDf$h2 == 0.3), ]
resHigh <- resFilDf[(resFilDf$h2 == 0.7), ]

### plot
source("scripts/function_scripts/1.2.function.GPAccuracy.R")
fileNameAccu4Data <- paste0(saveDir, scriptID, "Accuracy_Model_Boxplot_Tukey.pdf")
pdf(fileNameAccu4Data, width = 12, height = 9)
p1 <- MakePlot(dat = resLowSb, plotTitle = "h2=0.3 / Sb")
p2 <- MakePlot(dat = resLowGm, plotTitle = "h2=0.3 / Gm")
p3 <- MakePlot(dat = resHighSb, plotTitle = "h2=0.7 / Sb")
p4 <- MakePlot(dat = resHighGm, plotTitle = "h2=0.7 / Gm")
p5 <- MakePlot(dat = resLow, plotTitle = "h2=0.3")
p6 <- MakePlot(dat = resHigh, plotTitle = "h2=0.7")
print(p1); print(p2); print(p5); print(p3); print(p4); print(p6)
dev.off()



# ### result only for Model3
# resM3Df <- resDf[resDf$Model == "Model3_GMA_SMA", ]
# resM3SbDf <- resM3Df[resM3Df$Crop == "Sb", ]
# resM3GmDf <- resM3Df[resM3Df$Crop == "Gm", ]
# 
# # resDfNow <- resM3SbDf ### Sorghum
# resDfNow <- resM3GmDf ### Soybean
# ggplot(data = resDfNow, aes(x = gamma, y = Accuracy, fill = genCor)) + 
#   stat_summary(fun = mean, 
#                geom = "bar", 
#                position = position_dodge(width = 0.9)) + 
#   stat_summary(fun.data = mean_se, 
#                geom = "errorbar", 
#                position = position_dodge(width = 0.9), 
#                width = 0.2) + 
#   facet_wrap(~h2) + 
#   scale_y_continuous(limits = c(0, 1)) + 
#   theme_bw() + 
#   theme(text = element_text(size = 17))
# 
# # mean(resDfNow$Accuracy[resDfNow$h2 == 0.3]) ### mean accuracy for h2=0.3
# # mean(resDfNow$Accuracy[resDfNow$h2 == 0.7]) ### mean accuracy for h2=0.7
