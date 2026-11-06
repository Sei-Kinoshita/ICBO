###### Optional running-mean checks reported in the reviewer response ######
### Run after new GP and BO fits; the public archive does not include all MCMC draws.
### These two illustrative checks are not additional manuscript figures.
rm(list = ls())
library(dplyr)
library(ggplot2)
source("scripts/function_scripts/2.1.function.BOEval.R")
source("scripts/function_scripts/1.2.function.GPAccuracy.R")
parameters <- readRDS("data/phenotype/0.3.parmeterSet.rds")
modelDirs <- c("Model1_DGE", "Model2_DGE_IGE", "Model3_GMA_SMA",
  "Model4_Multi_DGE_IGE", "Model5_Multi_GMA_SMA")
gp <- EvaluateGpRunningMeans(parameters, modelDirs, "midstream/1.1.GP_CV/")
gp$plotData <- filter(gp$plotData, Model == "Model3", trainingDesign == "All")
PlotMcmcRunningMeans(gp, "midstream/ResponseFigure1_GP_RunningMeans.pdf", "GP")
bo <- EvaluateBoRunningMeans(parameters, "midstream/2.0.BO/Model3_GMA_SMA/", cycles = c(1L, 3L))
bo$plotData <- filter(bo$plotData, method == "EI", h2 == .7, gamma == .2)
PlotMcmcRunningMeans(bo, "midstream/ResponseFigure2_BO_RunningMeans.pdf", "BO")
