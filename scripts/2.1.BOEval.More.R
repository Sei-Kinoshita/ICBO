###### 0. Settings ######
##### 0.1. Remove objects #####
rm(list = ls())
gc(reset = TRUE);gc(reset = TRUE)

##### 0.2. Load Packages and Functions #####
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(grid)
### function for plot
source("scripts/function_scripts/2.1.function.BOEval.R")

##### 0.3. Define Parameters #####
scriptID <- "2.1."

nRep <- 10 ### phenotype rep
nRepBo <- 10
kappa <- round(qnorm(p = 0.95), digits = 2)
# indexName <- c("Mu", "PI", "EI", paste0("UCB_", kappa))
indexName <- c("Mu", "PI", "EI", "UCB")
nIndex <- length(indexName)

minBoForPlot <- 20   ### no plot for BoNo<minBoForPlot
ciLevel <- 0.95      ### confidence interval for plot

##### 0.4. Define Directories #####
dataDir <- "data/"
phenoDir <- paste0(dataDir, "phenotype/")
boDir <- "midstream/2.0.BO/Model3_GMA_SMA/"

tarDir <- "midstream/2.1.BO_Evaluation/MoreInterpretation/"
saveDir <- paste0(tarDir, "Plot/")
if (!dir.exists(saveDir)) {
  dir.create(saveDir, recursive = TRUE)
}


###### 1. Load Data ######
### scenario parameter set
parSetDf <- readRDS(file = "data/phenotype/0.3.parmeterSet.rds")
nPar <- nrow(parSetDf)

parUniqVec <- (1:nPar) %% 18
parUniqVec[parUniqVec == 0] <- 18
nParUniq <- length(unique(parUniqVec))

### initial pair No
fileNameInitPariNo <- paste0(boDir, "2.0.BO_initPairNoLst.rds")
initPairNoLst <- readRDS(file = fileNameInitPariNo)



###### 2. Calculate Metrics from BO result ######
pbmcapply::pbmclapply(X = 1:nParUniq, FUN = function(parUniqNo) {
  
  genCorNow <- parSetDf$genCor[parUniqNo]
  gammaNow <- parSetDf$gamma[parUniqNo]
  h2Now <- parSetDf$h2[parUniqNo]
  # #############################
  # genCorNow <- 0.7; gammaNow <- 0.4; h2Now <- 0.3
  # #############################
  parNoVec <- which((parSetDf$genCor == genCorNow) & (parSetDf$gamma == gammaNow) 
                    & (parSetDf$h2 == h2Now))
  
  ### load phenotype to calc true genotypic values
  fileNamePheno <- paste0(phenoDir, "0.3.allSampledRes_phenotype_genCor_", genCorNow, 
                          "_h2_", h2Now, "_gamma_", gammaNow, "_nRep_", nRep, ".rds")
  ### yTrue
  yTrueLst <- lapply(X = 1:nRep, FUN = function(repNo) {
    phenoLst <- readRDS(file = fileNamePheno)[[repNo]]
    yTrue <- scale(phenoLst$yA) + scale(phenoLst$yB)
    return(yTrue)
  })
  ### gTrue
  gTrueLst <- lapply(X = 1:nRep, FUN = function(repNo) {
    phenoLst <- readRDS(file = fileNamePheno)[[repNo]]
    ### true best pair
    gTrueA <- scale(phenoLst$dA + phenoLst$sBtoA + phenoLst$vA)
    gTrueB <- scale(phenoLst$dB + phenoLst$sAtoB + phenoLst$vB)
    gTrue <- gTrueA + gTrueB
    return(gTrue)
  })
  
  ### calculate correlation for each acquisition function
  lapply(X = parNoVec, FUN = function(parNo) {
    ###
    listNo <- which(parNoVec == parNo)
    ###
    gTrue <- gTrueLst[[listNo]]
    yTrue <- yTrueLst[[listNo]]
    gTopTrue <- order(gTrue, decreasing = TRUE)[1:56]
    yTopTrue <- order(yTrue, decreasing = TRUE)[1:56]
    ### directory of each index's result
    parDir <- paste0(boDir, "parNo_", parNo, "/")
    
    ### for each BO repitition
    lapply(X = 1:nRepBo, FUN = function(boRepNo) {
      
      ### initial selected pairNo for BO
      initPairNo <- initPairNoLst[[boRepNo]]
      ### for each acquisition function
      indexLst <- lapply(X = c("Mu", "PI", "EI", "UCB"), FUN = function(indexNow) {
        
        boResDir <- paste0(parDir, indexNow, "/BO_result/")
        gpResDir <- paste0(parDir, indexNow, "/BGLR_result/")
        
        ### selected pairs' No
        fileNameSelPair <- paste0(boResDir, "2.0.selectedPairNo_allIteration_BORepNo_", 
                                  boRepNo, ".rds")
        selPairNoLst <- readRDS(file = fileNameSelPair)
        ### gBest for all iterations
        fileNameGBest <- paste0(boResDir, "2.0.g_best_allIteration_BORepNo_", 
                                boRepNo, ".rds")
        gBestLst <- readRDS(file = fileNameGBest)
        
        ### all iterations' gpRes files
        fileNameLstGpRes <- list.files(path = gpResDir, 
                                       pattern = paste0("2.0.fitting_result_BORepNo_", boRepNo, "_iterNo_"))
        nIterVec <- sort(as.numeric(sub(".*iterNo_([0-9]+)\\.rds$", "\\1", fileNameLstGpRes)))
        ### for each iteration
        fileNameSaveResEachIndex <- paste0(tarDir, scriptID, "allIterLst_parNo_", 
                                           parNo, "_boRepNo_", boRepNo, "_", 
                                           indexNow, ".rds")
        if (!file.exists(fileNameSaveResEachIndex)) {
          
          allResIterLst <- lapply(X = nIterVec, FUN = function(iterNo) {
            
            ### gpRes
            fileNameGpRes <- paste0(gpResDir, "2.0.fitting_result_BORepNo_", boRepNo, 
                                    "_iterNo_", iterNo, ".rds")
            gpResLst <- readRDS(file = fileNameGpRes)
            if (iterNo == 1) {
              selctedPairNoSum <- initPairNo
            } else {
              selctedPairNoSum <- initPairNo
              for (x in 2:iterNo) {
                selctedPairNoSum <- c(selctedPairNoSum, selPairNoLst[[x]])
              }
            }
            
            ### posterior for all pairs
            postMean <- gpResLst$Sb$mean_A_B_AB + gpResLst$Gm$mean_A_B_AB
            postSd <- gpResLst$Sb$sd_A_B_AB + gpResLst$Gm$sd_A_B_AB
            ### posterior for selected 56 pairs
            allNoVec <- rep(NA, length(gTrue))
            ### order in 5214 (initial)
            allNoVec[-selctedPairNoSum] <- seq_len(length(gTrue) - length(selctedPairNoSum))
            ### order in other cycles
            newSelPairNo <- allNoVec[selPairNoLst[[iterNo + 1]]]
            postMeanTop <- postMean[newSelPairNo]
            postSdTop <- postSd[newSelPairNo]
            ### top 56 for yTrue in the current cycle (order in new)
            gTrueTopNow <- order(gTrue[-selctedPairNoSum], decreasing = TRUE)[1:56]
            
            ###
            ### all pairs' posterior sd v.s. gTrue
            sdTrueCor <- cor(x = postSdTop,
                             y = gTrue[selPairNoLst[[iterNo + 1]]],
                             method = "spearman")
            ###
            ### all pairs' posterior sd v.s. (yTrue - posterior mean) ###
            sdErrorCor <- cor(x = postSdTop,
                              y = (abs(yTrue[selPairNoLst[[iterNo + 1]]] - postMeanTop)),
                              method = "spearman")
            ###
            ### for selected 56 pairs: posterior mean v.s. gTrue
            meanTrueCor <- cor(x = postMeanTop, 
                               y = gTrue[selPairNoLst[[iterNo + 1]]],
                               method = "spearman")
            ###
            ### for selected 56 pairs: how many top 56 pairs match for posterior mean and gTrue
            topMatch <- length(intersect(newSelPairNo, 
                                         gTrueTopNow))
            
            ### result for return
            corRes <- c(postSd_gTrue = sdTrueCor, 
                        postSd_error = sdErrorCor, 
                        postMean_gTrue = meanTrueCor, 
                        nTopMatch = topMatch)
            
            ### only for PI and EI
            if (indexNow %in% c("PI", "EI")) {
              z <- (postMeanTop - gBestLst[1]) / postSdTop
              denz <- dnorm(z)
              proz <- pnorm(z)
              if (indexNow == "PI") {
                ### dPI/dMu
                dbyMu <- denz / postSdTop
                dbySigma <- -denz * z / postSdTop
              } else {
                dbyMu <- proz
                dbySigma <- -denz
              }
              ###
              ### mu & sigma ration for top 56 pairs
              muSigmaRat <- dbyMu / dbySigma
              
              resLst <- list(Correlation = corRes, Mu_Sigma_Ratio = muSigmaRat)
              return(resLst)
            }
            
            return(corRes)
          })
          
          ### save
          if (length(allResIterLst) == 0) {
            saveRDS(object = paste0("Encounter Best Value at initial pairs"), 
                    file = fileNameSaveResEachIndex)
          } else {
            names(allResIterLst) <- paste0("Iter_", 1:length(allResIterLst))
            saveRDS(object = allResIterLst, file = fileNameSaveResEachIndex)
          }
        } else {
          allResIterLst <- readRDS(fileNameSaveResEachIndex)
        }
      })
    })
  })
}, mc.cores = nParUniq)



###### 3.1 Plot for each BORep ######
metricOrder <- c("postSd_gTrue", "postSd_error", "postMean_gTrue_top56", "nTop_Match")
### y-axis range for each metric across all nParUniq
allResDf_yRange <- do.call(what = rbind, args = lapply(X = 1:nParUniq, FUN = function(parUniqNo) {
  
  genCorNow <- parSetDf$genCor[parUniqNo]
  gammaNow <- parSetDf$gamma[parUniqNo]
  h2Now <- parSetDf$h2[parUniqNo]
  
  parNoVec <- which((parSetDf$genCor == genCorNow) & (parSetDf$gamma == gammaNow) 
                    & (parSetDf$h2 == h2Now))
  
  allResDf <- do.call(what = rbind, args = lapply(X = parNoVec, FUN = function(parNo) {
    
    genCorNow <- parSetDf$genCor[parNo]
    gammaNow <- parSetDf$gamma[parNo]
    h2Now <- parSetDf$h2[parNo]
    qtnRepNoNow <- parSetDf$rep[parNo]
    
    oneParDf <- do.call(what = rbind, args = lapply(X = 1:nRepBo, FUN = function(boRepNo) {
      
      oneBoRepDf <- do.call(what = rbind, args = lapply(X = indexName, FUN = function(indexNameNow) {
        
        tarFileNameNow <- paste0(tarDir, "2.1.allIterLst_parNo_", parNo, 
                                 "_boRepNo_", boRepNo, "_", indexNameNow, ".rds")
        inputNow <- readRDS(file = tarFileNameNow)
        
        oneIndexDf <- do.call(what = rbind, 
                              args = lapply(X = 1:length(inputNow), FUN = function(cycleNo) {
                                
                                corVec <- inputNow[[cycleNo]]
                                if (indexNameNow %in% c("PI", "EI")) { 
                                  corVec <- inputNow[[cycleNo]][[1]]
                                }
                                boNoNow <- ((qtnRepNoNow - 1) * nRepBo) + boRepNo
                                data.frame(GenCor = genCorNow, 
                                           Gamma = gammaNow, 
                                           h2 = h2Now, 
                                           QtnRep = qtnRepNoNow, 
                                           BoRep = boRepNo, 
                                           Index = indexNameNow, 
                                           BoNo = boNoNow, 
                                           Cycle = cycleNo, 
                                           postSd_gTrue = corVec[1], 
                                           postSd_error = corVec[2],
                                           postMean_gTrue_top56 = corVec[3],
                                           nTop_Match = corVec[4])
                              }))
        
        return(oneIndexDf)
      }))
      
      return(oneBoRepDf)
    }))
    
    return(oneParDf)
  }))
  
  rownames(allResDf) <- 1:nrow(allResDf)
  
  allResDf <- allResDf %>%
    mutate(
      Index = factor(Index, levels = indexName),
      Scenario = paste0("GenCor=", GenCor, ", Gamma=", Gamma, ", h2=", h2)
    )
  
  ### same summary as SaveBoMeanPlotPdf
  d_long <- allResDf %>%
    pivot_longer(
      cols = all_of(metricOrder),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    mutate(
      Metric = factor(Metric, levels = metricOrder),
      Value = as.numeric(Value)
    )
  
  alphaNow <- 1 - ciLevel
  
  d_sum <- d_long %>%
    group_by(Scenario, Index, Metric, Cycle) %>%
    summarise(
      nBo = n_distinct(BoNo),
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      SE = SD / sqrt(nBo),
      tcrit = qt(1 - alphaNow / 2, df = pmax(nBo - 1, 1)),
      Lower = Mean - tcrit * SE,
      Upper = Mean + tcrit * SE,
      .groups = "drop"
    ) %>%
    mutate(
      Lower = ifelse(is.na(Lower), Mean, Lower),
      Upper = ifelse(is.na(Upper), Mean, Upper)
    )
  
  d_sum_plot <- d_sum %>%
    filter(nBo >= minBoForPlot) %>%
    filter(!(Metric %in% c("postSd_gTrue", "postSd_error") & Index == "Mu"))
  
  return(d_sum_plot)
}))

metricYRangeDf <- allResDf_yRange %>%
  pivot_longer(
    cols = c(Mean, Lower, Upper),
    names_to = "Stat",
    values_to = "AxisValue"
  ) %>%
  group_by(Metric) %>%
  summarise(
    ymin = min(AxisValue, na.rm = TRUE),
    ymax = max(AxisValue, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Metric = factor(Metric, levels = metricOrder)
  )

### plot
pbmcapply::pbmclapply(X = 1:nParUniq, FUN = function(parUniqNo) {
  
  genCorNow <- parSetDf$genCor[parUniqNo]
  gammaNow <- parSetDf$gamma[parUniqNo]
  h2Now <- parSetDf$h2[parUniqNo]
  # #############################
  # genCorNow <- 0.7; gammaNow <- 0.4; h2Now <- 0.3
  # #############################
  parNoVec <- which((parSetDf$genCor == genCorNow) & (parSetDf$gamma == gammaNow) 
                    & (parSetDf$h2 == h2Now))
  
  ### each scenario and QTN rep
  allResDf <- do.call(what = rbind, args = lapply(X = parNoVec, FUN = function(parNo) {
    
    # parNo <- parNoVec[1]
    
    genCorNow <- parSetDf$genCor[parNo]
    gammaNow <- parSetDf$gamma[parNo]
    h2Now <- parSetDf$h2[parNo]
    qtnRepNoNow <- parSetDf$rep[parNo]
    
    ### each BO rep
    oneParDf <- do.call(what = rbind, args = lapply(X = 1:nRepBo, FUN = function(boRepNo) {
      
      # boRepNo <- 1
      
      ### each index
      oneBoRepDf <- do.call(what = rbind, args = lapply(X = indexName, FUN = function(indexNameNow) {
        
        # indexNameNow <- indexName[1]
        
        tarFileNameNow <- paste0(tarDir, "2.1.allIterLst_parNo_", parNo, 
                                 "_boRepNo_", boRepNo, "_", indexNameNow, ".rds")
        inputNow <- readRDS(file = tarFileNameNow)
        
        ### each cycle
        oneIndexDf <- do.call(what = rbind, 
                              args = lapply(X = 1:length(inputNow), FUN = function(cycleNo) {
                                
                                # cycleNo <- 1
                                corVec <- inputNow[[cycleNo]]
                                if (indexNameNow %in% c("PI", "EI")) { 
                                  corVec <- inputNow[[cycleNo]][[1]]
                                }
                                boNoNow <- ((qtnRepNoNow - 1) * nRepBo) + boRepNo
                                resDf <- data.frame(GenCor = genCorNow, 
                                                    Gamma = gammaNow, 
                                                    h2 = h2Now, 
                                                    QtnRep = qtnRepNoNow, 
                                                    BoRep = boRepNo, 
                                                    Index = indexNameNow, 
                                                    BoNo = boNoNow, 
                                                    Cycle = cycleNo, 
                                                    postSd_gTrue = corVec[1], 
                                                    postSd_error = corVec[2],
                                                    postMean_gTrue_top56 = corVec[3],
                                                    nTop_Match = corVec[4])
                                return(resDf)
                              }))
        
        return(oneIndexDf)
      }))
      
      return(oneBoRepDf)
    }))
    
    return(oneParDf)
  }))
  rownames(allResDf) <- 1:nrow(allResDf)
  
  ### for plot
  allResDf <- allResDf %>%
    mutate(
      Index = factor(Index, levels = indexName),
      Scenario = paste0("GenCor=", GenCor, ", Gamma=", Gamma, ", h2=", h2)
    )
  
  ### save and plot
  fileNamePlotNow <- paste0(saveDir, scriptID, "MetricChange_BoNo_GenCor_", 
                            genCorNow, "_Gamma_", gammaNow, "_h2_", h2Now, ".pdf")
  if (!file.exists(fileNamePlotNow)) {
    SaveBoPlotPdf(
      df = allResDf,
      pdf_file = fileNamePlotNow,
      GenCor_value = genCorNow,
      Gamma_value = gammaNow,
      h2_value = h2Now
    )
  }
  
  ### Plot to summarize BoRep and PhenoRep
  ### save and plot
  fileNamePlotNow <- paste0(saveDir, scriptID, "MetricChange_CI_GenCor_",
                            genCorNow, "_Gamma_", gammaNow, "_h2_", h2Now, ".pdf")
  if (!file.exists(fileNamePlotNow)) {
    SaveBoMeanPlotPdf(
      df = allResDf,
      pdf_file = fileNamePlotNow,
      GenCor_value = genCorNow,
      Gamma_value = gammaNow,
      h2_value = h2Now,
      minBoForPlot = minBoForPlot, 
      ciLevel = ciLevel,
      metricYRangeDf = metricYRangeDf
    )
  }
}, mc.cores = nParUniq)

##### 3.2. 9 cells plot #####
### data for 9-panel plot
allResDf_panel <- do.call(what = rbind, args = lapply(X = 1:nParUniq, FUN = function(parUniqNo) {
  
  genCorNow <- parSetDf$genCor[parUniqNo]
  gammaNow <- parSetDf$gamma[parUniqNo]
  h2Now <- parSetDf$h2[parUniqNo]
  
  parNoVec <- which((parSetDf$genCor == genCorNow) & (parSetDf$gamma == gammaNow) 
                    & (parSetDf$h2 == h2Now))
  
  allResDf <- do.call(what = rbind, args = lapply(X = parNoVec, FUN = function(parNo) {
    
    genCorNow <- parSetDf$genCor[parNo]
    gammaNow <- parSetDf$gamma[parNo]
    h2Now <- parSetDf$h2[parNo]
    qtnRepNoNow <- parSetDf$rep[parNo]
    
    oneParDf <- do.call(what = rbind, args = lapply(X = 1:nRepBo, FUN = function(boRepNo) {
      
      oneBoRepDf <- do.call(what = rbind, args = lapply(X = indexName, FUN = function(indexNameNow) {
        
        tarFileNameNow <- paste0(tarDir, "2.1.allIterLst_parNo_", parNo, 
                                 "_boRepNo_", boRepNo, "_", indexNameNow, ".rds")
        inputNow <- readRDS(file = tarFileNameNow)
        
        oneIndexDf <- do.call(what = rbind, 
                              args = lapply(X = 1:length(inputNow), FUN = function(cycleNo) {
                                
                                corVec <- inputNow[[cycleNo]]
                                if (indexNameNow %in% c("PI", "EI")) { 
                                  corVec <- inputNow[[cycleNo]][[1]]
                                }
                                boNoNow <- ((qtnRepNoNow - 1) * nRepBo) + boRepNo
                                data.frame(GenCor = genCorNow, 
                                           Gamma = gammaNow, 
                                           h2 = h2Now, 
                                           QtnRep = qtnRepNoNow, 
                                           BoRep = boRepNo, 
                                           Index = indexNameNow, 
                                           BoNo = boNoNow, 
                                           Cycle = cycleNo, 
                                           postSd_gTrue = corVec[1], 
                                           postSd_error = corVec[2],
                                           postMean_gTrue_top56 = corVec[3],
                                           nTop_Match = corVec[4])
                              }))
        
        return(oneIndexDf)
      }))
      
      return(oneBoRepDf)
    }))
    
    return(oneParDf)
  }))
  
  rownames(allResDf) <- 1:nrow(allResDf)
  
  allResDf <- allResDf %>%
    mutate(
      Index = factor(Index, levels = indexName),
      Scenario = paste0("GenCor=", GenCor, ", Gamma=", Gamma, ", h2=", h2)
    )
  
  return(allResDf)
}))

allResDf_panel_sum <- allResDf_panel %>%
  pivot_longer(
    cols = all_of(metricOrder),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(Metric, levels = metricOrder),
    Value = as.numeric(Value)
  )

alphaNow <- 1 - ciLevel

allResDf_panel_sum <- allResDf_panel_sum %>%
  group_by(GenCor, Gamma, h2, Scenario, Index, Metric, Cycle) %>%
  summarise(
    nBo = n_distinct(BoNo),
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    SE = SD / sqrt(nBo),
    tcrit = qt(1 - alphaNow / 2, df = pmax(nBo - 1, 1)),
    Lower = Mean - tcrit * SE,
    Upper = Mean + tcrit * SE,
    .groups = "drop"
  ) %>%
  mutate(
    Lower = ifelse(is.na(Lower), Mean, Lower),
    Upper = ifelse(is.na(Upper), Mean, Upper)
  ) %>%
  filter(nBo >= minBoForPlot) %>%
  filter(!(Metric %in% c("postSd_gTrue", "postSd_error") & Index == "Mu"))

metricYRangeDf_panel <- allResDf_panel_sum %>%
  pivot_longer(
    cols = c(Mean, Lower, Upper),
    names_to = "Stat",
    values_to = "AxisValue"
  ) %>%
  group_by(Metric) %>%
  summarise(
    ymin = min(AxisValue, na.rm = TRUE),
    ymax = max(AxisValue, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Metric = factor(Metric, levels = metricOrder)
  )

fileNamePlotNow <- paste0(saveDir, scriptID, "MetricChange_CI_9Panel_h2_0.3.pdf")
SaveBoMeanPlotPdf_9Panel(
  df = allResDf_panel,
  pdf_file = fileNamePlotNow,
  h2_value = 0.3,
  minBoForPlot = minBoForPlot,
  ciLevel = ciLevel,
  metricYRangeDf = metricYRangeDf_panel
)

fileNamePlotNow <- paste0(saveDir, scriptID, "MetricChange_CI_9Panel_h2_0.7.pdf")
SaveBoMeanPlotPdf_9Panel(
  df = allResDf_panel,
  pdf_file = fileNamePlotNow,
  h2_value = 0.7,
  minBoForPlot = minBoForPlot,
  ciLevel = ciLevel,
  metricYRangeDf = metricYRangeDf_panel
)
