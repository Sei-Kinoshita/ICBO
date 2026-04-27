###### 0. Settings ######
##### 0.1. Remove objects #####
rm(list = ls())
gc(reset = TRUE);gc(reset = TRUE)

##### 0.2. Load Packages and Functions #####
require(dplyr)
require(tidyr)
require(ggplot2)
require(pracma)
require(multcompView)

##### 0.3. Define Parameters #####
scriptID <- "2.1."

nRep <- 10

kappa <- round(qnorm(p = 0.95), digits = 2) ### parameter for UCB
nRepBo <- 10 ### repetition of BO initial value
indexName <- c("Mu", "PI", "EI", "UCB")
nIndex <- length(indexName)

##### 0.4. Define Directories #####
dataDir <- "data/"
phenoDir <- paste0(dataDir, "phenotype/")
boDir <- "midstream/2.0.BO/Model3_GMA_SMA/"

saveDir <- paste0("midstream/", scriptID, "BO_Evaluation/")
if (!dir.exists(saveDir)) {
  dir.create(saveDir, recursive = TRUE)
}


###### 1. Load Data ######
### scenario parameter set
fileNameParSet <- paste0(phenoDir, "0.3.parmeterSet.rds")
parSetDf <- readRDS(file = fileNameParSet)
nPar <- nrow(parSetDf)

parUniqVec <- (1:nPar) %% 18
parUniqVec[parUniqVec == 0] <- 18
nParUniq <- length(unique(parUniqVec))

### initial best value
fileNameInitBest <- paste0(boDir, "2.0.BO_initPair_bestValue.rds")
initBestLst <- readRDS(file = fileNameInitBest)
### initial pair No
fileNameInitPariNo <- paste0(boDir, "2.0.BO_initPairNoLst.rds")
initPairNoLst <- readRDS(file = fileNameInitPariNo)

### eigen decomposition of GbyG mat (not ZGZ^T)
fileNameGbyG <- "data/genotype/0.3.Kronecker_GRM_AMat_twoCrop_nDim_5214.rds"
gbygMat <- readRDS(file = fileNameGbyG)
fileNameGbyGEig <- "data/genotype/0.3.Eigen_Kronecker_GRM_AMat_twoCrop.rds"
eigGbyG <- readRDS(file = fileNameGbyGEig)
### PCA for GbyG matrix
pc1 <- eigGbyG$vectors[, 1] * sqrt(eigGbyG$values[1])
pc2 <- eigGbyG$vectors[, 2] * sqrt(eigGbyG$values[2])
pcaDf <- data.frame(PC1 = pc1, PC2 = pc2)

### intercrop pair name
fileNamePairDf <- paste0(phenoDir, "0.3.pairNameDf_5214Pair.rds")
pairDf <- readRDS(file = fileNamePairDf)
nPair <- nrow(pairDf)
pairName <- apply(X = pairDf, MARGIN = 1, FUN = paste0, collapse = "_")




###### 2. Evaluate BO Result ######
##### 2.1. Best Value for each iteration #####
fileNameMaxResLst <- paste0(saveDir, scriptID, "BO_allRes_MaxforEachIteartion.rds")
if (!file.exists(fileNameMaxResLst)) {
  boResLst <- lapply(X = 1:nParUniq, FUN = function(parUniqNo) {
    
    parNoNow <- (1:nPar)[parUniqVec == parUniqNo]
    
    oneLst <- lapply(X = parNoNow, FUN = function(no) {
      
      ### directory of each index's result
      parDir <- paste0(boDir, "parNo_", no, "/")
      muDir <- paste0(parDir, "Mu/BO_result/")
      piDir <- paste0(parDir, "PI/BO_result/")
      eiDir <- paste0(parDir, "EI/BO_result/")
      # ucbDir <- paste0(parDir, "UCB_", kappa, "/BO_result/")
      ucbDir <- paste0(parDir, "UCB/BO_result/")
      
      ### load BO result
      oneParLst <- lapply(X = 1:nRepBo, FUN = function(boRepNo) {
        ### each index
        fileNameBoRes <- paste0("2.0.g_best_allIteration_BORepNo_", boRepNo, ".rds")
        boResMu <- readRDS(file = paste0(muDir, fileNameBoRes))
        boResPi <- readRDS(file = paste0(piDir, fileNameBoRes))
        boResEi <- readRDS(file = paste0(eiDir, fileNameBoRes))
        boResUcb <- readRDS(file = paste0(ucbDir, fileNameBoRes))
        len <- c(length(boResMu), length(boResPi), length(boResEi), length(boResUcb))
        ### to same length (nIter)
        if (length(unique(len)) != 1) {
          maxLen <- max(len)
          boResMu <- c(boResMu, rep(tail(boResMu, 1), maxLen - length(boResMu)))
          boResPi <- c(boResPi, rep(tail(boResPi, 1), maxLen - length(boResPi)))
          boResEi <- c(boResEi, rep(tail(boResEi, 1), maxLen - length(boResEi)))
          boResUcb <- c(boResUcb, rep(tail(boResUcb, 1), maxLen - length(boResUcb)))
        }
        ### initial best value
        initBest <- initBestLst[[boRepNo]][no]
        ### result for all index
        bindRes <- rbind(c(initBest, boResMu), c(initBest, boResPi), 
                         c(initBest, boResEi), c(initBest, boResUcb))
        nIter <- ncol(bindRes) - 1
        resDf <- data.frame(bindRes, 
                            Index = indexName, 
                            BoRepNo = rep(boRepNo, nIndex), 
                            ParNo = no)
        names(resDf) <- c(paste0("Iter_", 0:nIter), "Index", "BoRepNo", "ParNo")
        
        return(resDf)
      })
      
      ### match iteration length
      maxLen <- max(sapply(oneParLst, ncol))
      for (i in 1:length(oneParLst)) {
        resDfNow <- oneParLst[[i]]
        k <- maxLen - ncol(resDfNow)
        nIter <- ncol(resDfNow) - 3
        if (k > 0) {
          fillCol <- resDfNow[, nIter, drop = FALSE]
          add <- as.data.frame(
            matrix(rep(as.matrix(fillCol), k), nrow = nrow(resDfNow), ncol = k)
          )
          colnames(add) <- paste0("Iter_", 
                                  (nIter):(nIter - 1 + k))
          oneParLst[[i]] <- cbind(resDfNow[, 1:(nIter)], 
                                  add, 
                                  resDfNow[, (nIter + 1):(nIter + 3)])
        }
      }
      
      return(do.call(what = rbind, args = oneParLst))
    })
    
    ### match iteration length
    maxLen <- max(sapply(oneLst, ncol))
    for (i in 1:length(oneLst)) {
      resDfNow <- oneLst[[i]]
      k <- maxLen - ncol(resDfNow)
      nIter <- ncol(resDfNow) - 3
      if (k > 0) {
        fillCol <- resDfNow[, nIter, drop = FALSE]
        add <- as.data.frame(
          matrix(rep(as.matrix(fillCol), k), nrow = nrow(resDfNow), ncol = k)
        )
        colnames(add) <- paste0("Iter_", 
                                (nIter):(nIter - 1 + k))
        oneLst[[i]] <- cbind(resDfNow[, 1:(nIter)], 
                             add, 
                             resDfNow[, (nIter + 1):(nIter + 3)])
      }
    }
    
    return(do.call(what = rbind, args = oneLst))
  })
  ### save
  saveRDS(object = boResLst, file = fileNameMaxResLst)
} else {
  boResLst <- readRDS(file = fileNameMaxResLst)
}


##### 2.2. AUC for each scenario and index #####
boAucLst <- lapply(X = boResLst, FUN = function(oneParDf) {
  ### number of iterations (include 0)
  nIterNow <- ncol(oneParDf) - 3
  maxAllIterMat <- as.matrix(oneParDf[, 1:nIterNow])
  ### AUC for each scenairo
  aucNow <- apply(X = maxAllIterMat, MARGIN = 1, FUN = function(maxVal) {
    endCycle <- which.max(maxVal)
    trapz(x = 0:(nIterNow - 1), y = maxVal[1:nIterNow])
  })
  oneParDf$AUC <- aucNow
  
  return(oneParDf[, (ncol(oneParDf) - 3):ncol(oneParDf)])
})


##### 2.3. Number of iteration to the best value #####
fileNameNiterResLst <- paste0(saveDir, scriptID, "BO_allRes_nIteartion.rds")
if (!file.exists(fileNameNiterResLst)) {
  boIterDf <- do.call(what = rbind, args = lapply(X = 1:nPar, FUN = function(parNo) {
    
    ### parmeter set
    genCor <- parSetDf$genCor[parNo]
    h2 <- parSetDf$h2[parNo]
    gamma <- parSetDf$gamma[parNo]
    repNo <- parSetDf$rep[parNo]
    
    ### directory of each index's result
    parDir <- paste0(boDir, "parNo_", parNo, "/")
    muDir <- paste0(parDir, "Mu/BO_result/")
    piDir <- paste0(parDir, "PI/BO_result/")
    eiDir <- paste0(parDir, "EI/BO_result/")
    # ucbDir <- paste0(parDir, "UCB_", kappa, "/BO_result/")
    ucbDir <- paste0(parDir, "UCB/BO_result/")
    
    ### load BO result
    oneParDf <- do.call(what = rbind, 
                        args = oneParLst <- lapply(X = 1:nRepBo, FUN = function(boRepNo) {
                          ### each index
                          fileNameBoRes <- paste0("2.0.g_best_allIteration_BORepNo_", boRepNo, ".rds")
                          boResMu <- readRDS(file = paste0(muDir, fileNameBoRes))
                          boResPi <- readRDS(file = paste0(piDir, fileNameBoRes))
                          boResEi <- readRDS(file = paste0(eiDir, fileNameBoRes))
                          boResUcb <- readRDS(file = paste0(ucbDir, fileNameBoRes))
                          len <- c(length(boResMu), length(boResPi), 
                                   length(boResEi), length(boResUcb))
                          ### output
                          resDf <- data.frame(nIter = len, Index = indexName, 
                                              GenCor = genCor, h2 = h2, Gamma = gamma, 
                                              BoRepNo = boRepNo, ParNo = parNo)
                          
                          return(resDf)
                        }))
    
    return(oneParDf)
  }))
  ### save
  saveRDS(object = boIterDf, file = fileNameNiterResLst)
} else {
  boIterDf <- readRDS(file = fileNameNiterResLst)
}




###### 3. Plot ######
##### 3.1. Selected Pair for each parNo #####
fileNameGrmPlot <- paste0(saveDir, scriptID, "BORes_GRM_selPairPlot.pdf")
if (!file.exists(fileNameGrmPlot)) {
  pdf(file = fileNameGrmPlot, width = 15, height = 12)
  grmPlotLst <- pbmcapply::pbmclapply(X = 1:nPar, FUN = function(parNo) {
    
    # parNo <- 1
    
    ### parmeter set
    genCor <- parSetDf$genCor[parNo]
    h2 <- parSetDf$h2[parNo]
    gamma <- parSetDf$gamma[parNo]
    repNo <- parSetDf$rep[parNo]
    
    ### directory of each index's result
    parDir <- paste0(boDir, "parNo_", parNo, "/")
    
    oneParLst <- lapply(X = 1:nRepBo, FUN = function(boRepNo) {
      
      # boRepNo <- 1
      
      fileNameBoRes <- paste0("2.0.selectedPairNo_allIteration_BORepNo_", boRepNo, ".rds")
      
      plotDfLst <- lapply(X = indexName, FUN = function(indexNameNow){
        
        indexDirNow <- paste0(parDir, indexNameNow, "/BO_result/")
        boResNow <- readRDS(file = paste0(indexDirNow, fileNameBoRes))
        ### unselected pair
        bgDf <- data.frame(PC1 = pcaDf$PC1, PC2 = pcaDf$PC2,
                           Index = indexNameNow, Type = "Background",
                           CycleNo = NA, Size = 0.8)
        ### initial selected pair
        initDf <- data.frame(PC1 = pcaDf$PC1[boResNow[[1]]],
                             PC2 = pcaDf$PC2[boResNow[[1]]],
                             Index = indexNameNow, Type = "Initial",
                             CycleNo = 1, Size = 2)
        ### divide selected pair by cycle number
        if (length(boResNow) >= 2) {
          cycleDfLst <- lapply(X = 2:length(boResNow), FUN = function(i){
            data.frame(PC1 = pcaDf$PC1[boResNow[[i]]],
                       PC2 = pcaDf$PC2[boResNow[[i]]],
                       Index = indexNameNow, Type = "Cycle",
                       CycleNo = i, Size = 2)
          })
          outDf <- do.call(rbind, c(list(bgDf, initDf), cycleDfLst))
        } else {
          outDf <- rbind(bgDf, initDf)
        }
        
        outDf
      })
      plotDf <- do.call(rbind, plotDfLst)
      plotDf$Index <- factor(plotDf$Index, levels = indexName)
      
      ### plot
      mainTitle <- paste0("genCor: ", genCor, "  h2: ", h2,
                          "  gamma: ", gamma, "  phenoRep: ", repNo,
                          "  boRep: ", boRepNo)
      ggplot() + 
        geom_point(data = subset(plotDf, Type == "Background"),
                   aes(x = PC1, y = PC2),
                   color = "grey85", size = 0.8) +
        geom_point(data = subset(plotDf, Type == "Initial"),
                   aes(x = PC1, y = PC2), color = "black", size = 2) +
        geom_point(data = subset(plotDf, Type == "Cycle"),
                   aes(x = PC1, y = PC2, color = CycleNo), size = 2, alpha = 0.6) +
        scale_color_viridis_c() +
        facet_wrap(~Index, ncol = 2) +
        labs(x = "PC1", y = "PC2", color = "Cycle", title = mainTitle) +
        theme_bw() +
        theme(axis.text = element_text(size = 14),
              axis.title = element_text(size = 16),
              plot.title = element_text(size = 16))
    })
  }, mc.cores = 30)
  dev.off()
}


##### 3.2. Change of Best Value #####
fileNameBoPlot <- paste0(saveDir, scriptID, "BORes_CompareIndex_plot.pdf")
if (!file.exists(fileNameBoPlot)) {
  pdf(file = fileNameBoPlot, width = 15, height = 12)
  lapply(X = 1:length(boResLst), FUN = function(resNo) {
    
    df <- boResLst[[resNo]]
    df$Index <- factor(df$Index, levels = indexName)
    
    timeCols <- grep("^Iter_\\d+$", names(df), value = TRUE)
    timeCols <- timeCols[order(as.integer(sub("^Iter_", "", timeCols)))]
    timeVec <- sort(as.integer(sub("^Iter_", "", timeCols)))
    
    long <- df %>%
      pivot_longer(cols = all_of(timeCols), names_to = "Iter", values_to = "y")
    ### mean over all ParRep and BoRep
    meanDf <- long %>%
      group_by(Index, ParNo, Iter) %>%
      summarise(m_par = mean(y, na.rm = TRUE), .groups = "drop") %>%
      group_by(Index, Iter) %>%
      summarise(mean = mean(m_par, na.rm = TRUE), .groups = "drop")
    ### for error bar over all ParRep and BoRep
    errDf <- long %>%
      group_by(Index, Iter) %>%
      summarise(
        se = sd(y, na.rm = TRUE) / sqrt(sum(!is.na(y))),
        .groups = "drop"
      )
    ### bind for plot
    sumDf <- meanDf %>%
      left_join(errDf, by = c("Index", "Iter")) %>%
      # mutate(Iter = factor(Iter, levels = timeCols))
      mutate(Iter = as.integer(sub("^Iter_", "", Iter)),
             Iter = factor(Iter, levels = timeVec))
    ### plot
    plotTitle <- paste0("genCor: ", parSetDf[resNo, 1], 
                        " h2: ", parSetDf[resNo, 2], 
                        " gamma: ", parSetDf[resNo, 3])
    g <- ggplot(sumDf, aes(Iter, mean, color = Index, group = Index)) +
      geom_line(linewidth = 1.7) +
      geom_point(size = 2) +
      geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.5) +
      scale_x_discrete(breaks = timeVec, labels = timeVec) + 
      labs(title = plotTitle, y = "Best Value") + 
      theme_bw() +
      theme(axis.text.y = element_text(size = 20), 
            axis.text.x = element_text(angle = 60, hjust = 1, size = 18))
    print(g)
    
    return(NULL)
  })
  dev.off()
}


##### 3.3. AUC for each scenario #####
fileNameAucPlot <- paste0(saveDir, scriptID, "BORes_AUC_CompareIndex_plot.pdf")
if (!file.exists(fileNameAucPlot)) {
  pdf(file = fileNameAucPlot, width = 15, height = 12)
  fileNameFrdRes <- paste0(saveDir, scriptID, "BORes_AUC_Friedman_Wilcoxon_test.txt")
  txtPath <- file(fileNameFrdRes, open = "wt")
  lapply(X = 1:length(boAucLst), FUN = function(resNo) {
    
    oneParDf <- boAucLst[[resNo]]
    oneParDf$Index <- factor(oneParDf$Index, levels = indexName)
    
    # ### mean AUC for each parameter
    # meanAuc <- tapply(X = oneParDf$AUC, INDEX = oneParDf$Index, FUN = mean)
    
    ### statistical test
    ### Friedman test for AUC for each index
    modiData <- data.frame(Index = oneParDf$Index, AUC = oneParDf$AUC, 
                           Rep = rep(1:(nrow(oneParDf) / 4), each = 4))
    frdRes <- friedman.test(AUC ~ Index | Rep, data = modiData)
    
    ### Wilcoxon test for each index pair
    wideDat <- matrix(oneParDf$AUC, ncol = nIndex, byrow = TRUE)
    colnames(wideDat) <- indexName
    wideDat <- as.data.frame(wideDat)
    pVals <- c(
      wilcox.test(wideDat$Mu, wideDat$PI, paired = TRUE, exact = FALSE)$p.value, 
      wilcox.test(wideDat$Mu, wideDat$EI, paired = TRUE, exact = FALSE)$p.value, 
      wilcox.test(wideDat$Mu, wideDat$UCB, paired = TRUE, exact = FALSE)$p.value, 
      wilcox.test(wideDat$PI, wideDat$EI, paired = TRUE, exact = FALSE)$p.value, 
      wilcox.test(wideDat$PI, wideDat$UCB, paired = TRUE, exact = FALSE)$p.value, 
      wilcox.test(wideDat$EI, wideDat$UCB, paired = TRUE, exact = FALSE)$p.value
    )
    pValAdj <- p.adjust(pVals, method = "bonferroni")
    pValAdjRes <- data.frame(Pair = c("Mu_PI", "Mu_EI", "Mu_UCB", "PI_EI", "PI_UCB", "EI_UCB"), 
                             pAdj = pValAdj)
    
    ### add letters
    names(pValAdj) <- c("Mu-PI", "Mu-EI", "Mu-UCB", "PI-EI", "PI-UCB", "EI-UCB")
    letVec <- multcompView::multcompLetters(pValAdj, threshold = 0.05)$Letters
    ypos <- tapply(oneParDf$AUC, oneParDf$Index, max, na.rm = TRUE)
    gap <- 0.05 * diff(range(oneParDf$AUC, na.rm = TRUE))
    labDf <- data.frame(Index = factor(names(letVec), levels = indexName),
                        Lab = letVec,
                        Y = ypos[names(letVec)] + gap)
    
    ### plot
    plotTitle <- paste0("genCor: ", parSetDf[resNo, 1], 
                        " h2: ", parSetDf[resNo, 2], 
                        " gamma: ", parSetDf[resNo, 3])
    
    g <- ggplot(data = oneParDf, aes(x = Index, y = AUC, fill = Index)) + 
      geom_boxplot() + 
      geom_text(data = labDf, aes(x = Index, y = Y, label = Lab),
                inherit.aes = FALSE, size = 6) + 
      expand_limits(y = max(labDf$Y) + gap) + 
      labs(title = plotTitle) + 
      theme_bw() + 
      theme(text = element_text(size = 17))
    print(g)
    
    out <- c(
      plotTitle, "", 
      "=== Friedman test ===", 
      capture.output(print(frdRes)), 
      "", 
      "=== Pairwise comparison by Wilcoxon", 
      capture.output(print(pValAdjRes)), 
      ""
    )
    
    writeLines(out, txtPath)
  })
  close(txtPath)
  dev.off()
}


##### 3.4. Mean Number of BO cycles #####
boIterDf$GenCor <- factor(boIterDf$GenCor)
boIterDf$h2 <- factor(boIterDf$h2)
boIterDf$Gamma <- factor(boIterDf$Gamma)
boIterDf$Index <- factor(boIterDf$Index, levels = c("Mu", "PI", "EI", "UCB"))
# ######################
# boIterDf$nIter <- 47.054 - boIterDf$nIter
# ######################
# boIterDf$Speed <- 1 / (boIterDf$nIter)
# boIterDf$Speed <- 47.054 / boIterDf$nIter

### mean iteration each parNo
fileNameBoIterPlot <- paste0(saveDir, scriptID, "BOnIter_CompareIndex_plot.pdf")
pdf(file = fileNameBoIterPlot, width = 15, height = 10)
meanIterDf <- sapply(X = 1:nParUniq, FUN = function(parUniqNo) {
  ### combine same parameter set
  parUniqNoNow <- which(parUniqVec == parUniqNo)
  ### parmeter set
  genCor <- parSetDf$genCor[parUniqNoNow[1]]
  h2 <- parSetDf$h2[parUniqNoNow[1]]
  gamma <- parSetDf$gamma[parUniqNoNow[1]]
  
  ### plot
  plotDf <- boIterDf[boIterDf$ParNo %in% parUniqNoNow, ]
  
  plotTitle <- paste0("genCor: ", plotDf$GenCor[1], " h2: ", plotDf$h2[1],
                      " gamma: ", plotDf$Gamma[1])
  meanIter <- tapply(X = plotDf$nIter, INDEX = plotDf$Index, FUN = mean)
  
  ### statistical test
  repid <- sort(unique(plotDf$ParNo))
  blkdf <- do.call(what = rbind, args = lapply(X = 1:nRepBo, FUN = function(bono) {
    tmp <- plotDf[plotDf$BoRepNo == bono, c("ParNo", "Index", "nIter")]
    tmp$Rep <- match(tmp$ParNo, repid) + (bono - 1) * length(repid)
    tmp
  }))
  frdres <- friedman.test(nIter ~ Index | Rep, data = blkdf)
  
  widedat <- reshape(blkdf[, c("Rep", "Index", "nIter")],
                     idvar = "Rep", timevar = "Index", direction = "wide")
  names(widedat) <- sub("^nIter\\.", "", names(widedat))
  pvals <- c(
    wilcox.test(widedat$Mu, widedat$PI, paired = TRUE, exact = FALSE)$p.value,
    wilcox.test(widedat$Mu, widedat$EI, paired = TRUE, exact = FALSE)$p.value,
    wilcox.test(widedat$Mu, widedat$UCB, paired = TRUE, exact = FALSE)$p.value,
    wilcox.test(widedat$PI, widedat$EI, paired = TRUE, exact = FALSE)$p.value,
    wilcox.test(widedat$PI, widedat$UCB, paired = TRUE, exact = FALSE)$p.value,
    wilcox.test(widedat$EI, widedat$UCB, paired = TRUE, exact = FALSE)$p.value
  )
  padj <- p.adjust(pvals, method = "bonferroni")
  names(padj) <- c("Mu-PI", "Mu-EI", "Mu-UCB", "PI-EI", "PI-UCB", "EI-UCB")
  letvec <- multcompView::multcompLetters(padj, threshold = 0.05)$Letters
  ypos <- tapply(plotDf$nIter, plotDf$Index, max, na.rm = TRUE)
  gap <- 0.05 * diff(range(plotDf$nIter, na.rm = TRUE))
  if (gap == 0) {
    gap <- 0.5
  }
  labdf <- data.frame(Index = factor(names(letvec), levels = indexName),
                      Lab = letvec,
                      Y = ypos[names(letvec)] + gap)
  
  g <- ggplot(data = plotDf, aes(x = Index, y = nIter, fill = Index)) +
    geom_boxplot() +
    geom_text(data = labdf, aes(x = Index, y = Y, label = Lab),
              inherit.aes = FALSE, size = 6) +
    expand_limits(y = max(labdf$Y) + gap) +
    labs(title = plotTitle, y = "# of iterations")
  print(g)
  
  return(meanIter)
})
colnames(meanIterDf) <- paste0("ParNo_", 1:18)

### mean iteration each scenario
ggplot(data = boIterDf, aes(x = Index, y = nIter, fill = h2)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~GenCor)

ggplot(data = boIterDf, aes(x = Index, y = nIter, fill = GenCor)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~h2)

ggplot(data = boIterDf, aes(x = Index, y = nIter, fill = h2)) +
  geom_bar(stat = "identity", position = "dodge")

ggplot(data = boIterDf, aes(x = Index, y = nIter, fill = Gamma)) +
  geom_bar(stat = "identity", position = "dodge")

ggplot(data = boIterDf, aes(x = Index, y = nIter, fill = GenCor)) +
  geom_bar(stat = "identity", position = "dodge")

### 9 cells barplot
ggplot(data = boIterDf, aes(x = h2, y = nIter, fill = Index)) +
  stat_summary(fun = mean,
               geom = "bar",
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = mean_se,
               geom = "errorbar",
               position = position_dodge(width = 0.9),
               width = 0.2) +
  facet_grid(GenCor ~ Gamma) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),
        strip.text = element_text(size = 17))

### 9 cells barplot + letters
plotDf <- boIterDf
plotDf$GenCor <- factor(plotDf$GenCor, levels = sort(unique(plotDf$GenCor)))
plotDf$Gamma  <- factor(plotDf$Gamma,  levels = sort(unique(plotDf$Gamma)))
plotDf$h2     <- factor(plotDf$h2,     levels = sort(unique(plotDf$h2)))
plotDf$Index  <- factor(plotDf$Index,  levels = indexName)

### summary for barplot
sumDf <- aggregate(nIter ~ GenCor + Gamma + h2 + Index,
                   data = plotDf,
                   FUN = mean)
seDf <- aggregate(nIter ~ GenCor + Gamma + h2 + Index,
                  data = plotDf,
                  FUN = function(x) sd(x) / sqrt(length(x)))
names(sumDf)[names(sumDf) == "nIter"] <- "Mean"
names(seDf)[names(seDf) == "nIter"] <- "SE"
sumDf <- merge(sumDf, seDf,
               by = c("GenCor", "Gamma", "h2", "Index"),
               sort = FALSE)

### --------------------------------
### statistical letters for each cell
### cell = one GenCor x Gamma x h2
### --------------------------------
labDf <- do.call(rbind, lapply(X = split(plotDf,
                                         list(plotDf$GenCor, plotDf$Gamma, plotDf$h2),
                                         drop = TRUE),
                               FUN = function(datNow) {
                                 
                                 datNow$Index <- factor(datNow$Index, levels = indexName)
                                 
                                 ### block for paired test:
                                 ### one block = same phenotype replicate (ParNo) and BO replicate (BoRepNo)
                                 repid <- sort(unique(datNow$ParNo))
                                 blkdf <- do.call(rbind, lapply(X = 1:nRepBo, FUN = function(bono) {
                                   tmp <- datNow[datNow$BoRepNo == bono, c("ParNo", "Index", "nIter")]
                                   tmp$Rep <- match(tmp$ParNo, repid) + (bono - 1) * length(repid)
                                   tmp
                                 }))
                                 
                                 ### Friedman test
                                 frdres <- friedman.test(nIter ~ Index | Rep, data = blkdf)
                                 
                                 ### paired Wilcoxon
                                 widedat <- reshape(blkdf[, c("Rep", "Index", "nIter")],
                                                    idvar = "Rep", timevar = "Index", direction = "wide")
                                 names(widedat) <- sub("^nIter\\.", "", names(widedat))
                                 
                                 pvals <- c(
                                   wilcox.test(widedat$Mu,  widedat$PI,  paired = TRUE, exact = FALSE)$p.value,
                                   wilcox.test(widedat$Mu,  widedat$EI,  paired = TRUE, exact = FALSE)$p.value,
                                   wilcox.test(widedat$Mu,  widedat$UCB, paired = TRUE, exact = FALSE)$p.value,
                                   wilcox.test(widedat$PI,  widedat$EI,  paired = TRUE, exact = FALSE)$p.value,
                                   wilcox.test(widedat$PI,  widedat$UCB, paired = TRUE, exact = FALSE)$p.value,
                                   wilcox.test(widedat$EI,  widedat$UCB, paired = TRUE, exact = FALSE)$p.value
                                 )
                                 padj <- p.adjust(pvals, method = "bonferroni")
                                 names(padj) <- c("Mu-PI", "Mu-EI", "Mu-UCB", "PI-EI", "PI-UCB", "EI-UCB")
                                 
                                 letvec <- multcompView::multcompLetters(padj, threshold = 0.05)$Letters
                                 
                                 ### y position
                                 ymaxDf <- sumDf[sumDf$GenCor == datNow$GenCor[1] &
                                                   sumDf$Gamma  == datNow$Gamma[1] &
                                                   sumDf$h2     == datNow$h2[1], ]
                                 ypos <- ymaxDf$Mean + ymaxDf$SE
                                 names(ypos) <- as.character(ymaxDf$Index)
                                 
                                 gap <- 0.05 * diff(range(plotDf$nIter, na.rm = TRUE))
                                 if (gap == 0) {
                                   gap <- 0.3
                                 }
                                 
                                 outDf <- data.frame(GenCor = datNow$GenCor[1],
                                                     Gamma  = datNow$Gamma[1],
                                                     h2     = datNow$h2[1],
                                                     Index  = factor(names(letvec), levels = indexName),
                                                     Lab    = letvec,
                                                     Y      = ypos[names(letvec)] + gap)
                                 
                                 return(outDf)
                               }))

### x position for letters
h2Lev <- levels(plotDf$h2)
idxLev <- levels(plotDf$Index)
nIdx <- length(idxLev)
dodgeWidth <- 0.9

labDf$x <- NA
for (i in 1:nrow(labDf)) {
  h2Pos <- match(as.character(labDf$h2[i]), h2Lev)
  idxPos <- match(as.character(labDf$Index[i]), idxLev)
  labDf$x[i] <- h2Pos + ((idxPos - (nIdx + 1) / 2) / nIdx) * dodgeWidth
}

### plot
ggplot(data = sumDf, aes(x = h2, y = Mean, fill = Index)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = dodgeWidth)) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                position = position_dodge(width = dodgeWidth),
                width = 0.2) +
  geom_text(data = labDf,
            aes(x = x, y = Y, label = Lab),
            inherit.aes = FALSE,
            size = 5) +
  facet_grid(GenCor ~ Gamma) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),
        axis.title = element_text(size = 17),
        strip.text = element_text(size = 17),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 16)) +
  labs(x = expression(h^2),
       y = "# of iterations")

dev.off()


##### 3.5. Min distance to observed set #####
diagVec <- diag(gbygMat)
distMat <- sqrt(outer(diagVec, diagVec, "+") - 2 * gbygMat)

fileNameDistPlot <- paste0(saveDir, scriptID, "BOPlot_Dist_gTrueCompare.pdf")
if (!file.exists(fileNameDistPlot)) {
  pdf(file = fileNameDistPlot, width = 15, height = 12)
  distPlotLst <- lapply(X = 1:nParUniq, FUN = function(parUniqNo) {
    
    parNoNow <- (1:nPar)[parUniqVec == parUniqNo]
    genCor <- parSetDf$genCor[parNoNow[1]]
    h2 <- parSetDf$h2[parNoNow[1]]
    gamma <- parSetDf$gamma[parNoNow[1]]
    
    ### load phenotype to calc true genotypic values
    fileNamePheno <- paste0(phenoDir, "0.3.allSampledRes_phenotype_genCor_", genCor, 
                            "_h2_", h2, "_gamma_", gamma, "_nRep_", nRep, ".rds")
    ### gTrue
    gTrueLst <- lapply(X = 1:nRep, FUN = function(repNo) {
      phenoLst <- readRDS(file = fileNamePheno)[[repNo]]
      ### true best pair
      gTrueA <- scale(phenoLst$dA + phenoLst$sBtoA + phenoLst$vA)
      gTrueB <- scale(phenoLst$dB + phenoLst$sAtoB + phenoLst$vB)
      gTrue <- gTrueA + gTrueB
      return(gTrue)
    })
    
    plotDfLst <- lapply(X = parNoNow, FUN = function(no) {
      
      parDir <- paste0(boDir, "parNo_", no, "/")
      repNoNow <- parSetDf$rep[no]
      
      ### gTrueNow
      trueGVec <- gTrueLst[[repNoNow]]
      
      do.call(what = rbind, args = lapply(X = indexName, FUN = function(indexNameNow){
        
        boDfLst <- lapply(X = 1:nRepBo, FUN = function(boRepNoNow){
          
          fileNameBoRes <- paste0("2.0.selectedPairNo_allIteration_BORepNo_", boRepNoNow, ".rds")
          indexDirNow <- paste0(parDir, indexNameNow, "/BO_result/")
          boResNow <- readRDS(file = paste0(indexDirNow, fileNameBoRes))
          
          cycleDfLst <- lapply(X = seq_along(boResNow)[-1], FUN = function(i){
            
            selNow <- boResNow[[i]]
            obsNow <- unique(unlist(boResNow[seq_len(i - 1)]))
            
            minDistNow <- apply(distMat[selNow, obsNow, drop = FALSE], 1, min)
            trueGMax <- max(trueGVec)
            trueGMin <- min(trueGVec)
            trueGNow <- (trueGMax - trueGVec[selNow]) / (trueGMax - trueGMin)
            
            data.frame(Index = indexNameNow, RepNo = repNoNow,
                       BoRepNo = boRepNoNow, CycleNo = i, PairNo = selNow,
                       MinDist = minDistNow, RelGap = trueGNow)
          })
          
          do.call(rbind, cycleDfLst)
        })
        
        do.call(rbind, boDfLst)
      }))
    })
    
    plotDf <- do.call(what = rbind, args = plotDfLst)
    plotDf$Index <- factor(plotDf$Index, levels = indexName)
    
    mainTitle <- paste0("genCor: ", genCor, "  h2: ", h2, "  gamma: ", gamma)
    
    ### min distance v.s. trueGV gap (all points)
    # p1 <- ggplot(plotDf, aes(x = MinDist, y = RelGap, color = CycleNo)) +
    #   geom_point(alpha = 0.15, size = 0.7) +
    #   scale_color_viridis_c() +
    #   facet_wrap(~Index, ncol = 2) +
    #   labs(x = "Min distance to observed set",
    #        y = "Relative gap from true maximum",
    #        color = "Cycle",
    #        title = mainTitle) +
    #   theme_bw() +
    #   theme(axis.text = element_text(size = 14),
    #         axis.title = element_text(size = 16),
    #         strip.text = element_text(size = 14),
    #         plot.title = element_text(size = 16))
    
    ### min distance v.s. trueGV gap (only median)
    repSumDf <- aggregate(cbind(MinDist, RelGap) ~ Index + RepNo + BoRepNo + CycleNo,
                          data = plotDf, FUN = median)
    allSumDf <- aggregate(cbind(MinDist, RelGap) ~ Index + CycleNo,
                          data = repSumDf, FUN = median)
    p1 <- ggplot() +
      geom_point(data = repSumDf,
                 aes(x = MinDist, y = RelGap, color = CycleNo),
                 alpha = 0.25, size = 1.2) +
      geom_path(data = allSumDf,
                aes(x = MinDist, y = RelGap, group = Index),
                color = "black", linewidth = 1) +
      geom_point(data = allSumDf,
                 aes(x = MinDist, y = RelGap, fill = CycleNo),
                 shape = 21, color = "black", size = 3) +
      scale_color_viridis_c() +
      scale_fill_viridis_c() +
      facet_wrap(~Index, ncol = 2) +
      labs(x = "Min distance to observed set",
           y = "Relative gap from true maximum",
           color = "Cycle", fill = "Cycle", title = mainTitle) +
      theme_bw() +
      theme(axis.text = element_text(size = 14),
            axis.title = element_text(size = 16),
            strip.text = element_text(size = 14))
    
    ### min distance in each cycle
    p2 <- ggplot(plotDf, aes(x = factor(CycleNo), y = MinDist, color = CycleNo)) +
      geom_jitter(width = 0.25, height = 0, alpha = 0.15, size = 0.7) +
      scale_color_viridis_c() +
      facet_wrap(~Index, ncol = 2) +
      labs(x = "Cycle",
           y = "Min distance to observed set",
           color = "Cycle", title = mainTitle) +
      theme_bw() +
      theme(axis.text = element_text(size = 14),
            axis.title = element_text(size = 16),
            strip.text = element_text(size = 14))
    
    print(p1)
    print(p2)
  })
  dev.off()
}