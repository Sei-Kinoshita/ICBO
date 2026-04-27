###### 0. Settings ######
##### 0.1. Remove objects #####
rm(list = ls())
gc(reset = TRUE);gc(reset = TRUE)

##### 0.2. Load Packages and Functions #####
require(BGLR)
source("scripts/function_scripts/2.0.function.BO.R")

##### 0.3. Define Parameters #####
scriptID <- "2.0."
cropName <- c("Sb", "Gm")
### scenario parameter
genCorVec <- c(-0.7, -0.1, 0.7)
h2Vec <- c(0.3, 0.7)
gammaVec <- c(0, 0.2, 0.4)
nRep <- 10 ### repetition of phenotype
nRepBo <- 10 ### repetition of BO initial value
### number of observed individuals at first
nInitSb <- 8
nInitGm <- 7
nInitComb <- nInitSb * nInitGm
### number of new observed individuals at each iteration
nPickSb <- nInitSb
nPickGm <- nInitGm
nPickComb <- nInitComb
### parameter for UCB
kappa <- round(qnorm(p = 0.95), digits = 2)
### parameter for BGLR
nIter <- 10000
burnIn <- 2000
thin <- 20
# nIter <- 10
# burnIn <- 2

nCore <- 40

##### 0.4. Define Directories #####
dataDir <- "data/"
genoDir <- paste0(dataDir, "genotype/")
phenoDir <- paste0(dataDir, "phenotype/")

saveDirRandom <- paste0("midstream/", scriptID, "BO/RandomSearch/")
if (!dir.exists(saveDirRandom)) {
  dir.create(saveDirRandom, recursive = TRUE)
}
saveDirM3 <- paste0("midstream/", scriptID, "BO/Model3_GMA_SMA/")
if (!dir.exists(saveDirM3)) {
  dir.create(saveDirM3, recursive = TRUE)
}
saveDirM5 <- paste0("midstream/", scriptID, "BO/Model5_Multi_GMA_SMA/")
if (!dir.exists(saveDirM5)) {
  dir.create(saveDirM5, recursive = TRUE)
}




###### 1. Load Data ######
### K mat (ZGZ^T)
fileNameKSb <- "midstream/1.0.GP_Fitting/1.0.Sb_KMat.rds"
fileNameKGm <- "midstream/1.0.GP_Fitting/1.0.Gm_KMat.rds"
fileNameKGbyG <- "midstream/1.0.GP_Fitting/1.0.GbyG_KMat.rds"
K1 <- readRDS(file = fileNameKSb)
K2 <- readRDS(file = fileNameKGm)
K3 <- readRDS(file = fileNameKGbyG)

### intercrop pair name
fileNamePairDf <- paste0(phenoDir, "0.3.pairNameDf_5214Pair.rds")
pairDf <- readRDS(file = fileNamePairDf)
nPair <- nrow(pairDf)
pairName <- apply(X = pairDf, MARGIN = 1, FUN = paste0, collapse = "_")

sbName <- unique(pairDf$Sb)
gmName <- unique(pairDf$Gm)
nSb <- length(sbName)
nGm <- length(gmName)

### scenario parameter set
fileNameParSet <- paste0(phenoDir, "0.3.parmeterSet.rds")
parSetDf <- readRDS(file = fileNameParSet)
nPar <- nrow(parSetDf)




###### 2. Prepare for BO ######
### number of iterations for BO
nIterBo <- floor((nPair - nInitComb) / nPickComb)

### set seed for initial pairs of BO
fileNameSeedBo <- paste0(saveDirM3, scriptID, "seed_BO_initPair.rds")
if (!file.exists(fileNameSeedBo)) {
  seed <- sample(x = 1e06, size = nRepBo, replace = FALSE)
  saveRDS(object = seed, file = fileNameSeedBo)
} else {
  seed <- readRDS(file = fileNameSeedBo)
}
fileNameInitPairNo <- paste0(saveDirM3, scriptID, "BO_initPairNoLst.rds")
if (!file.exists(fileNameInitPairNo)) {
  initPairNoLst <- lapply(X = 1:nRepBo, FUN = function(boRepNo) {
    set.seed(seed[boRepNo])
    sample(x = 1:nPair, size = nInitComb, replace = FALSE)
  })
  saveRDS(object = initPairNoLst, file = fileNameInitPairNo)
} else {
  initPairNoLst <- readRDS(file = fileNameInitPairNo)
}

### save best value of initial BO pairs
fileNameBoInitPair <- paste0(saveDirM3, scriptID, "BO_initPair_bestValue.rds")
if (!file.exists(fileNameBoInitPair)) {
  ### for best pair in initial pair
  initPairUpdated <- FALSE
  
  bestInitPairLst <- lapply(X = 1:nRepBo, FUN = function(boRepNo) {
    unlist(lapply(X = 1:nPar, FUN = function(parNo) {
      ### scenario parameter
      genCor <- parSetDf$genCor[parNo]
      h2 <- parSetDf$h2[parNo]
      gamma <- parSetDf$gamma[parNo]
      repNo <- parSetDf$rep[parNo]
      ### load phenotype
      fileNamePheno <- paste0(phenoDir, "0.3.allSampledRes_phenotype_genCor_", genCor, 
                              "_h2_", h2, "_gamma_", gamma, "_nRep_", nRep, ".rds")
      phenoLst <- readRDS(file = fileNamePheno)[[repNo]]
      ### genotypic value and scaling
      gA <- phenoLst$dA + phenoLst$sBtoA + phenoLst$vA
      gB <- phenoLst$dB + phenoLst$sAtoB + phenoLst$vB
      gAScaled <- scale(gA)
      gBScaled <- scale(gB)
      ### initial pairNo
      initPairNo <- initPairNoLst[[boRepNo]]
      ### true best pair
      gBestTrue <- max(gAScaled[initPairNo] + gBScaled[initPairNo])
      bestPairNo <- which.max((gAScaled + gBScaled) == gBestTrue)
      bestPairName <- pairName[bestPairNo]
      names(gBestTrue) <- bestPairName
      ### resample if bestPairNo is included
      tryNo <- 0
      while (bestPairNo %in% initPairNo) {
        tryNo <- tryNo + 1
        set.seed(seed[boRepNo] + tryNo)
        initPairNo <- sample(x = 1:nPair, size = nInitComb, replace = FALSE)
      }
      if (tryNo > 0) {
        initPairNoLst[[boRepNo]] <- initPairNo
        seed[boRepNo] <- seed[boRepNo] + tryNo
        initPairUpdated <- TRUE
      }
      
      return(gBestTrue)
    }))
  })
  
  ### save updated initPairNoLst
  if (initPairUpdated) {
    saveRDS(object = initPairNoLst, file = fileNameInitPairNo)
    saveRDS(object = seed, file = fileNameSeedBo)
  }
  ### save initial best pair
  saveRDS(object = bestInitPairLst, file = fileNameBoInitPair)
}




###### 3. Perform BO ######
# ##### 3.1. Use Model 3 #####
# ### perform BO (change initial value for nRepBo times)
# for (boRepNo in 1:nRepBo) {
#   cat("---------- BORep: ", boRepNo, " / ", nRepBo, " ----------\n")
#   ### multicore for parameter
#   PerformBOM3 <- function(parNo) {
#     # parNo <- 1
#     ### scenario parameter
#     genCor <- parSetDf$genCor[parNo]
#     h2 <- parSetDf$h2[parNo]
#     gamma <- parSetDf$gamma[parNo]
#     repNo <- parSetDf$rep[parNo]
#     ### for saving BO result
#     ##############'#####################
#     dirNow <- paste0(saveDirM3, "parNo_", parNo, "/UCB/")
#     ##############'#####################
#     dirBORes <- paste0(dirNow, "BO_result/")
#     if (!dir.exists(dirBORes)) {
#       dir.create(dirBORes, recursive = TRUE)
#     }
#     ### load phenotype
#     fileNamePheno <- paste0(phenoDir, "0.3.allSampledRes_phenotype_genCor_", genCor,
#                             "_h2_", h2, "_gamma_", gamma, "_nRep_", nRep, ".rds")
#     phenoLst <- readRDS(file = fileNamePheno)[[repNo]]
#     ### genotypic value and scaling
#     gA <- phenoLst$dA + phenoLst$sBtoA + phenoLst$vA
#     gB <- phenoLst$dB + phenoLst$sAtoB + phenoLst$vB
#     gAScaled <- scale(gA)
#     gBScaled <- scale(gB)
#     
#     ### true best pair in GV
#     gBestTrue <- max(gAScaled + gBScaled, na.rm = TRUE)
#     ### for observed phenotype (update by each iteration)
#     yObsSb <- rep(NA, nPair)
#     yObsGm <- rep(NA, nPair)
#     yObsPair <- rep(NA, nPair)
#     gObsPair <- rep(NA, nPair) ### for true genotypic value
#     names(yObsSb) <- names(yObsGm) <- names(yObsPair) <- names(gObsPair) <- 
#       pairName
#     
#     ### save for each iteration's best individual
#     gBestAllIter <- NULL
#     selPairAllIterLst <- list()
#     
#     for (iterNo in 1:nIterBo) {
#       
#       cat("Iteration: ", iterNo, "\n")
#       
#       if (iterNo == 1) {
#         ### pick the initial observe individuals
#         trainPairNo <- initPairNoLst[[boRepNo]]
#         ### save observe pair number
#         selPairAllIterLst <- c(selPairAllIterLst, list(trainPairNo))
#       } else {
#         ### pick the next observe individuals
#         ##############'#####################
#         newObsDf <- evalDf[order(evalDf$UCB, decreasing = TRUE)[1:nPickComb], ]
#         ##############'#####################
#         trainPairNo <- c(trainPairNo, newObsDf$ID)
#         ### save observe pair number
#         selPairAllIterLst <- c(selPairAllIterLst, list(newObsDf$ID))
#       }
#       ### observe individuals
#       yObsSb[trainPairNo] <- phenoLst$scaled_yA[trainPairNo]
#       yObsGm[trainPairNo] <- phenoLst$scaled_yB[trainPairNo]
#       yObsPair[trainPairNo] <- phenoLst$scaled_yA[trainPairNo] + 
#         phenoLst$scaled_yB[trainPairNo]
#       ### save current best pair in GV
#       gObsPair[trainPairNo] <- gAScaled[trainPairNo] + gBScaled[trainPairNo]
#       gBestNow <- sort(gObsPair, decreasing = TRUE, na.last = TRUE)[1]
#       gBestAllIter <- c(gBestAllIter, gBestNow)
#       ### convergence rule
#       if (gBestNow == gBestTrue) {
#         break
#       }
#       ### define BGLR save directory
#       dirBGLRPar <- paste0(dirNow, "BGLR_parameter/")
#       if (!dir.exists(dirBGLRPar)) {
#         dir.create(dirBGLRPar, recursive = TRUE)
#       }
#       dirBGLRRes <- paste0(dirNow, "BGLR_result/")
#       if (!dir.exists(dirBGLRRes)) {
#         dir.create(dirBGLRRes, recursive = TRUE)
#       }
#       fileNameFitRes <- paste0(dirBGLRRes, scriptID, "fitting_result_BORepNo_",
#                                boRepNo, "_iterNo_", iterNo, ".rds")
#       if (!file.exists(fileNameFitRes)) {
#         ### fitting for all individuals
#         fitResLst <- lapply(X = 1:2, FUN = function(cropNo) {
#           ### for BGLR result save
#           fileNameParameter <- paste0(dirBGLRPar, scriptID, cropName[cropNo],
#                                       "_BORepNo_", boRepNo, "_iterNo_", iterNo, "_")
#           ### prepare phenotype (only observed)
#           if (cropNo == 1) {
#             y <- yObsSb[trainPairNo]
#           } else {
#             y <- yObsGm[trainPairNo]
#           }
#           ### BGLR
#           ETA <-  list(Sb = list(K = K1[trainPairNo, trainPairNo],
#                                  model = "RKHS"),
#                        Gm = list(K = K2[trainPairNo, trainPairNo],
#                                  model = "RKHS"),
#                        GbyG = list(K = K3[trainPairNo, trainPairNo],
#                                    model = "RKHS"))
#           fit <- BGLR::BGLR(y = y, ETA = ETA, nIter = nIter, burnIn = burnIn,
#                             saveAt = fileNameParameter, verbose = FALSE)
#           ### value for calculate acquisition function
#           CalcPostMeanVar(yObs = y, obs_idx = trainPairNo,
#                           mu_hat = fit$mu, sigma2_e = fit$varE,
#                           K_A = K1, K_B = K2, K_AB = K3,
#                           sigma2_A = fit$ETA$Sb$varU,
#                           sigma2_B = fit$ETA$Gm$varU,
#                           sigma2_AB = fit$ETA$GbyG$varU)
#           
#         })
#         names(fitResLst) <- cropName
#         saveRDS(object = fitResLst, file = fileNameFitRes)
#       } else {
#         fitResLst <- readRDS(file = fileNameFitRes)
#       }
#       ### value for calculate acquisition function
#       muCand <- fitResLst$Sb$mean_A_B_AB + fitResLst$Gm$mean_A_B_AB
#       sdCand <- fitResLst$Sb$sd_A_B_AB + fitResLst$Gm$sd_A_B_AB
#       names(muCand) <- names(sdCand) <- pairName[-trainPairNo]
#       ### candidate name for BO evaluation
#       candIdx <- which(is.na(yObsPair))
#       ### get acquisition function
#       evalDf <- cal.AF.BGLR(gBest = gBestNow, candIdx = candIdx,
#                             muCand = muCand, sdCand = sdCand, kappa = kappa)
#       gc(reset = TRUE);gc(reset = TRUE)
#     }
#     
#     ### save all iteration BO result
#     fileNameBOyBest <- paste0(dirBORes, scriptID, "g_best_allIteration_BORepNo_",
#                               boRepNo, ".rds")
#     saveRDS(object = gBestAllIter, file = fileNameBOyBest)
#     fileNameBOSelPair <- paste0(dirBORes, scriptID, "selectedPairNo_allIteration_BORepNo_",
#                                 boRepNo, ".rds")
#     saveRDS(object = selPairAllIterLst, file = fileNameBOSelPair)
#   }
#   pbmcapply::pbmclapply(X = 1:nPar, FUN = PerformBOM3, mc.cores = nCore)
# }




##### 3.2. Use Model 5 #####
### perform BO (change initial value for nRepBo times)
for (boRepNo in 1:nRepBo) {
  cat("---------- BORep: ", boRepNo, " / ", nRepBo, " ----------\n")
  ### multicore for parameter
  PerformBOM5 <- function(parNo) {
    # parNo <- 1
    ### scenario parameter
    genCor <- parSetDf$genCor[parNo]
    h2 <- parSetDf$h2[parNo]
    gamma <- parSetDf$gamma[parNo]
    repNo <- parSetDf$rep[parNo]
    ### for saving BO result
    ##############'#####################
    dirNow <- paste0(saveDirM5, "parNo_", parNo, "/Mu/")
    ##############'#####################
    dirBORes <- paste0(dirNow, "BO_result/")
    if (!dir.exists(dirBORes)) {
      dir.create(dirBORes, recursive = TRUE)
    }
    ### load phenotype
    fileNamePheno <- paste0(phenoDir, "0.3.allSampledRes_phenotype_genCor_", genCor,
                            "_h2_", h2, "_gamma_", gamma, "_nRep_", nRep, ".rds")
    phenoLst <- readRDS(file = fileNamePheno)[[repNo]]
    ### genotypic value and scaling
    gA <- phenoLst$dA + phenoLst$sBtoA + phenoLst$vA
    gB <- phenoLst$dB + phenoLst$sAtoB + phenoLst$vB
    gAScaled <- scale(gA)
    gBScaled <- scale(gB)
    
    ### true best pair in GV
    gBestTrue <- max(gAScaled + gBScaled, na.rm = TRUE)
    ### for observed phenotype (update by each iteration)
    yObsSb <- rep(NA, nPair)
    yObsGm <- rep(NA, nPair)
    yObsPair <- rep(NA, nPair)
    gObsPair <- rep(NA, nPair) ### for true genotypic value
    names(yObsSb) <- names(yObsGm) <- names(yObsPair) <- names(gObsPair) <- 
      pairName
    
    ### save for each iteration's best individual
    gBestAllIter <- NULL
    selPairAllIterLst <- list()
    
    for (iterNo in 1:nIterBo) {
      
      cat("Iteration: ", iterNo, "\n")
      
      if (iterNo == 1) {
        ### pick the initial observe individuals
        trainPairNo <- initPairNoLst[[boRepNo]]
        ### save observe pair number
        selPairAllIterLst <- c(selPairAllIterLst, list(trainPairNo))
      } else {
        ### pick the next observe individuals
        ##############'#####################
        newObsDf <- evalDf[order(evalDf$Mu, decreasing = TRUE)[1:nPickComb], ]
        ##############'#####################
        trainPairNo <- c(trainPairNo, newObsDf$ID)
        ### save observe pair number
        selPairAllIterLst <- c(selPairAllIterLst, list(newObsDf$ID))
      }
      ## observe individuals
      yObsSb[trainPairNo] <- phenoLst$scaled_yA[trainPairNo]
      yObsGm[trainPairNo] <- phenoLst$scaled_yB[trainPairNo]
      yObsPair[trainPairNo] <- phenoLst$scaled_yA[trainPairNo] +
        phenoLst$scaled_yB[trainPairNo]
      ### save current best pair in GV
      gObsPair[trainPairNo] <- gAScaled[trainPairNo] + gBScaled[trainPairNo]
      gBestNow <- sort(gObsPair, decreasing = TRUE, na.last = TRUE)[1]
      gBestAllIter <- c(gBestAllIter, gBestNow)
      ### convergence rule
      if (gBestNow == gBestTrue) {
        break
      }
      ### define BGLR save directory
      dirBGLRPar <- paste0(dirNow, "BGLR_parameter/")
      if (!dir.exists(dirBGLRPar)) {
        dir.create(dirBGLRPar, recursive = TRUE)
      }
      dirBGLRRes <- paste0(dirNow, "BGLR_result/")
      if (!dir.exists(dirBGLRRes)) {
        dir.create(dirBGLRRes, recursive = TRUE)
      }
      fileNameFitRes <- paste0(dirBGLRRes, scriptID, "fitting_result_BORepNo_",
                               boRepNo, "_iterNo_", iterNo, ".rds")
      if (!file.exists(fileNameFitRes)) {
        ### fitting for all individuals
        ### for BGLR result save
        fileNameParameter <- paste0(dirBGLRPar, scriptID, "Multitrait_BORepNo_",
                                    boRepNo, "_iterNo_", iterNo, "_")
        ### prepare phenotype (only observed)
        yNow <- cbind(yObsSb[trainPairNo], yObsGm[trainPairNo])
        ### for model
        ETA <- list(GenoSb = list(K = K1[trainPairNo, trainPairNo],
                                  model = "RKHS", Cov = list(type = "UN")),
                    GenoGm = list(K = K2[trainPairNo, trainPairNo],
                                  model = "RKHS", Cov = list(type = "UN")),
                    GenoGbyG = list(K = K3[trainPairNo, trainPairNo],
                                    model = "RKHS", Cov = list(type = "DIAG")))
        resCov <- list(S0 = diag(2), type = "UN")
        ### fitting
        fitResM5 <- Multitrait(y = yNow, resCov = resCov,
                               saveAt = fileNameParameter, ETA = ETA,
                               nIter = nIter, burnIn = burnIn, thin = thin,
                               verbose = FALSE)
        
        
        ### value for calculate acquisition function
        resLst <- CalcPostMeanVar_MT(yObs = yNow, obs_idx = trainPairNo,
                                     mu_hat = fitResM5$mu, K_A = K1, K_B = K2, K_AB = K3,
                                     G1 = fitResM5$ETA$GenoSb$Cov$Omega,
                                     G2 = fitResM5$ETA$GenoGm$Cov$Omega,
                                     G3 = fitResM5$ETA$GenoGbyG$Cov$Omega,
                                     R0 = fitResM5$resCov$R)
        ### save results of unobserved pairs' posterior
        saveRDS(object = resLst, file = fileNameFitRes)
      } else {
        resLst <- readRDS(file = fileNameFitRes)
      }
      
      ### value for calculate acquisition function
      muCand <- resLst$mean_uSum
      sdCand <- sqrt(resLst$cov_uSum)
      names(muCand) <- names(sdCand) <- pairName[-trainPairNo]
      ### candidate name for BO evaluation
      candIdx <- which(is.na(yObsPair))
      ### get acquisition function
      evalDf <- cal.AF.BGLR(gBest = gBestNow, candIdx = candIdx,
                            muCand = muCand, sdCand = sdCand, kappa = kappa)
      gc(reset = TRUE);gc(reset = TRUE)
    }
    
    ### save all iteration BO result
    fileNameBOyBest <- paste0(dirBORes, scriptID, "y_best_allIteration_BORepNo_",
                              boRepNo, ".rds")
    saveRDS(object = yBestAllIter, file = fileNameBOyBest)
    fileNameBOSelPair <- paste0(dirBORes, scriptID, "selectedPairNo_allIteration_BORepNo_",
                                boRepNo, ".rds")
    saveRDS(object = selPairAllIterLst, file = fileNameBOSelPair)
  }
  pbmcapply::pbmclapply(X = 1:nPar, FUN = PerformBOM5, mc.cores = 20)
}




# ###### 4. Random Search ######
# ### set seed for initial pairs of BO
# fileNameSeedRs <- paste0(saveDirRandom, scriptID, "seed_RandomSearch_AllIterPairNo.rds")
# if (!file.exists(fileNameSeedRs)) {
#   seedRs <- sample(x = 1e06, size = nRepBo, replace = FALSE)
#   saveRDS(object = seedRs, file = fileNameSeedRs)
# } else {
#   seedRs <- readRDS(file = fileNameSeedRs)
# }
# 
# ##### 4.1. PairSet for Each Cycle #####
# nIterRandom <- ceiling(nPair / nPickComb)
# randomIterNoLst <- lapply(X = 1:nRepBo, FUN = function(boRepNo) {
#   initPairNo <- initPairNoLst[[boRepNo]]
#   pairNoRemain <- (1:nPair)[-initPairNo]
#   ### randomize
#   set.seed(seedRs[boRepNo])
#   pairNoSamp <- sample(pairNoRemain)
#   ### divide
#   idxMid <- head(pairNoSamp, nPickComb * (nIterRandom - 2))
#   idxLast <- tail(pairNoSamp, (nPair %% nPickComb))
#   
#   idxRes <- 1:nPair
#   idxRes[initPairNo] <- 1 ### initial pairs
#   idxRes[idxMid] <- rep(2:(nIterRandom - 1), each = nPickComb) ### middle iteration pair
#   idxRes[idxLast] <- nIterRandom ### last pairs
#   
#   return(idxRes)
# })
# 
# ##### 4.2. Perform Random Search #####
# for (boRepNo in 1:nRepBo) {
#   cat("---------- BORep: ", boRepNo, " / ", nRepBo, " ----------\n")
#   ### multicore for parameter
#   PerformRS <- function(parNo) {
#     # parNo <- 1
#     ### scenario parameter
#     genCor <- parSetDf$genCor[parNo]
#     h2 <- parSetDf$h2[parNo]
#     gamma <- parSetDf$gamma[parNo]
#     repNo <- parSetDf$rep[parNo]
#     ### for saving BO result
#     dirNow <- paste0(saveDirRandom, "parNo_", parNo, "/RandomSearch/")
#     dirBORes <- paste0(dirNow, "BO_result/")
#     if (!dir.exists(dirBORes)) {
#       dir.create(dirBORes, recursive = TRUE)
#     }
#     ### load phenotype
#     fileNamePheno <- paste0(phenoDir, "0.3.allSampledRes_phenotype_genCor_", genCor,
#                             "_h2_", h2, "_gamma_", gamma, "_nRep_", nRep, ".rds")
#     phenoLst <- readRDS(file = fileNamePheno)[[repNo]]
#     
#     ### true best pair
#     yBestTrue <- max(phenoLst$yA + phenoLst$yB, na.rm = TRUE)
#     ### for observed phenotype (update by each iteration)
#     yObsSb <- rep(NA, nPair)
#     yObsGm <- rep(NA, nPair)
#     yObsPair <- rep(NA, nPair)
#     names(yObsSb) <- names(yObsGm) <- names(yObsPair) <- pairName
#     
#     ### save for each iteration's best individual
#     yBestAllIter <- NULL
#     selPairAllIterLst <- list()
#     
#     for (iterNo in 1:nIterRandom) {
#       
#       cat("Iteration: ", iterNo, "\n")
#       
#       if (iterNo == 1) {
#         ### pick the initial observe individuals
#         trainPairNo <- which(randomIterNoLst[[boRepNo]] == iterNo)
#         ### save observe pair number
#         selPairAllIterLst <- c(selPairAllIterLst, list(trainPairNo))
#       } else {
#         ### pick the next observe individuals
#         newTrainPairNo <- which(randomIterNoLst[[boRepNo]] == iterNo)
#         trainPairNo <- c(trainPairNo, newTrainPairNo)
#         ### save observe pair number
#         selPairAllIterLst <- c(selPairAllIterLst, list(newTrainPairNo))
#       }
#       ### observe individuals
#       yObsSb[trainPairNo] <- phenoLst$yA[trainPairNo]
#       yObsGm[trainPairNo] <- phenoLst$yB[trainPairNo]
#       yObsPair[trainPairNo] <- phenoLst$yA[trainPairNo] + phenoLst$yB[trainPairNo]
#       ### save current best pair
#       yBestNow <- sort(yObsPair, decreasing = TRUE, na.last = TRUE)[1]
#       yBestAllIter <- c(yBestAllIter, yBestNow)
#       ### convergence rule
#       if (yBestNow == yBestTrue) {
#         break
#       }
#     }
#     
#     ### save all iteration BO result
#     fileNameBOyBest <- paste0(dirBORes, scriptID, "y_best_allIteration_BORepNo_",
#                               boRepNo, ".rds")
#     saveRDS(object = yBestAllIter, file = fileNameBOyBest)
#     fileNameBOSelPair <- paste0(dirBORes, scriptID, "selectedPairNo_allIteration_BORepNo_",
#                                 boRepNo, ".rds")
#     saveRDS(object = selPairAllIterLst, file = fileNameBOSelPair)
#   }
#   pbmcapply::pbmclapply(X = 1:nPar, FUN = PerformRS, mc.cores = nCore)
# }
