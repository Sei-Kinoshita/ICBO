###### 1. Settings ######
##### 1.1. Remove objects #####
rm(list = ls())
gc(reset = TRUE);gc(reset = TRUE)

##### 1.2. Load Packages and Functions #####
require(BGLR)
source("scripts/function_scripts/2.0.function.BO.R")

##### 1.3. Define Parameters #####
scriptID <- "2.0."
indexName <- c("Mu", "PI", "EI", "UCB")
essThreshold <- 400
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

nCore <- 40

##### 1.4. Define Directories #####
dataDir <- "data/"
genoDir <- paste0(dataDir, "genotype/")
phenoDir <- paste0(dataDir, "phenotype/")

saveDirM3 <- paste0("midstream/2.0.BO/Model3_GMA_SMA/")
if (!dir.exists(saveDirM3)) {
  dir.create(saveDirM3, recursive = TRUE)
}

###### 2. Load Data ######
### K mat (ZGZ^T)
fileNameKSb <- "midstream/1.0.GP_Fitting/1.0.Sb_KMat.rds"
fileNameKGm <- "midstream/1.0.GP_Fitting/1.0.Gm_KMat.rds"
fileNameKGbyG <- "midstream/1.0.GP_Fitting/1.0.GbyG_KMat.rds"
K1 <- readRDS(file = fileNameKSb)
K2 <- readRDS(file = fileNameKGm)
K3 <- readRDS(file = fileNameKGbyG)
kernelInputInfo <- file.info(c(fileNameKSb, fileNameKGm, fileNameKGbyG))[, c("size", "mtime")]

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
parNoVec <- seq_len(nPar)

### Existing Kronecker GRM; distances are computed only for each selected batch.
grmPair <- K3

###### 3. Prepare for BO ######
### number of iterations for BO
nIterBo <- floor((nPair - nInitComb) / nPickComb)

### set seed for initial pairs of BO
fileNameSeedBo <- paste0(saveDirM3, scriptID, "seed_BO_initPair.rds")
if (!file.exists(fileNameSeedBo)) {
  set.seed(20260928)
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

### Exclude the scenario's true optimum using one common initial batch per run.
fileNameBoInitPair <- paste0(saveDirM3, scriptID, "BO_initPair_bestValue.rds")
if (!file.exists(fileNameBoInitPair) || is.null(attr(initPairNoLst, "byScenario"))) {
  initPairByScenario <- vector("list", nPar)
  bestInitPairLst <- lapply(seq_len(nRepBo), function(i) rep(NA_real_, nPar))
  for (parNo in parNoVec) {
    genCor <- parSetDf$genCor[parNo]
    h2 <- parSetDf$h2[parNo]
    gamma <- parSetDf$gamma[parNo]
    repNo <- parSetDf$rep[parNo]
    fileNamePheno <- paste0(phenoDir, "0.3.",
                            "allSampledRes_phenotype_genCor_", genCor,
                            "_h2_", h2, "_gamma_", gamma, "_nRep_", nRep, ".rds")
    phenoLst <- readRDS(fileNamePheno)[[repNo]]
    gTrue <- as.numeric(scale(phenoLst$dA + phenoLst$sBtoA + phenoLst$vA) +
                          scale(phenoLst$dB + phenoLst$sAtoB + phenoLst$vB))
    bestPairNo <- which.max(gTrue)
    initPairByScenario[[parNo]] <- vector("list", nRepBo)
    for (boRepNo in seq_len(nRepBo)) {
      initPairNo <- initPairNoLst[[boRepNo]]
      tryNo <- 0L
      while (bestPairNo %in% initPairNo) {
        tryNo <- tryNo + 1L
        set.seed(seed[boRepNo] + tryNo)
        initPairNo <- sample(seq_len(nPair), nInitComb, replace = FALSE)
      }
      initPairByScenario[[parNo]][[boRepNo]] <- initPairNo
      bestInitPairLst[[boRepNo]][parNo] <- max(gTrue[initPairNo])
    }
  }
  attr(initPairNoLst, "byScenario") <- initPairByScenario
  saveRDS(initPairNoLst, fileNameInitPairNo)
  saveRDS(bestInitPairLst, fileNameBoInitPair)
}

###### 4. Perform BO ######
##### 4.1. Use Model 3 #####
### perform BO (change initial value for nRepBo times)
for (indexNow in indexName) {
for (boRepNo in 1:nRepBo) {
  cat("---------- BORep: ", boRepNo, " / ", nRepBo, " ----------\n")
  ### multicore for parameter
  PerformBOM3 <- function(parNo) {
    runStart <- proc.time()[["elapsed"]]
    set.seed(seed[boRepNo] + 1000L * parNo)
    ### scenario parameter
    genCor <- parSetDf$genCor[parNo]
    h2 <- parSetDf$h2[parNo]
    gamma <- parSetDf$gamma[parNo]
    repNo <- parSetDf$rep[parNo]
    ### for saving BO result
    dirNow <- paste0(saveDirM3, "parNo_", parNo, "/", indexNow, "/")
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
    batchDiversityLst <- cycleTimingLst <- mcmcDiagnosticsLst <- list()
    
    for (iterNo in 1:nIterBo) {
      cycleStart <- proc.time()[["elapsed"]]
      selectionStart <- cycleStart
      
      cat("Iteration: ", iterNo, "\n")
      
      if (iterNo == 1) {
        ### pick the initial observe individuals
        trainPairNo <- attr(initPairNoLst, "byScenario")[[parNo]][[boRepNo]]
        ### save observe pair number
        selPairAllIterLst <- c(selPairAllIterLst, list(trainPairNo))
      } else {
        ### pick the next observe individuals
        newObsDf <- evalDf[order(evalDf[[indexNow]], decreasing = TRUE)[1:nPickComb], ]
        trainPairNo <- c(trainPairNo, newObsDf$ID)
        ### save observe pair number
        selPairAllIterLst <- c(selPairAllIterLst, list(newObsDf$ID))
      }
      batchSelectionTime <- proc.time()[["elapsed"]] - selectionStart
      batchDiversityLst[[iterNo]] <- cbind(cycle = iterNo,
        BatchDiversity(tail(selPairAllIterLst, 1)[[1]], pairDf, grmPair))
      ### observe individuals
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
        cycleTimingLst[[iterNo]] <- data.frame(cycle = iterNo, nEvaluatedPair = length(trainPairNo),
          gpFittingPredictionTime = 0, acquisitionTime = 0, batchSelectionTime = batchSelectionTime,
          totalCycleTime = proc.time()[["elapsed"]] - cycleStart, fitReused = FALSE, gpFitPerformed = FALSE)
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
      ### Verify the inputs and restore the RNG stream when reusing a completed fit.
      fitContext <- list(trainingPairNo = trainPairNo,
        ySb = yObsSb[trainPairNo], yGm = yObsGm[trainPairNo],
        nIter = nIter, burnIn = burnIn, thin = thin,
        objectiveSDDefinition = "sqrt_sum_crop_variances",
        kernelInputInfo = kernelInputInfo,
        rngStateBeforeFit = get(".Random.seed", envir = .GlobalEnv))
      if (indexNow %in% c("PI", "EI")) {
        fitContext$incumbentDefinition <- "max_observed_posterior_mean"
      }
      fitReused <- file.exists(fileNameFitRes)
      if (fitReused) {
        fitResLst <- readRDS(file = fileNameFitRes)
        if (!identical(attr(fitResLst, "fitContext"), fitContext) ||
            is.null(attr(fitResLst, "rngStateAfterFit"))) {
          stop("Cached BO fit is incompatible or lacks restart metadata: ", fileNameFitRes,
               ". Recompute this BO run after removing its incompatible cached fit files.")
        }
      }
      if (!fitReused) {
        ### fitting for all individuals
        fitResLst <- lapply(X = 1:2, FUN = function(cropNo) {
          ### for BGLR result save
          fileNameParameter <- paste0(dirBGLRPar, scriptID, cropName[cropNo],
                                      "_BORepNo_", boRepNo, "_iterNo_", iterNo, "_")
          ### prepare phenotype (only observed)
          if (cropNo == 1) {
            y <- yObsSb[trainPairNo]
          } else {
            y <- yObsGm[trainPairNo]
          }
          ### BGLR
          ETA <-  list(Sb = list(K = K1[trainPairNo, trainPairNo],
                                 model = "RKHS"),
                       Gm = list(K = K2[trainPairNo, trainPairNo],
                                 model = "RKHS"),
                       GbyG = list(K = K3[trainPairNo, trainPairNo],
                                   model = "RKHS"))
          fitStart <- proc.time()[["elapsed"]]
          fit <- BGLR::BGLR(y = y, ETA = ETA, nIter = nIter, burnIn = burnIn, thin = thin, 
                            saveAt = fileNameParameter, verbose = FALSE)
          ### value for calculate acquisition function
          postRes <- CalcPostMeanVar(yObs = y, obs_idx = trainPairNo,
                          mu_hat = fit$mu, sigma2_e = fit$varE,
                          K_A = K1, K_B = K2, K_AB = K3,
                          sigma2_A = fit$ETA$Sb$varU,
                          sigma2_B = fit$ETA$Gm$varU,
                          sigma2_AB = fit$ETA$GbyG$varU,
                          predictObserved = indexNow %in% c("PI", "EI"))
          postRes$fitPredictionTime <- proc.time()[["elapsed"]] - fitStart
          postRes$mcmcDiagnostics <- SingleChainDiagnostics(fileNameParameter, cropName[cropNo],
            burnIn = burnIn, thin = fit$thin, essThreshold = essThreshold)
          postRes$mcmcSettings <- list(nIter = fit$nIter, burnIn = fit$burnIn, thin = fit$thin)
          return(postRes)
          
        })
        names(fitResLst) <- cropName
        attr(fitResLst, "fitContext") <- fitContext
        attr(fitResLst, "rngStateAfterFit") <- get(".Random.seed", envir = .GlobalEnv)
        ### Publish only a complete checkpoint; an interrupted write is not reused.
        saveRDS(object = fitResLst, file = paste0(fileNameFitRes, ".tmp"))
        if (!file.rename(paste0(fileNameFitRes, ".tmp"), fileNameFitRes)) {
          stop("Could not finalize the BO checkpoint: ", fileNameFitRes)
        }
      } else {
        assign(".Random.seed", attr(fitResLst, "rngStateAfterFit"), envir = .GlobalEnv)
      }
      mcmcDiagnosticsLst[[iterNo]] <- cbind(cycle = iterNo,
        do.call(rbind, lapply(fitResLst, function(x) x$mcmcDiagnostics)))
      acquisitionStart <- proc.time()[["elapsed"]]
      ### value for calculate acquisition function
      muCand <- fitResLst$Sb$mean_A_B_AB + fitResLst$Gm$mean_A_B_AB
      ### Model 3 fits crops independently; f = g_Sb + g_Gm (manuscript Sections 2.4 and 2.6).
      sdCand <- sqrt(fitResLst$Sb$sd_A_B_AB^2 + fitResLst$Gm$sd_A_B_AB^2)
      names(muCand) <- names(sdCand) <- pairName[-trainPairNo]
      ### candidate name for BO evaluation
      candIdx <- which(is.na(yObsPair))
      ### EI/PI use the current fit's latent objective means over observed pairs only.
      ### One incumbent is shared by the next batch; recompute after its observations enter the fit.
      fBest <- NA_real_
      if (indexNow %in% c("PI", "EI")) {
        muObserved <- fitResLst$Sb$meanObserved_A_B_AB + fitResLst$Gm$meanObserved_A_B_AB
        stopifnot(length(muObserved) == length(trainPairNo), all(is.finite(muObserved)))
        fBest <- max(muObserved)
      }
      ### get acquisition function
      evalDf <- cal.AF.BGLR(gBest = fBest, candIdx = candIdx,
                            muCand = muCand, sdCand = sdCand, kappa = kappa)
      acquisitionTime <- proc.time()[["elapsed"]] - acquisitionStart
      cycleTimingLst[[iterNo]] <- data.frame(cycle = iterNo, nEvaluatedPair = length(trainPairNo),
        gpFittingPredictionTime = if (fitReused) 0 else sum(vapply(fitResLst, function(x) x$fitPredictionTime, numeric(1))),
        acquisitionTime = acquisitionTime, batchSelectionTime = batchSelectionTime,
        totalCycleTime = proc.time()[["elapsed"]] - cycleStart, fitReused = fitReused, gpFitPerformed = !fitReused)
      gc(reset = TRUE);gc(reset = TRUE)
      cycleTimingLst[[iterNo]]$totalCycleTime <- proc.time()[["elapsed"]] - cycleStart
    }
    
    attr(selPairAllIterLst, "batchDiversity") <- do.call(rbind, batchDiversityLst)
    attr(selPairAllIterLst, "timing") <- do.call(rbind, cycleTimingLst)
    attr(selPairAllIterLst, "mcmcDiagnostics") <- do.call(rbind, mcmcDiagnosticsLst)
    attr(selPairAllIterLst, "runTiming") <- data.frame(genCor = genCor, h2 = h2, gamma = gamma,
      phenotypeRep = repNo, method = indexNow, boRep = boRepNo,
      totalOptimizationTime = proc.time()[["elapsed"]] - runStart,
      optimumFound = gBestNow == gBestTrue)
    ### save all iteration BO result
    fileNameBOyBest <- paste0(dirBORes, scriptID, "g_best_allIteration_BORepNo_",
                              boRepNo, ".rds")
    saveRDS(object = gBestAllIter, file = fileNameBOyBest)
    fileNameBOSelPair <- paste0(dirBORes, scriptID, "selectedPairNo_allIteration_BORepNo_",
                                boRepNo, ".rds")
    saveRDS(object = selPairAllIterLst, file = fileNameBOSelPair)
  }
  pbmcapply::pbmclapply(X = parNoVec, FUN = PerformBOM3, mc.cores = nCore)
}
}
