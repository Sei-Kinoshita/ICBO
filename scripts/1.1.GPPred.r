###### 0. Settings ######
##### 0.1. Remove objects #####
rm(list = ls())
gc(reset = TRUE);gc(reset = TRUE)

##### 0.2. Load packages #####
require(BGLR)

##### 0.3. Set parameters #####
scriptID <- "1.1."

cropName <- c("Sb", "Gm")

nRep <- 10 ### repetition for phenotype

# nCores <- 30
nCores <- 40
nCV <- 10
### for BGLR
nIter <- 10000
burnIn <- 2000
thin <- 20
# nIter <- 10
# burnIn <- 2
# thin <- 1

##### 0.4. Define directories #####
dirData <- "data/"
dirGeno <- paste0(dirData, "genotype/")
dirPheno <- paste0(dirData, "phenotype/")
dirComb <- paste0(dirData, "selTrainData_GP/")
dirFit <- "midstream/1.0.GP_Fitting/"

dirScr <- paste0("midstream/", scriptID, "GP_CV/") 
if (!dir.exists(dirScr)) {
  dir.create(dirScr, recursive = TRUE)
}
#### 0.4.1. Model1 directories ####
dirM1 <- paste0(dirScr, "Model1_DGE/")
if (!dir.exists(dirM1)) {
  dir.create(dirM1, recursive = TRUE)
}
#### 0.4.2. Model2 directories ####
dirM2 <- paste0(dirScr, "Model2_DGE_IGE/")
if (!dir.exists(dirM2)) {
  dir.create(dirM2, recursive = TRUE)
}
#### 0.4.3. Model3 directories ####
dirM3 <- paste0(dirScr, "Model3_GMA_SMA/")
if (!dir.exists(dirM3)) {
  dir.create(dirM3, recursive = TRUE)
}
#### 0.4.4. Model4 directories ####
dirM4 <- paste0(dirScr, "Model4_Multi_DGE_IGE/")
if (!dir.exists(dirM4)) {
  dir.create(dirM4, recursive = TRUE)
}
#### 0.4.5. Model5 directories ####
dirM5 <- paste0(dirScr, "Model5_Multi_GMA_SMA/")
if (!dir.exists(dirM5)) {
  dir.create(dirM5, recursive = TRUE)
}



###### 1. Load Data ######
##### 1.1. Parameter Set #####
fileNameParSet <- paste0(dirPheno, "0.3.parmeterSet.rds")
parSetDf <- readRDS(file = fileNameParSet)
nPar <- nrow(parSetDf)

##### 1.2. Pair Name (fixed order) #####
fileNamePairName <- paste0(dirPheno, "0.3.pairNameDf_5214Pair.rds")
pairNameDf <- readRDS(file = fileNamePairName)

##### 1.3. K matrix #####
fileNameK1 <- paste0(dirFit, "1.0.Sb_KMat.rds")
fileNameK2 <- paste0(dirFit, "1.0.Gm_KMat.rds")
fileNameK3 <- paste0(dirFit, "1.0.GbyG_KMat.rds")
K1 <- readRDS(file = fileNameK1)
K2 <- readRDS(file = fileNameK2)
K3 <- readRDS(file = fileNameK3)

##### 1.4. Selected combination for CV #####
fileNameSelCombLst1 <- paste0(dirComb, "0.4.trainData_selCombName_random.rds")
fileNameSelCombLst2 <- paste0(dirComb, "0.4.trainData_selCombName_allComb.rds")
fileNameSelCombLst3 <- paste0(dirComb, "0.4.trainData_selCombName_useOnce.rds")
selCombRandomLst <- readRDS(file = fileNameSelCombLst1)
selCombAllLst <- readRDS(file = fileNameSelCombLst2)
selCombOnceLst <- readRDS(file = fileNameSelCombLst3)

cvCombLst <- list(Random = selCombRandomLst, All = selCombAllLst, Once = selCombOnceLst)
cvMethod <- names(cvCombLst)
rm(selCombRandomLst);rm(selCombAllLst);rm(selCombOnceLst)




###### 2. Model1 (only DGE) ######
##### 2.1. Prediction #####
for (cvMethodNo in 1:3) {
  for (cvNo in 1:nCV) {
    cat("---------- Model: 1 / CV:", cvNo," / cvMethod:", cvMethod[cvMethodNo], "----------\n")

    cvCombLstNow <- cvCombLst[[cvMethodNo]][[cvNo]]

    PerformM1 <- function(parNo) {

      ### parmeter set
      genCor <- parSetDf$genCor[parNo]
      h2 <- parSetDf$h2[parNo]
      gamma <- parSetDf$gamma[parNo]
      repNo <- parSetDf$rep[parNo]

      ### prepare saving directory
      ### model parameter by BGLR
      dirM1Par <- paste0(dirM1, cvMethod[cvMethodNo], "/CV", cvNo, "/modelPar/")
      if (!dir.exists(dirM1Par)) {
        dir.create(dirM1Par, recursive = TRUE)
      }
      ### all result by BGLR
      dirM1Fit <- paste0(dirM1, cvMethod[cvMethodNo], "/CV", cvNo, "/allRes/")
      if (!dir.exists(dirM1Fit)) {
        dir.create(dirM1Fit, recursive = TRUE)
      }

      ### load phenotype data
      fileNamePhenoINow <- paste0(dirPheno, "0.3.scaled_phenotype_genCor_", genCor,
                                  "_h2_", h2, "_gamma_", gamma, "_nRep_", nRep, ".rds")
      phenoLst <- readRDS(fileNamePhenoINow)

      ### fitting for each crop
      lapply(X = 1:2, FUN = function(cropNo) {
        ### intercropping phenotype for each crop
        phenoVecNow <- phenoLst[[repNo]][[cropNo]]
        ### for save
        fileNameParameter <- paste0(dirM1Par, scriptID, cropName[cropNo],
                                    "_parNo_", parNo, "_", cvMethod[cvMethodNo],
                                    "_CV_", cvNo, "_")
        fileNameFitting <- paste0(dirM1Fit, scriptID, cropName[cropNo],
                                  "_fitting_allRes_parNo_", parNo, "_",
                                  cvMethod[cvMethodNo], "_CV_", cvNo, ".rds")
        if (!file.exists(fileNameFitting)) {
          ### extract selected individuals only
          yNow <- rep(NA, length(phenoVecNow))
          names(yNow) <- names(phenoVecNow)
          ### training data for monocrop
          trainPos <- match(paste(cvCombLstNow$Sb, cvCombLstNow$Gm),
                            paste(pairNameDf$Sb, pairNameDf$Gm))
          yNow[trainPos] <- phenoVecNow[trainPos]
          ### for model
          if (cropNo == 1) {
            ETA <- list(GenoRKHS = list(K = K1, model = "RKHS"))
          } else {
            ETA <- list(GenoRKHS = list(K = K2, model = "RKHS"))
          }
          ### fitting
          fmResM1 <- BGLR(y = yNow,
                          saveAt = fileNameParameter, ETA = ETA,
                          nIter = nIter, burnIn = burnIn, thin = thin,
                          verbose = FALSE)
          ### save all fitting results of BayesB
          saveRDS(object = fmResM1, file = fileNameFitting)
        }
      })
    }
    pbmcapply::pbmclapply(X = 1:nPar, FUN = PerformM1, mc.cores = nCores)
    gc(reset = TRUE);gc(reset = TRUE)
  }
}




###### 3. Model2 (DGE + IGE) ######
##### 3.1. Prediction #####
for (cvMethodNo in 1:3) {
  for (cvNo in 1:nCV) {
    cat("---------- Model: 2 / CV:", cvNo," / cvMethod:", cvMethod[cvMethodNo], "----------\n")

    cvCombLstNow <- cvCombLst[[cvMethodNo]][[cvNo]]

    PerformM2 <- function(parNo) {

      ### parmeter set
      genCor <- parSetDf$genCor[parNo]
      h2 <- parSetDf$h2[parNo]
      gamma <- parSetDf$gamma[parNo]
      repNo <- parSetDf$rep[parNo]

      ### prepare saving directory
      ### model parameter by BGLR
      dirM2Par <- paste0(dirM2, cvMethod[cvMethodNo], "/CV", cvNo, "/modelPar/")
      if (!dir.exists(dirM2Par)) {
        dir.create(dirM2Par, recursive = TRUE)
      }
      ### all result by BGLR
      dirM2Fit <- paste0(dirM2, cvMethod[cvMethodNo], "/CV", cvNo, "/allRes/")
      if (!dir.exists(dirM2Fit)) {
        dir.create(dirM2Fit, recursive = TRUE)
      }

      ### load phenotype data
      fileNamePhenoINow <- paste0(dirPheno, "0.3.scaled_phenotype_genCor_", genCor,
                                  "_h2_", h2, "_gamma_", gamma, "_nRep_", nRep, ".rds")
      phenoLst <- readRDS(fileNamePhenoINow)

      ### fitting for each crop
      lapply(X = 1:2, FUN = function(cropNo) {
        ### intercropping phenotype for each crop
        phenoVecNow <- phenoLst[[repNo]][[cropNo]]
        ### for save
        fileNameParameter <- paste0(dirM2Par, scriptID, cropName[cropNo],
                                    "_parNo_", parNo, "_", cvMethod[cvMethodNo],
                                    "_CV_", cvNo, "_")
        fileNameFitting <- paste0(dirM2Fit, scriptID, cropName[cropNo],
                                  "_fitting_allRes_parNo_", parNo, "_",
                                  cvMethod[cvMethodNo], "_CV_", cvNo, ".rds")
        if (!file.exists(fileNameFitting)) {
          ### extract selected individuals only
          yNow <- rep(NA, length(phenoVecNow))
          names(yNow) <- names(phenoVecNow)
          ### training data for monocrop
          trainPos <- match(paste(cvCombLstNow$Sb, cvCombLstNow$Gm),
                            paste(pairNameDf$Sb, pairNameDf$Gm))
          yNow[trainPos] <- phenoVecNow[trainPos]
          ### for model
          ETA <- list(GenoSb = list(K = K1, model = "RKHS"),
                      GenoGm = list(K = K2, model = "RKHS"))
          ### fitting
          fmResM2 <- BGLR(y = yNow,
                          saveAt = fileNameParameter, ETA = ETA,
                          nIter = nIter, burnIn = burnIn, thin = thin,
                          verbose = FALSE)
          ### save all fitting results of BayesB
          saveRDS(object = fmResM2, file = fileNameFitting)
        }
      })
    }
    pbmcapply::pbmclapply(X = 1:nPar, FUN = PerformM2, mc.cores = nCores)
    gc(reset = TRUE);gc(reset = TRUE)
  }
}




###### 4. Model3 (DGE + IGE + SMA) ######
##### 4.1. Prediction #####
for (cvMethodNo in 1:3) {
  for (cvNo in 1:nCV) {
    cat("---------- Model: 3 / CV:", cvNo," / cvMethod:", cvMethod[cvMethodNo], "----------\n")

    cvCombLstNow <- cvCombLst[[cvMethodNo]][[cvNo]]

    PerformM3 <- function(parNo) {

      ### parmeter set
      genCor <- parSetDf$genCor[parNo]
      h2 <- parSetDf$h2[parNo]
      gamma <- parSetDf$gamma[parNo]
      repNo <- parSetDf$rep[parNo]

      ### prepare saving directory
      ### model parameter by BGLR
      dirM3Par <- paste0(dirM3, cvMethod[cvMethodNo], "/CV", cvNo, "/modelPar/")
      if (!dir.exists(dirM3Par)) {
        dir.create(dirM3Par, recursive = TRUE)
      }
      ### all result by BGLR
      dirM3Fit <- paste0(dirM3, cvMethod[cvMethodNo], "/CV", cvNo, "/allRes/")
      if (!dir.exists(dirM3Fit)) {
        dir.create(dirM3Fit, recursive = TRUE)
      }

      ### load phenotype data
      fileNamePhenoINow <- paste0(dirPheno, "0.3.scaled_phenotype_genCor_", genCor,
                                  "_h2_", h2, "_gamma_", gamma, "_nRep_", nRep, ".rds")
      phenoLst <- readRDS(fileNamePhenoINow)

      ### fitting for each crop
      lapply(X = 1:2, FUN = function(cropNo) {
        ### intercropping phenotype for each crop
        phenoVecNow <- phenoLst[[repNo]][[cropNo]]
        ### for save
        fileNameParameter <- paste0(dirM3Par, scriptID, cropName[cropNo],
                                    "_parNo_", parNo, "_", cvMethod[cvMethodNo],
                                    "_CV_", cvNo, "_")
        fileNameFitting <- paste0(dirM3Fit, scriptID, cropName[cropNo],
                                  "_fitting_allRes_parNo_", parNo, "_",
                                  cvMethod[cvMethodNo], "_CV_", cvNo, ".rds")
        if (!file.exists(fileNameFitting)) {
          ### training data for monocrop
          trainPos <- match(paste(cvCombLstNow$Sb, cvCombLstNow$Gm),
                            paste(pairNameDf$Sb, pairNameDf$Gm))
          ### extract selected individuals only
          yNow <- phenoVecNow[trainPos]
          ### for model
          ETA <- list(GenoSb = list(K = K1[trainPos, trainPos], model = "RKHS"),
                      GenoGm = list(K = K2[trainPos, trainPos], model = "RKHS"),
                      GenoGbyG = list(K = K3[trainPos, trainPos], model = "RKHS"))
          ### fitting
          fmResM3 <- BGLR(y = yNow,
                          saveAt = fileNameParameter, ETA = ETA,
                          nIter = nIter, burnIn = burnIn, thin = thin,
                          verbose = FALSE)
          ### save all fitting results of BayesB
          saveRDS(object = fmResM3, file = fileNameFitting)

          ### predict new individuals
          mu <- fmResM3$mu
          sigma2e <- fmResM3$varE
          sigmaK1 <- fmResM3$ETA[[1]]$varU
          sigmaK2 <- fmResM3$ETA[[2]]$varU
          sigmaK3 <- fmResM3$ETA[[3]]$varU
          kMixoo <- sigmaK1 * K1[trainPos, trainPos] +
            sigmaK2 * K2[trainPos, trainPos] +
            sigmaK3 * K3[trainPos, trainPos] ### observed inds' K mat
          kMixno <- sigmaK1 * K1[-trainPos, trainPos] +
            sigmaK2 * K2[-trainPos, trainPos] +
            sigmaK3 * K3[-trainPos, trainPos] ### new inds' K mat
          cMat <- kMixoo + sigma2e * diag(length(trainPos))
          yNew <- mu + kMixno %*% solve(cMat, (yNow - mu))
          ### save
          fileNamePred <- paste0(dirM3Fit, scriptID, cropName[cropNo],
                                 "_predicted_yNew_parNo_", parNo, "_",
                                 cvMethod[cvMethodNo], "_CV_", cvNo, ".rds")
          saveRDS(object = yNew, file = fileNamePred)
        }
      })
    }
    pbmcapply::pbmclapply(X = 1:nPar, FUN = PerformM3, mc.cores = nCores)
    gc(reset = TRUE);gc(reset = TRUE)
  }
}




###### 5. Model4 (DGE + IGE + Multitrait) ######
##### 5.1. Prediction #####
for (cvMethodNo in 1:3) {
  for (cvNo in 1:nCV) {
    cat("---------- Model: 4 / CV:", cvNo," / cvMethod:", cvMethod[cvMethodNo], "----------\n")

    cvCombLstNow <- cvCombLst[[cvMethodNo]][[cvNo]]

    PerformM4 <- function(parNo) {

      ### parmeter set
      genCor <- parSetDf$genCor[parNo]
      h2 <- parSetDf$h2[parNo]
      gamma <- parSetDf$gamma[parNo]
      repNo <- parSetDf$rep[parNo]

      ### prepare saving directory
      ### model parameter by BGLR
      dirM4Par <- paste0(dirM4, cvMethod[cvMethodNo], "/CV", cvNo, "/modelPar/")
      if (!dir.exists(dirM4Par)) {
        dir.create(dirM4Par, recursive = TRUE)
      }
      ### all result by BGLR
      dirM4Fit <- paste0(dirM4, cvMethod[cvMethodNo], "/CV", cvNo, "/allRes/")
      if (!dir.exists(dirM4Fit)) {
        dir.create(dirM4Fit, recursive = TRUE)
      }

      ### load phenotype data
      fileNamePhenoINow <- paste0(dirPheno, "0.3.scaled_phenotype_genCor_", genCor,
                                  "_h2_", h2, "_gamma_", gamma, "_nRep_", nRep, ".rds")
      phenoLst <- readRDS(fileNamePhenoINow)

      ### fitting
      ### intercropping phenotype for two crops
      phenoMatNow <- cbind(phenoLst[[repNo]]$Sb, phenoLst[[repNo]]$Gm)
      ### for save
      fileNameParameter <- paste0(dirM4Par, scriptID, "Multitrait_parNo_", parNo,
                                  "_", cvMethod[cvMethodNo], "_CV_", cvNo, "_")
      fileNameFitting <- paste0(dirM4Fit, scriptID,
                                "Multitrait_fitting_allRes_parNo_", parNo, "_",
                                cvMethod[cvMethodNo], "_CV_", cvNo, ".rds")
      if (!file.exists(fileNameFitting)) {
        ### extract selected individuals only
        yNow <- matrix(NA, nrow = nrow(phenoMatNow), ncol = 2)
        rownames(yNow) <- rownames(phenoMatNow)
        ### training data for monocrop
        trainPos <- match(paste(cvCombLstNow$Sb, cvCombLstNow$Gm),
                          paste(pairNameDf$Sb, pairNameDf$Gm))
        yNow[trainPos, ] <- phenoMatNow[trainPos, ]
        ### for model
        ETA <- list(GenoSb = list(K = K1, model = "RKHS", Cov = list(type = "UN")),
                    GenoGm = list(K = K2, model = "RKHS", Cov = list(type = "UN")))
        resCov <- list(S0 = diag(2), type = "UN")
        ### fitting
        fmResM4 <- Multitrait(y = yNow, resCov = resCov,
                              saveAt = fileNameParameter, ETA = ETA,
                              nIter = nIter, burnIn = burnIn, thin = thin,
                              verbose = FALSE)
        ### save all fitting results of BayesB
        saveRDS(object = fmResM4, file = fileNameFitting)
      }
    }
    pbmcapply::pbmclapply(X = 1:nPar, FUN = PerformM4, mc.cores = nCores)
    gc(reset = TRUE);gc(reset = TRUE)
  }
}




###### 6. Model5 (DGE + IGE + SMA + Multitrait) ######
##### 6.1. Prediction #####
for (cvMethodNo in 1:3) {
  for (cvNo in 1:nCV) {
    cat("---------- Model: 5 / CV:", cvNo," / cvMethod:", cvMethod[cvMethodNo], "----------\n")
    
    cvCombLstNow <- cvCombLst[[cvMethodNo]][[cvNo]]
    
    PerformM5 <- function(parNo) {
      
      ### parmeter set
      genCor <- parSetDf$genCor[parNo]
      h2 <- parSetDf$h2[parNo]
      gamma <- parSetDf$gamma[parNo]
      repNo <- parSetDf$rep[parNo]
      
      ### prepare saving directory
      ### model parameter by BGLR
      dirM5Par <- paste0(dirM5, cvMethod[cvMethodNo], "/CV", cvNo, "/modelPar/")
      if (!dir.exists(dirM5Par)) {
        dir.create(dirM5Par, recursive = TRUE)
      }
      ### all result by BGLR
      dirM5Fit <- paste0(dirM5, cvMethod[cvMethodNo], "/CV", cvNo, "/allRes/")
      if (!dir.exists(dirM5Fit)) {
        dir.create(dirM5Fit, recursive = TRUE)
      }
      
      ### load phenotype data
      fileNamePhenoINow <- paste0(dirPheno, "0.3.scaled_phenotype_genCor_", genCor,
                                  "_h2_", h2, "_gamma_", gamma, "_nRep_", nRep, ".rds")
      phenoLst <- readRDS(fileNamePhenoINow)
      
      ### fitting
      ### intercropping phenotype for two crops
      phenoMatNow <- cbind(phenoLst[[repNo]]$Sb, phenoLst[[repNo]]$Gm)
      ### for save
      fileNameParameter <- paste0(dirM5Par, scriptID, "Multitrait_parNo_", parNo, 
                                  "_", cvMethod[cvMethodNo], "_CV_", cvNo, "_")
      fileNameFitting <- paste0(dirM5Fit, scriptID, 
                                "Multitrait_fitting_allRes_parNo_", parNo, "_", 
                                cvMethod[cvMethodNo], "_CV_", cvNo, ".rds")
      if (!file.exists(fileNameFitting)) {
        ### training data for monocrop
        trainPos <- match(paste(cvCombLstNow$Sb, cvCombLstNow$Gm),
                          paste(pairNameDf$Sb, pairNameDf$Gm))
        ### extract selected individuals only
        yNow <- phenoMatNow[trainPos, ]
        ### for model
        ETA <- list(GenoSb = list(K = K1[trainPos, trainPos], model = "RKHS", 
                                  Cov = list(type = "UN")), 
                    GenoGm = list(K = K2[trainPos, trainPos], model = "RKHS", 
                                  Cov = list(type = "UN")), 
                    GenoGbyG = list(K = K3[trainPos, trainPos], model = "RKHS", 
                                    Cov = list(type = "DIAG")))
        resCov <- list(S0 = diag(2), type = "UN")
        ### fitting
        fmResM5 <- Multitrait(y = yNow, resCov = resCov, 
                              saveAt = fileNameParameter, ETA = ETA,
                              nIter = nIter, burnIn = burnIn, thin = thin,
                              verbose = FALSE)
        ### save all fitting results of BayesB
        saveRDS(object = fmResM5, file = fileNameFitting)
        
        ### predict new individuals
        muVec <- fmResM5$mu
        u1o <- fmResM5$ETA$GenoSb$u ### trained u1
        u2o <- fmResM5$ETA$GenoGm$u ### trained u2
        u3o <- fmResM5$ETA$GenoGbyG$u ### trained u3
        ### prediction of new u1
        k1oo <- K1[trainPos, trainPos] ### K matrix for observed lines
        k1no <- K1[-trainPos, trainPos] ### K matrix for new & observed lines
        diag(k1oo) <- diag(k1oo) + 1e-08
        u1new <- k1no %*% solve(k1oo) %*% u1o
        ### prediction of new u2
        k2oo <- K2[trainPos, trainPos] ### K matrix for observed lines
        k2no <- K2[-trainPos, trainPos] ### K matrix for new & observed lines
        diag(k2oo) <- diag(k2oo) + 1e-08
        u2new <- k2no %*% solve(k2oo) %*% u2o
        ### prediction of new u3
        k3oo <- K3[trainPos, trainPos] ### K matrix for observed lines
        k3no <- K3[-trainPos, trainPos] ### K matrix for new & observed lines
        diag(k3oo) <- diag(k3oo) + 1e-08
        u3new <- k3no %*% solve(k3oo) %*% u3o
        ### prediction of phenotypic value for new individuals
        yNew <- matrix(rep(muVec, each = nrow(u1new)), ncol = 2) + 
          u1new + u2new + u3new
        colnames(yNew) <- cropName
        ### save
        fileNamePred <- paste0(dirM5Fit, scriptID, 
                               "Multitrait_predicted_yNew_parNo_", parNo, "_", 
                               cvMethod[cvMethodNo], "_CV_", cvNo, ".rds")
        saveRDS(object = yNew, file = fileNamePred)
      }
    }
    pbmcapply::pbmclapply(X = 1:nPar, FUN = PerformM5, mc.cores = nCores)
    gc(reset = TRUE);gc(reset = TRUE)
  }
}