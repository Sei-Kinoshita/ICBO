###### 0. Settings ######
##### 0.1. Remove objects #####
rm(list = ls())
gc(reset = TRUE);gc(reset = TRUE)

##### 0.2. Load packages #####
require(BGLR)

##### 0.3. Set parameters #####
scriptID <- "1.0."

cropName <- c("Sb", "Gm")
### total number of QTNs
nQtnLst <- list(Sb = 3000, Gm = 2000)
### number of QTNs for each chr
nQtnEachChrLst <- list(Sb = 300, Gm = 100)
### genetic correlation of two traits (DGE and IGE)
genCorVec <- c(-0.7, -0.1, 0.7)
### heritability (DGE+ IGE, only GMA)
h2Vec <- c(0.3, 0.7)
### SMA / GMA
gammaVec <- c(0, 0.2, 0.4)

nRep <- 10 ### repetition of phenotype
nCores <- 30

### for BGLR
nIter <- 10000
burnIn <- 2000
thin <- 20
# nIter <- 10
# burnIn <- 2
# thin <- 1


##### 0.4. Define directories #####
dirGeno <- "data/genotype/"
dirPheno <- "data/phenotype/"

dirScr <- paste0("midstream/", scriptID, "GP_Fitting/") 
if (!dir.exists(dirScr)) {
  dir.create(dirScr, recursive = TRUE)
}

#### 0.4.1. RKHS directories ####
dirM1 <- paste0(dirScr, "Model1_DGE/")
if (!dir.exists(dirM1)) {
  dir.create(dirM1, recursive = TRUE)
}

#### 0.4.2. BLUP directories ####
dirM2 <- paste0(dirScr, "Model2_DGE_IGE/")
if (!dir.exists(dirM2)) {
  dir.create(dirM2, recursive = TRUE)
}

#### 0.4.3. BLUP GbyG directories ####
dirM3 <- paste0(dirScr, "Model3_GMA_SMA/")
if (!dir.exists(dirM3)) {
  dir.create(dirM3, recursive = TRUE)
}

#### 0.4.4. Model 2 + Multitrait directories ####
dirM4 <- paste0(dirScr, "Model4_Multi_DGE_IGE/")
if (!dir.exists(dirM4)) {
  dir.create(dirM4, recursive = TRUE)
}

#### 0.4.5. Model 3 + Multitrait directories ####
dirM5 <- paste0(dirScr, "Model5_Multi_GMA_SMA/")
if (!dir.exists(dirM5)) {
  dir.create(dirM5, recursive = TRUE)
}




###### 1. Load Data ######
# ##### 1.1. Genotype Data #####
# fileNameGenoSb <- paste0(dirGeno, "0.1.Sorghum79_geno_Maf_0.1_LD_0.6.rds")
# genoSbMat <- readRDS(file = fileNameGenoSb) - 1
# fileNameGenoGm <- paste0(dirGeno, "0.1.Soybean66_geno_Maf_0.1_LD_0.6.rds")
# genoGmMat <- readRDS(file = fileNameGenoGm) - 1
# 
# genoLst <- list(Sb = genoSbMat, Gm = genoGmMat)
# nIndVec <- c(nrow(genoSbMat), nrow(genoGmMat))

##### 1.2. GRM #####
fileNameGrmSb <- paste0(dirGeno, "0.3.Sb_GRM_AMat_nInd_79.rds")
grmSb <- readRDS(file = fileNameGrmSb)
fileNameGrmGm <- paste0(dirGeno, "0.3.Gm_GRM_AMat_nInd_66.rds")
grmGm <- readRDS(file = fileNameGrmGm)
fileNameGrmGbyG <- paste0(dirGeno, "0.3.Kronecker_GRM_AMat_twoCrop_nDim_5214.rds")
grmGbyG <- readRDS(file = fileNameGrmGbyG)

##### 1.3. Parameter Set #####
fileNameParSet <- paste0(dirPheno, "0.3.parmeterSet.rds")
parSetDf <- readRDS(file = fileNameParSet)
nPar <- nrow(parSetDf)

##### 1.5. Pair Name (fixed order) #####
fileNamePairName <- paste0(dirPheno, "0.3.pairNameDf_5214Pair.rds")
pairNameDf <- readRDS(file = fileNamePairName)

##### 1.6. design matrix for model 2&3 #####
pairNameDf$Sb <- factor(pairNameDf$Sb, levels = rownames(grmSb))
pairNameDf$Gm <- factor(pairNameDf$Gm, levels = rownames(grmGm))
zMatSb <- model.matrix(~Sb - 1, pairNameDf)
zMatGm <- model.matrix(~Gm - 1, pairNameDf)
zMatGbyG <- diag(nrow(grmGbyG))

### K matrix for model 2&3
fileNameK1 <- paste0(dirScr, scriptID, "Sb_KMat.rds")
fileNameK2 <- paste0(dirScr, scriptID, "Gm_KMat.rds")
fileNameK3 <- paste0(dirScr, scriptID, "GbyG_KMat.rds")
if (all(file.exists(file1 = fileNameK1, file2 = fileNameK2, file3 = fileNameK3))) {
  K1 <- readRDS(file = fileNameK1)
  K2 <- readRDS(file = fileNameK2)
  K3 <- readRDS(file = fileNameK3)
} else {
  K1 <- zMatSb %*% grmSb %*% t(zMatSb)
  K2 <- zMatGm %*% grmGm %*% t(zMatGm)
  K3 <- grmGbyG ### because zMatGbyG is identity matrix
  saveRDS(object = K1, file = fileNameK1)
  saveRDS(object = K2, file = fileNameK2)
  saveRDS(object = K3, file = fileNameK3)
}



rm(grmSb, grmGm, grmGbyG, zMatSb, zMatGm, zMatGbyG)
gc(reset = TRUE);gc(reset = TRUE)

###### 2. GP by RKHS (Model1 only DGE) ######
##### 2.1. Fitting #####
cat("---------- Model: 1 Fitting ----------\n")
PerformM1 <- function(parNo) {

  ### parmeter set
  genCor <- parSetDf$genCor[parNo]
  h2 <- parSetDf$h2[parNo]
  gamma <- parSetDf$gamma[parNo]
  repNo <- parSetDf$rep[parNo]

  ### prepare saving directory
  dirM1Par <- paste0(dirM1, "modelPar/") ### model parameter by BGLR
  if (!dir.exists(dirM1Par)) {
    dir.create(dirM1Par, recursive = TRUE)
  }
  dirM1Fit <- paste0(dirM1, "allRes/") ### all fitting result by BGLR
  if (!dir.exists(dirM1Fit)) {
    dir.create(dirM1Fit, recursive = TRUE)
  }

  ### load phenotype data
  fileNamePhenoINow <- paste0(dirPheno, "0.3.scaled_phenotype_genCor_", genCor,
                              "_h2_", h2, "_gamma_", gamma, "_nRep_", nRep, ".rds")
  phenoLst <- readRDS(fileNamePhenoINow)

  lapply(X = 1:2, FUN = function(cropNo) {

    phenoLstNow <- phenoLst[[repNo]]
    if (cropNo == 1) {
      yNow <- phenoLstNow$Sb
    } else {
      yNow <- phenoLstNow$Gm
    }
    fileNameParameter <- paste0(dirM1Par, scriptID, cropName[cropNo], "_",
                                parNo, "_")
    fileNameFitting <- paste0(dirM1Fit, scriptID, cropName[cropNo],
                              "_parNo_", parNo, "_fitting_allRes.rds")
    if (!file.exists(fileNameFitting)) {
      ### for model
      if (cropNo == 1) {
        ETA <- list(GenoRKHS = list(K = K1, model = "RKHS"))
      } else {
        ETA <- list(GenoRKHS = list(K = K2, model = "RKHS"))
      }
      ### fitting
      fitResM1 <- BGLR(y = yNow, saveAt = fileNameParameter, ETA = ETA,
                       nIter = nIter, burnIn = burnIn, thin = thin,
                       verbose = FALSE)
      saveRDS(object = fitResM1, file = fileNameFitting)
      gc(reset = TRUE);gc(reset = TRUE)
    }
  })
}
pbmcapply::pbmclapply(X = 1:nPar, FUN = PerformM1, mc.cores = nCores)
gc(reset = TRUE);gc(reset = TRUE)




###### 3. GP (Model2 DGE + IGE) ######
##### 3.1. Fitting #####
cat("---------- Model: 2 Fitting ----------\n")
PerformM2 <- function(parNo) {

  ### parmeter set
  genCor <- parSetDf$genCor[parNo]
  h2 <- parSetDf$h2[parNo]
  gamma <- parSetDf$gamma[parNo]
  repNo <- parSetDf$rep[parNo]

  ### prepare saving directory
  dirM2Par <- paste0(dirM2, "modelPar/") ### model parameter by BGLR
  if (!dir.exists(dirM2Par)) {
    dir.create(dirM2Par, recursive = TRUE)
  }
  dirM2Fit <- paste0(dirM2, "allRes/") ### all fitting result by BGLR
  if (!dir.exists(dirM2Fit)) {
    dir.create(dirM2Fit, recursive = TRUE)
  }

  ### load phenotype data
  fileNamePhenoINow <- paste0(dirPheno, "0.3.scaled_phenotype_genCor_", genCor,
                              "_h2_", h2, "_gamma_", gamma, "_nRep_", nRep, ".rds")
  phenoLst <- readRDS(fileNamePhenoINow)

  lapply(X = 1:2, FUN = function(cropNo) {

    phenoLstNow <- phenoLst[[repNo]]
    if (cropNo == 1) {
      yNow <- phenoLstNow$Sb
    } else {
      yNow <- phenoLstNow$Gm
    }
    fileNameParameter <- paste0(dirM2Par, scriptID, cropName[cropNo], "_",
                                parNo, "_")
    fileNameFitting <- paste0(dirM2Fit, scriptID, cropName[cropNo],
                              "_parNo_", parNo, "_fitting_allRes.rds")
    if (!file.exists(fileNameFitting)) {
      ### for model
      ETA <- list(GenoSb = list(K = K1, model = "RKHS"),
                  GenoGm = list(K = K2, model = "RKHS"))
      ### fitting
      fitResM2 <- BGLR(y = yNow, saveAt = fileNameParameter, ETA = ETA,
                       nIter = nIter, burnIn = burnIn, thin = thin,
                       verbose = FALSE)
      saveRDS(object = fitResM2, file = fileNameFitting)
    }
  })
  gc(reset = TRUE);gc(reset = TRUE)
}
pbmcapply::pbmclapply(X = 1:nPar, FUN = PerformM2, mc.cores = nCores)
gc(reset = TRUE);gc(reset = TRUE)




###### 4. GP (Model3 DGE + IGE + SMA) ######
##### 4.1. Fitting #####
cat("---------- Model: 3 Fitting ----------\n")
PerformM3 <- function(parNo) {

  ### parmeter set
  genCor <- parSetDf$genCor[parNo]
  h2 <- parSetDf$h2[parNo]
  gamma <- parSetDf$gamma[parNo]
  repNo <- parSetDf$rep[parNo]

  ### prepare saving directory
  dirM3Par <- paste0(dirM3, "modelPar/") ### model parameter by BGLR
  if (!dir.exists(dirM3Par)) {
    dir.create(dirM3Par, recursive = TRUE)
  }
  dirM3Fit <- paste0(dirM3, "allRes/") ### all fitting result by BGLR
  if (!dir.exists(dirM3Fit)) {
    dir.create(dirM3Fit, recursive = TRUE)
  }

  ### load phenotype data
  fileNamePhenoINow <- paste0(dirPheno, "0.3.scaled_phenotype_genCor_", genCor,
                              "_h2_", h2, "_gamma_", gamma, "_nRep_", nRep, ".rds")
  phenoLst <- readRDS(fileNamePhenoINow)

  lapply(X = 1:2, FUN = function(cropNo) {

    phenoLstNow <- phenoLst[[repNo]]
    if (cropNo == 1) {
      yNow <- phenoLstNow$Sb
    } else {
      yNow <- phenoLstNow$Gm
    }
    fileNameParameter <- paste0(dirM3Par, scriptID, cropName[cropNo], "_",
                                parNo, "_")
    fileNameFitting <- paste0(dirM3Fit, scriptID, cropName[cropNo],
                              "_parNo_", parNo, "_fitting_allRes.rds")
    if (!file.exists(fileNameFitting)) {
      ### for model
      ETA <- list(GenoSb = list(K = K1, model = "RKHS"),
                  GenoGm = list(K = K2, model = "RKHS"),
                  GenoGbyG = list(K = K3, model = "RKHS"))
      ### fitting
      fitResM3 <- BGLR(y = yNow, saveAt = fileNameParameter, ETA = ETA,
                       nIter = nIter, burnIn = burnIn, thin = thin,
                       verbose = FALSE)
      saveRDS(object = fitResM3, file = fileNameFitting)
    }
  })
  gc(reset = TRUE);gc(reset = TRUE)
}
pbmcapply::pbmclapply(X = 1:nPar, FUN = PerformM3, mc.cores = nCores)
gc(reset = TRUE);gc(reset = TRUE)




###### 5. GP (Model4 DGE + IGE + Multitrait) ######
##### 5.1. Fitting #####
cat("---------- Model: 4 Fitting ----------\n")
PerformM4 <- function(parNo) {

  ### parmeter set
  genCor <- parSetDf$genCor[parNo]
  h2 <- parSetDf$h2[parNo]
  gamma <- parSetDf$gamma[parNo]
  repNo <- parSetDf$rep[parNo]

  ### prepare saving directory
  dirM4Par <- paste0(dirM4, "modelPar/") ### model parameter by BGLR
  if (!dir.exists(dirM4Par)) {
    dir.create(dirM4Par, recursive = TRUE)
  }
  dirM4Fit <- paste0(dirM4, "allRes/") ### all fitting result by BGLR
  if (!dir.exists(dirM4Fit)) {
    dir.create(dirM4Fit, recursive = TRUE)
  }

  ### load phenotype data
  fileNamePhenoINow <- paste0(dirPheno, "0.3.scaled_phenotype_genCor_", genCor,
                              "_h2_", h2, "_gamma_", gamma, "_nRep_", nRep, ".rds")
  phenoLst <- readRDS(fileNamePhenoINow)
  phenoLstNow <- phenoLst[[repNo]]

  ### Multitrait BGLR
  fileNameParameter <- paste0(dirM4Par, scriptID, "Multitrait_",
                              parNo, "_")
  fileNameFitting <- paste0(dirM4Fit, scriptID,
                            "Multitrait_parNo_", parNo, "_fitting_allRes.rds")
  if (!file.exists(fileNameFitting)) {
    yNow <- cbind(phenoLstNow$Sb, phenoLstNow$Gm)
    ### for model
    ETA <- list(GenoSb = list(K = K1, model = "RKHS", Cov = list(type = "UN")),
                GenoGm = list(K = K2, model = "RKHS", Cov = list(type = "UN")))
    resCov <- list(S0 = diag(2), type = "UN")
    ### fitting
    fitResM4 <- Multitrait(y = yNow, resCov = resCov,
                           saveAt = fileNameParameter, ETA = ETA,
                           nIter = nIter, burnIn = burnIn, thin = thin,
                           verbose = FALSE)
    saveRDS(object = fitResM4, file = fileNameFitting)
  }
  gc(reset = TRUE);gc(reset = TRUE)
}
pbmcapply::pbmclapply(X = 1:nPar, FUN = PerformM4, mc.cores = nCores)
gc(reset = TRUE);gc(reset = TRUE)




###### 6. GP (Model5 DGE + IGE + SMA + Multitrait) ######
##### 6.1. Fitting #####
cat("---------- Model: 5 Fitting ----------\n")
PerformM5 <- function(parNo) {
  
  ### parmeter set
  genCor <- parSetDf$genCor[parNo]
  h2 <- parSetDf$h2[parNo]
  gamma <- parSetDf$gamma[parNo]
  repNo <- parSetDf$rep[parNo]
  
  ### prepare saving directory
  dirM5Par <- paste0(dirM5, "modelPar/") ### model parameter by BGLR
  if (!dir.exists(dirM5Par)) {
    dir.create(dirM5Par, recursive = TRUE)
  }
  dirM5Fit <- paste0(dirM5, "allRes/") ### all fitting result by BGLR
  if (!dir.exists(dirM5Fit)) {
    dir.create(dirM5Fit, recursive = TRUE)
  }
  
  ### load phenotype data
  fileNamePhenoINow <- paste0(dirPheno, "0.3.scaled_phenotype_genCor_", genCor, 
                              "_h2_", h2, "_gamma_", gamma, "_nRep_", nRep, ".rds")
  phenoLst <- readRDS(fileNamePhenoINow)
  phenoLstNow <- phenoLst[[repNo]]
  
  ### Multitrait BGLR
  fileNameParameter <- paste0(dirM5Par, scriptID, "Multitrait_", parNo, "_")
  fileNameFitting <- paste0(dirM5Fit, scriptID, 
                            "Multitrait_parNo_", parNo, "_fitting_allRes.rds")
  if (!file.exists(fileNameFitting)) {
    yNow <- cbind(phenoLstNow$Sb, phenoLstNow$Gm)
    ### for model
    ETA <- list(GenoSb = list(K = K1, model = "RKHS", Cov = list(type = "UN")), 
                GenoGm = list(K = K2, model = "RKHS", Cov = list(type = "UN")), 
                GenoGbyG = list(K = K3, model = "RKHS", Cov = list(type = "DIAG")))
    resCov <- list(S0 = diag(2), type = "UN")
    ### fitting
    fitResM5 <- Multitrait(y = yNow, resCov = resCov, 
                           saveAt = fileNameParameter, ETA = ETA,
                           nIter = nIter, burnIn = burnIn, thin = thin,
                           verbose = FALSE)
    saveRDS(object = fitResM5, file = fileNameFitting)
  }
  gc(reset = TRUE);gc(reset = TRUE)
}
pbmcapply::pbmclapply(X = 1:nPar, FUN = PerformM5, mc.cores = nCores)
gc(reset = TRUE);gc(reset = TRUE)