###### 0. Settings ######
##### 0.1. Remove objects #####
rm(list = ls())
gc(reset = TRUE);gc(reset = TRUE)

##### 0.2. Load Packages and Functions #####
require(rrBLUP)
# require(mvtnorm)

##### 0.3. Set scenario parameter #####
### total number of QTNs
nQTNSb <- 3000
nQTNGm <- 2000
nQtnLst <- list(Sb = nQTNSb, Gm = nQTNGm)
### number of QTNs for each chr
nQTNEachChrSb <- 300
nQTNEachChrGm <- 100
nQtnEachChrLst <- list(Sb = nQTNEachChrSb, Gm = nQTNEachChrGm)
### genetic correlation of two traits (DGE and IGE)
genCorVec <- c(-0.7, -0.1, 0.7)
### heritability (DGE+ IGE, only GMA)
h2Vec <- c(0.3, 0.7)
### SMA / GMA
gammaVec <- c(0, 0.2, 0.4)

##### 0.4. Define Parameters #####
scriptID <- "0.3."

cropName <- c("Sb", "Gm")
nCrop <- 2
### number of repetition for one scenario
nRep <- 10
### number of seed for sampling phenotype
nSeed <- length(genCorVec) * length(h2Vec) * length(gammaVec) * nRep

##### 0.5. Define Directories #####
dataDir <- "data/"
genoDir <- paste0(dataDir, "genotype/")
mrkEffDir <- paste0(dataDir, "trueMrkEff/")

saveDir <- paste0(dataDir, "phenotype/")
if (!dir.exists(saveDir)) {
  dir.create(saveDir, recursive = TRUE)
}




###### 1. Load Data ######
### Genotype Data
fileNameGenoSb <- paste0(genoDir, "0.1.Sorghum79_geno_Maf_0.1_LD_0.6.rds")
genoSbMat <- readRDS(file = fileNameGenoSb) - 1
fileNameGenoGm <- paste0(genoDir, "0.1.Soybean66_geno_Maf_0.1_LD_0.6.rds")
genoGmMat <- readRDS(file = fileNameGenoGm) - 1

genoLst <- list(Sb = genoSbMat, Gm = genoGmMat)
nIndVec <- c(nrow(genoSbMat), nrow(genoGmMat))

### Chromosome Information
fileNameGenoInfoSb <- paste0(genoDir, "0.1.Sorghum79_markerInfo_Maf_0.1_LD_0.6.rds")
mrkInfoSb <- readRDS(file = fileNameGenoInfoSb)
fileNameGenoInfoGm <- paste0(genoDir, "0.1.Soybean66_markerInfo_Maf_0.1_LD_0.6.rds")
mrkInfoGm <- readRDS(file = fileNameGenoInfoGm)

mrkInfoLst <- list(Sb = mrkInfoSb, Gm = mrkInfoGm)

nChrVec <- c(length(unique(mrkInfoSb$Chr)), length(unique(mrkInfoGm$Chr)))




###### 2. Modify Data ######
##### 2.1. Calculate GRM (-1/0/1) #####
### sorghum
fileNameGrmSb <- paste0(genoDir, scriptID, "Sb_GRM_AMat_nInd_", nIndVec[1], ".rds")
if (!file.exists(fileNameGrmSb)) {
  grmSb <- A.mat(X = genoSbMat)
  saveRDS(object = grmSb, file = fileNameGrmSb)
} else {
  grmSb <- readRDS(file = fileNameGrmSb)
}
### soybean
fileNameGrmGm <- paste0(genoDir, scriptID, "Gm_GRM_AMat_nInd_", nIndVec[2], ".rds")
if (!file.exists(fileNameGrmGm)) {
  grmGm <- A.mat(X = genoGmMat)
  saveRDS(object = grmGm, file = fileNameGrmGm)
} else {
  grmGm <- readRDS(file = fileNameGrmGm)
}
### Kronecker product
fileNameGrmKron <- paste0(genoDir, scriptID, "Kronecker_GRM_AMat_twoCrop_nDim_", 
                          nIndVec[1] * nIndVec[2], ".rds")
if (!file.exists(fileNameGrmKron)) {
  grmKron <- kronecker(X = grmSb, Y = grmGm) ### sorghum x soybean
  saveRDS(object = grmKron, file = fileNameGrmKron)
} else {
  grmKron <- readRDS(file = fileNameGrmKron)
}


##### 2.2. for sampling raw SMA #####
fileNameEigKron <- paste0(genoDir, scriptID, "Eigen_Kronecker_GRM_AMat_twoCrop.rds")
if (!file.exists(fileNameEigKron)) {
  eig <- eigen(grmKron, symmetric = TRUE)
  saveRDS(object = eig, file = fileNameEigKron)
} else {
  eig <- readRDS(file = fileNameEigKron)
}
lambda <- eig$values
Q <- eig$vectors
lambda[lambda < 0] <- 0


##### 2.3. Fix intercrop pair #####
fileNamePairName <- paste0(saveDir, scriptID, "pairNameDf_5214Pair.rds")
if (!file.exists(fileNamePairName)) {
  ### 5214 pair names (sb1 x gm1 ~ gm66, sb2 x gm1 ~ gm66, ..., same as grmKron)
  pairNameDf <- data.frame(Sb = rep(rownames(grmSb), each = nIndVec[2]), 
                           Gm = rep(rownames(grmGm), nIndVec[1]))
  rownames(grmKron) <- colnames(grmKron) <- 
    apply(X = pairNameDf, MARGIN = 1, FUN = paste0, collapse = "_")
  ### save
  saveRDS(object = pairNameDf, file = fileNamePairName)
} else {
  pairNameDf <- readRDS(file = fileNamePairName)
}
nPair <- nrow(pairNameDf)


##### 2.4. Parameter Set #####
fileNameParSet <- paste0(saveDir, scriptID, "parmeterSet.rds")
if (!file.exists(fileNameParSet)) {
  parDf <- expand.grid(list(genCorVec, h2Vec, gammaVec, 1:nRep))
  colnames(parDf) <- c("genCor", "h2", "gamma", "rep")
  saveRDS(object = parDf, file = fileNameParSet)
} else {
  parDf <- readRDS(file = fileNameParSet)
}
nPar <- nrow(parDf)


##### 2.5. Set seed for sampling #####
fileNameSeed <- paste0(saveDir, scriptID, "seed.rds")
if (!file.exists(fileNameSeed)) {
  ### set seed for each variable
  seedLst <- lapply(X = 1:2, FUN = function(x) {
    ### for each crop
    seedNow <- lapply(X = 1:nCrop, FUN = function(cropNo) {
      seed <- sample(x = 1e6, size = nSeed, replace = FALSE)
      ### for each scenario
      seedDf <- data.frame(seed = seed, parDf)
      return(seedDf)
    })
    names(seedNow) <- c("Sb", "Gm")
    return(seedNow)
  })
  names(seedLst) <- c("u", "e")
  saveRDS(object = seedLst, file = fileNameSeed)
} else {
  seedLst <- readRDS(file = fileNameSeed)
}




###### 3. Generate Phenotype ######
for (genCor in genCorVec) {
  for (h2 in h2Vec) {
    for (gamma in gammaVec) {
      
      fileNamePhenoRes <- paste0(saveDir, scriptID, "phenotype_allRes_genCor_", 
                                 genCor, "_h2_", h2, "_gamma_", gamma, 
                                 "_nRep_", nRep, ".rds")
      if (!file.exists(fileNamePhenoRes)) {
        ### result for one scenario
        resLst <- lapply(X = 1:nRep, FUN = function(repNo) {
          
          # repNo <- 1
          
          ### prepare for simulating phenotype
          phenoNo <- which(genCorVec == genCor) * which(h2Vec == h2) * repNo
          fileNameMrkEffLstSb <- paste0(mrkEffDir, "0.2.Sb_trueMrkEffLst_nQTN_3000_genCor_", 
                                        genCor, "_nRep_60.rds")
          mrkEffMatSb <- readRDS(file = fileNameMrkEffLstSb)[[phenoNo]]
          fileNameMrkEffLstGm <- paste0(mrkEffDir, "0.2.Gm_trueMrkEffLst_nQTN_2000_genCor_", 
                                        genCor, "_nRep_60.rds")
          mrkEffMatGm <- readRDS(file = fileNameMrkEffLstGm)[[phenoNo]]
          
          seedIdx <- which(seedLst$u$Sb$genCor == genCor & 
                             seedLst$u$Sb$h2 == h2 & 
                             seedLst$u$Sb$gamma == gamma & 
                             seedLst$u$Sb$rep == repNo)
          
          ### calculate dA, dB (DGE)
          gmaSbMat <- tcrossprod(genoSbMat, t(mrkEffMatSb))
          gmaGmMat <- tcrossprod(genoGmMat, t(mrkEffMatGm))
          dA <- gmaSbMat[pairNameDf$Sb, "DGE"] ### A: sorghum
          dB <- gmaGmMat[pairNameDf$Gm, "DGE"] ### B: soybean
          
          ### calculate sAB, sBA (IGE)
          sAtoB <- gmaSbMat[pairNameDf$Sb, "IGE"] ### IGE: sorghum to soybean
          sBtoA <- gmaGmMat[pairNameDf$Gm, "IGE"] ### IGE: soybean to sorghum
          
          ### scaling IGE to be DGE:IGE = 1:1
          sAtoB <- sAtoB * sqrt(var(dB) / var(sAtoB))
          sBtoA <- sBtoA * sqrt(var(dA) / var(sBtoA))
          
          ### sample uA, uB (raw SMA)
          set.seed(seed = seedLst$u$Sb$seed[seedIdx])
          # uA <- rmvnorm(n = 1, mean = rep(0, nPair), sigma = grmKron)
          # uA <- rep(0, nPair) + 
          #   as.vector(t(grmKronCholesky) %*% rnorm(nPair, mean = 0, sd = 1))
          uA <- as.vector(rep(0, nPair) + 
                            Q %*% (sqrt(lambda) * rnorm(nPair, mean = 0, sd = 1)))
          set.seed(seed = seedLst$u$Gm$seed[seedIdx])
          uB <- as.vector(rep(0, nPair) + 
                            Q %*% (sqrt(lambda) * rnorm(nPair, mean = 0, sd = 1)))
          
          ### scaled SMA (=vA, vB) based on gammaA, B
          ### cA, cB: scaling coefficient
          cA <- sqrt((gamma * var(dA + sBtoA)) / var(uA))
          cB <- sqrt((gamma * var(dB + sAtoB)) / var(uB))
          vA <- cA * uA
          vB <- cB * uB
          
          ### sample eA, eB (error term)
          eAVar <- (var(dA + sBtoA) / h2) - var(dA + sBtoA) - var(vA)
          eBVar <- (var(dB + sAtoB) / h2) - var(dB + sAtoB) - var(vB)
          set.seed(seed = seedLst$e$Sb$seed[seedIdx])
          eA <- rnorm(n = nPair, mean = 0, sd = sqrt(eAVar))
          set.seed(seed = seedLst$e$Gm$seed[seedIdx])
          eB <- rnorm(n = nPair, mean = 0, sd = sqrt(eBVar))
          
          ### generate yA, yB (phenotype)
          yA <- dA + sBtoA + vA + eA
          yB <- dB + sAtoB + vB + eB
          
          # ### check heritability
          # var(dA + sBtoA) / var(yA)
          # var(dB + sAtoB) / var(yB)
          
          ### scaled yA, yB
          scaled_yA <- scale(yA)
          scaled_yB <- scale(yB)
          
          ### save results
          res <- list(dA = dA, dB = dB, ### DGE
                      sAtoB = sAtoB, sBtoA = sBtoA, ### IGE
                      uA = uA, uB = uB, ### raw SMA
                      vA = vA, vB = vB, ### scaled SMA
                      eA = eA, eB = eB, ### error
                      yA = yA, yB = yB, ### phenotype
                      scaled_yA = scaled_yA, scaled_yB = scaled_yB) ### scaled phenotype
          return(res)
        })
        names(resLst) <- paste0("Rep_", 1:nRep)
        
        ### save all sampled result
        fileNameAllPheno <- paste0(saveDir, scriptID, "allSampledRes_phenotype_genCor_", 
                                   genCor, "_h2_", h2, "_gamma_", gamma, 
                                   "_nRep_", nRep, ".rds")
        saveRDS(object = resLst, file = fileNameAllPheno)
        ### save only scaled phenotype
        phenoLst <- lapply(X = resLst, FUN = function(x) {
          list(Sb = x$scaled_yA, Gm = x$scaled_yB)
        })
        names(phenoLst) <- paste0("Rep_", 1:nRep)
        fileNameOnlyPheno <- paste0(saveDir, scriptID, "scaled_phenotype_genCor_", 
                                    genCor, "_h2_", h2, "_gamma_", gamma, 
                                    "_nRep_", nRep, ".rds")
        saveRDS(object = phenoLst, file = fileNameOnlyPheno)
      } 
    }
  }
}