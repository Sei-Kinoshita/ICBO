###### 1.0. Prepare genomic kernels for prediction and optimization ######
### The historical filename is retained; full-panel diagnostic fits are excluded.
rm(list = ls())
gc(reset = TRUE)
dirGeno <- "data/genotype/"
dirPheno <- "data/phenotype/"
dirScr <- "midstream/1.0.GP_Fitting/"
dir.create(dirScr, recursive = TRUE, showWarnings = FALSE)

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
parNoVec <- seq_len(nPar)

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
fileNameK1 <- paste0(dirScr, "1.0.Sb_KMat.rds")
fileNameK2 <- paste0(dirScr, "1.0.Gm_KMat.rds")
fileNameK3 <- paste0(dirScr, "1.0.GbyG_KMat.rds")
if (all(file.exists(c(fileNameK1, fileNameK2, fileNameK3)))) {
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

cat("Prepared the three kernels used by scripts 1.1 and 2.0.\n")
