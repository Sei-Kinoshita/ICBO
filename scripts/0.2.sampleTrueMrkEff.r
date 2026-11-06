###### 1. Settings ######
##### 1.1. Remove objects #####
rm(list = ls())

##### 1.2. Load packages #####
require(pbmcapply)

### The optional local archive preserves the manuscript effects; it is not distributed.
if (file.exists("data/0.2.genotypeEffects_used.rds")) {
  message("Using archived genotype effects; no new QTN realization is sampled.")
} else {
  ### Seeds below define new realizations, not the unrecorded historical QTN RNG states.
  qtnSeedBase <- 20261106L

  ##### 1.3. Set scenario parameter #####
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

  ##### 1.4. Define parameters #####
  scriptID <- "0.2."
  ### marker position for tight linkage
  tightLinkThre <- 5
  ### for generating mrkEff
  nMrkEffRep <- 60 ### nHeritability (=2) * nGamma (=3) * 10 replications
  genCorDiffThre <- 0.1 ### allow difference with ideal genetic correlation

  cropName <- c("Sb", "Gm")
  nCrop <- 2

  nCores <- 30

  ##### 1.5. Define directories #####
  dataDir <- "data/"
  genoDir <- paste0(dataDir, "genotype/")

  dirMrkEff <- paste0(dataDir, "trueMrkEff/")
  if (!dir.exists(dirMrkEff)) {
    dir.create(dirMrkEff, recursive = TRUE)
  }

  ###### 2. Load data into R ######
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

  ###### 3. Generate Marker Effects ######
  ##### 3.1. Func: sample QTN positions #####
  SampleQtnPosTightLink <- function(linkageInfoDf, nChr, nQTNEachChr) {

    qtnNameEachChrLst <- sapply(X = 1:nChr, FUN = function(chrNo) {

      nMrkNow <- sum(linkageInfoDf$Chr == chrNo)
      mrkNoNow <- 1:nMrkNow
      linkageInfoNow <- linkageInfoDf[linkageInfoDf$Chr == chrNo, ]
      mapDisNow <- linkageInfoNow$MapDis[1:nMrkNow]
      mrkNameNow <- linkageInfoNow$MrkName[1:nMrkNow]

      qtnPosT1Now <- sample(x = 1:nMrkNow, size = nQTNEachChr, replace = FALSE)

      ### sample QTN of Trait2
      qtnPosT2Now <- NULL
      for (qtnNo in 1:length(qtnPosT1Now)) {
        qtnNoNow <- qtnPosT1Now[qtnNo]
        qtnT1MapDisNow <- mapDisNow[qtnNoNow]
        ### difference of map distance between qtnT1 and other markers
        mapDisDiff <- abs(rep(qtnT1MapDisNow, length(mapDisNow)) - mapDisNow)
        qtnCandNoT2 <- which(mapDisDiff < tightLinkThre)
        qtnCandNoT2 <- qtnCandNoT2[!(qtnCandNoT2 %in% qtnPosT1Now)] ### remove QTN T1

        if (length(qtnCandNoT2) != 0) {
          qtnNoT2 <- sample(x = qtnCandNoT2, size = 1)
        } else {
          ### if there are no markers in 5cM of QTN of Trait1
          if (qtnNoNow == 1) {
            ### if QTN of Trait1 is the start of the chromosome
            qtnCandNoT2 <- mrkNoNow[-c(qtnNoNow, qtnPosT1Now)]
            qtnNoT2 <- qtnCandNoT2[1]
          } else if (qtnNoNow == nMrkNow) {
            ### if QTN of Trait1 is the end of the chromosome
            qtnCandNoT2 <- mrkNoNow[-c(qtnNoNow, qtnPosT1Now)]
            qtnNoT2 <- rev(qtnCandNoT2)[1]
          } else {
            qtnCandNoT2 <- mrkNoNow[which(!(mrkNoNow %in% qtnPosT1Now[which(!(qtnPosT1Now %in% qtnNoNow))]))]
            qtnNoInCandNoT2 <- which(qtnCandNoT2 == qtnNoNow)
            mrkNameMatchCandNoT2 <- mrkNameNow[qtnCandNoT2]
            if (qtnNoInCandNoT2 == 1) {
              qtnNoT2InCandNoT2 <- qtnNoInCandNoT2 + 1
              qtnNoT2 <- qtnCandNoT2[qtnNoT2InCandNoT2]
            } else if (qtnNoInCandNoT2 == length(qtnCandNoT2)) {
              qtnNoT2InCandNoT2 <- qtnNoInCandNoT2 - 1
              qtnNoT2 <- qtnCandNoT2[qtnNoT2InCandNoT2]
            } else {
              ### For 2 markers adjacent to QTN T1, choose the one with higher LD
              qtnCandT2Cor1 <- cor(x = genoMat[, mrkNameMatchCandNoT2[qtnNoInCandNoT2]],
                                   y = genoMat[, mrkNameMatchCandNoT2[qtnNoInCandNoT2 + 1]]) ^ 2
              qtnCandT2Cor2 <- cor(x = genoMat[, mrkNameMatchCandNoT2[qtnNoInCandNoT2]],
                                   y = genoMat[, mrkNameMatchCandNoT2[qtnNoInCandNoT2 - 1]]) ^ 2
              qtnNoT2InCandNoT2 <- c((qtnNoInCandNoT2 + 1), (qtnNoInCandNoT2 - 1))[which.max(c(qtnCandT2Cor1, qtnCandT2Cor2))]
              qtnNoT2 <- qtnCandNoT2[qtnNoT2InCandNoT2]
            }
          }
        }

        qtnPosT2Now <- c(qtnPosT2Now, qtnNoT2)
      }

      qtnNameT1Now <- mrkNameNow[qtnPosT1Now]
      qtnNameT2Now <- mrkNameNow[qtnPosT2Now]
      qtnNameNowMat <- rbind(T1 = qtnNameT1Now, T2 = qtnNameT2Now)
      return(qtnNameNowMat)

    }, simplify = FALSE)

    qtnName2TraitMat <- do.call(what = cbind, args = qtnNameEachChrLst)
    return(qtnName2TraitMat)
  }

  ##### 3.2. Generate marker effects #####
  for (genCor in genCorVec) {

    ### parameters for generic correlation
    idealGenVarMat <- matrix(data = c(1, genCor, genCor, 1), nrow = 2, byrow = TRUE)
    cholGenVarMat <- chol(idealGenVarMat) ### Cholesky decomposition

    for (cropNo in 1:nCrop) {

      ### parameters for generic correlation
      ### a value
      aValue <- (1 - nQtnLst[[cropNo]]) / (1 + nQtnLst[[cropNo]])
      ### geno information for current crop
      mrkName <- colnames(genoLst[[cropNo]])
      nMrk <- length(mrkName)
      nQTN <- nQtnLst[[cropNo]]
      nQTNEachChr <- nQtnEachChrLst[[cropNo]]
      linkageInfoDf <- mrkInfoLst[[cropNo]]
      nChr <- length(unique(linkageInfoDf$Chr))
      genoMat <- genoLst[[cropNo]]

      fileNameMrkEffNow <- paste0(dirMrkEff, scriptID, cropName[cropNo],
                                  "_trueMrkEffLst_nQTN_", nQTN, "_genCor_", genCor,
                                  "_nRep_", nMrkEffRep, ".rds")
      if (!file.exists(fileNameMrkEffNow)) {
        mrkEffTightLinkLst <- pbmclapply(X = 1:nMrkEffRep, FUN = function(repNo) {
          set.seed(qtnSeedBase + match(genCor, genCorVec) * 10000L + cropNo * 1000L + repNo)
          loopTightLinkVec <- 0
          ### generate QTN effect
          while (loopTightLinkVec == 0) {

            mrkEffTightLinkMat <- matrix(data = 0, nrow = nMrk, ncol = 2)
            dimnames(mrkEffTightLinkMat) <- list(mrkName, c("DGE", "IGE"))

            ### generate QTN effect for each trait
            qtnMrkEffTightLinkT1 <- sapply(X = 1:nQTN, FUN = function(qtnNo) {
              (aValue ^ qtnNo) * sample(x = c(-1, 1), size = 1)
            })
            qtnMrkEffTightLinkT2 <- sapply(X = 1:nQTN, FUN = function(qtnNo) {
              (aValue ^ qtnNo) * sample(x = c(-1, 1), size = 1)
            })
            qtnMrkEffTightLinkAMat <- cbind(qtnMrkEffTightLinkT1, qtnMrkEffTightLinkT2)
            qtnMrkEffTightLinkAMatAdj <- qtnMrkEffTightLinkAMat %*% cholGenVarMat

            ### define QTN positions
            qtnMrkNameMat <- SampleQtnPosTightLink(linkageInfoDf = linkageInfoDf,
                                                   nChr = nChr,
                                                   nQTNEachChr = nQTNEachChr)
            ### calculate LD to adjust marker effect of Trait2
            qtnCorVec <- apply(X = qtnMrkNameMat, MARGIN = 2, FUN = function(x) {
              cor(x = genoMat[, x[1]], y = genoMat[, x[2]])
            })
            qtnMrkEffTightLinkResVec <- qtnMrkEffTightLinkAMatAdj[, 2] * qtnCorVec

            mrkEffTightLinkMat[qtnMrkNameMat[1, ], 1] <- qtnMrkEffTightLinkAMatAdj[, 1]
            mrkEffTightLinkMat[qtnMrkNameMat[2, ], 2] <- qtnMrkEffTightLinkResVec

            gvT1 <- genoMat %*% mrkEffTightLinkMat[, 1, drop = FALSE] ### genotypic value
            gvT2 <- genoMat %*% mrkEffTightLinkMat[, 2, drop = FALSE] ### genotypic value

            genCorNow <- cor(gvT1, gvT2)
            genCorDiff <- abs(genCorNow - genCor)
            if (genCorDiff <= genCorDiffThre) {
              loopTightLinkVec <- 1
            } else {
              set.seed(qtnSeedBase + match(genCor, genCorVec) * 10000L + cropNo * 1000L + repNo)
              loopTightLinkVec <- 0
            }
          }

          return(mrkEffTightLinkMat)

        }, mc.cores = nCores)
        ### save
        saveRDS(object = mrkEffTightLinkLst, file = fileNameMrkEffNow)
      }
      gc(reset = TRUE); gc(reset = TRUE)
    }
  }

}
