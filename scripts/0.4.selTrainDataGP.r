###### 0. Settings ######
##### 0.1. Remove objects #####
rm(list = ls())
gc(reset = TRUE);gc(reset = TRUE)

##### 0.2. Load packages #####

##### 0.3. Set parameters #####
scriptID <- "0.4."

cropName <- c("Sb", "Gm")

nCV <- 10
nComb <- 520

nSbM2 <- 26 ### number of selected Sb in method 2
nGmM2 <- 20 ### number of selected Gm in method 2

##### 0.4. Define directories #####
dataDir <- "data/"
genoDir <- paste0(dataDir, "genotype/")

saveDir <- paste0(dataDir, "selTrainData_GP/")
if (!dir.exists(saveDir)) {
  dir.create(saveDir, recursive = TRUE)
}




###### 1. Load Data ######
### Selected Individuals
fileNameSelIndLst <- paste0(dataDir, "SelectedIndList_kMedoids_nIndSb_79_nIndGm_66.rds")
selIndLst <- readRDS(file = fileNameSelIndLst)

### selected individuals' name
sbName <- selIndLst[["Sb"]]
gmName <- selIndLst[["Gm"]]

# ### distance matrix
# fileNameDistSb <- paste0(dirGeno, "0.4.DistMat_euclidean_Sb.rds")
# fileNameDistGm <- paste0(dirGeno, "0.4.DistMat_euclidean_Gm.rds")
# genoDistSb <- readRDS(file = fileNameDistSb)
# genoDistGm <- readRDS(file = fileNameDistGm)




###### 2. Select Training Data ######
##### 2.1. Completely Random #####
fileNameSelCombRandom <- paste0(saveDir, scriptID, 
                                "trainData_selCombName_random.rds")
if (!file.exists(fileNameSelCombRandom)) {
  selCombRandomLst <- lapply(X = 1:nCV, FUN = function(cvNo) {
    
    allCombDf <- expand.grid(sbName, gmName, stringsAsFactors = FALSE)
    colnames(allCombDf) <- c("Sb", "Gm")
    
    combPos <- sample(x = 1:nrow(allCombDf), size = nComb, replace = FALSE)
    finCombDf <- allCombDf[combPos, ]
    
    return(finCombDf)
  })
  names(selCombRandomLst) <- paste0("CV_", 1:nCV)
  saveRDS(object = selCombRandomLst, file = fileNameSelCombRandom)
} else {
  selCombRandomLst <- readRDS(file = fileNameSelCombRandom)
}

##### 2.2. All combinations for selected lines #####
fileNameSelCombAll <- paste0(saveDir, scriptID, "trainData_selCombName_allComb.rds")
if (!file.exists(fileNameSelCombAll)) {
  selCombAllLst <- lapply(X = 1:nCV, FUN = function(cvNo) {
    
    selSbName <- sample(x = sbName, size = nSbM2, replace = FALSE)
    selGmName <- sample(x = gmName, size = nGmM2, replace = FALSE)
    finCombDf <- expand.grid(selSbName, selGmName, stringsAsFactors = FALSE)
    colnames(finCombDf) <- c("Sb", "Gm")
    
    return(finCombDf)
  })
  names(selCombAllLst) <- paste0("CV_", 1:nCV)
  saveRDS(object = selCombAllLst, file = fileNameSelCombAll)
} else {
  selCombAllLst <- readRDS(file = fileNameSelCombAll)
}

##### 2.3. All lines used at least one time #####
fileNameSelCombOnce <- paste0(saveDir, scriptID, "trainData_selCombName_useOnce.rds")
if (!file.exists(fileNameSelCombOnce)) {
  selCombOnceLst <- lapply(X = 1:nCV, FUN = function(cvNo) {
    
    finCombDf <- as.data.frame(matrix(data = NA, nrow = nComb, ncol = 2))
    colnames(finCombDf) <- c("Sb", "Gm")
    
    ### sorghum
    while (any(is.na(finCombDf$Sb))) {
      sbPosVec <- which(is.na(finCombDf$Sb))
      basicSbPos <- sample(x = sbPosVec, size = min(length(sbName), length(sbPosVec)), 
                           replace = FALSE)
      if (length(sbPosVec) < length(sbName)) {
        finCombDf[basicSbPos, "Sb"] <- sample(x = sbName, size = length(sbPosVec), replace = FALSE)
      } else {
        finCombDf[basicSbPos, "Sb"] <- sbName
      }
    }
    
    ### soybean
    while(any(is.na(finCombDf$Gm))) {
      gmPosVec <- which(is.na(finCombDf$Gm))
      basicGmPos <- sample(x = gmPosVec, size = min(length(gmName), length(gmPosVec)), 
                           replace = FALSE)
      if (length(gmPosVec) < length(gmName)) {
        finCombDf[basicGmPos, "Gm"] <- sample(x = gmName, size = length(gmPosVec), replace = FALSE)
      } else {
        finCombDf[basicGmPos, "Gm"] <- gmName
      }
    }
    return(finCombDf)
  })
  names(selCombOnceLst) <- paste0("CV_", 1:nCV)
  saveRDS(object = selCombOnceLst, file = fileNameSelCombOnce)
} else {
  selCombOnceLst <- readRDS(file = fileNameSelCombOnce)
}


# lapply(X = selCombOnceLst, FUN = function(oneCvMat) {
#   table(oneCvMat$Gm)
# })