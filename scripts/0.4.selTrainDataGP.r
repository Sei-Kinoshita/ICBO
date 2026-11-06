###### 0. Settings ######
##### 0.1. Remove objects #####
rm(list = ls())
gc(reset = TRUE);gc(reset = TRUE)

##### 0.2. Load packages #####
source("scripts/function_scripts/2.0.function.BO.R")

##### 0.3. Set parameters #####
scriptID <- "0.4."

cropName <- c("Sb", "Gm")

nCV <- 10
nComb <- 520
cv2Seed <- 20260928

nSbCV1 <- 26 ### number of selected Sb in CV1
nGmCV1 <- 20 ### number of selected Gm in CV1

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




###### 2. Select Training Data ######
##### 2.1. Design-D: all combinations for selected lines #####
fileNameSelCombAll <- paste0(saveDir, scriptID, "trainData_selCombName_allComb.rds")
if (!file.exists(fileNameSelCombAll)) {
  ### Historical Design-D sets may be supplied locally; this seed is for new sets only.
  set.seed(20261106L)
  selCombAllLst <- lapply(X = 1:nCV, FUN = function(cvNo) {
    
    selSbName <- sample(x = sbName, size = nSbCV1, replace = FALSE)
    selGmName <- sample(x = gmName, size = nGmCV1, replace = FALSE)
    finCombDf <- expand.grid(selSbName, selGmName, stringsAsFactors = FALSE)
    colnames(finCombDf) <- c("Sb", "Gm")
    
    return(finCombDf)
  })
  names(selCombAllLst) <- paste0("CV_", 1:nCV)
  saveRDS(object = selCombAllLst, file = fileNameSelCombAll)
} else {
  selCombAllLst <- readRDS(file = fileNameSelCombAll)
}

##### 2.2. Design-B: balanced coverage of all lines with unique pairs #####
fileNameSelCombOnce <- paste0(saveDir, scriptID,
                             "trainData_selCombName_useOnce.rds")
selCombOnceLst <- if (file.exists(fileNameSelCombOnce)) readRDS(fileNameSelCombOnce) else NULL
### Replace the legacy design under the original filename instead of reusing duplicate pairs.
if (!identical(attr(selCombOnceLst, "design"), "balanced_unique_pairs")) {
  allCombDf <- expand.grid(Sb = sbName, Gm = gmName, stringsAsFactors = FALSE)
  selCombOnceLst <- lapply(X = 1:nCV, FUN = function(cvNo) {
    ### Balance both crops' marginal counts while selecting every pair at most once.
    set.seed(cv2Seed + cvNo)
    combPos <- SelectBroadCoverage(pairDf = allCombDf, trainingSize = nComb)
    finCombDf <- allCombDf[combPos, ]
    rownames(finCombDf) <- NULL
    return(finCombDf)
  })
  names(selCombOnceLst) <- paste0("CV_", 1:nCV)
  attr(selCombOnceLst, "seed") <- setNames(cv2Seed + seq_len(nCV), names(selCombOnceLst))
  attr(selCombOnceLst, "design") <- "balanced_unique_pairs"
  saveRDS(object = selCombOnceLst, file = fileNameSelCombOnce)
} else {
  selCombOnceLst <- readRDS(file = fileNameSelCombOnce)
}

### Check generated or reused training sets against the CV2 design.
stopifnot(length(selCombOnceLst) == nCV)
for (finCombDf in selCombOnceLst) {
  sbCount <- table(factor(finCombDf$Sb, levels = sbName))
  gmCount <- table(factor(finCombDf$Gm, levels = gmName))
  stopifnot(nrow(finCombDf) == nComb,
            !anyDuplicated(finCombDf[, c("Sb", "Gm")]),
            all(finCombDf$Sb %in% sbName), all(finCombDf$Gm %in% gmName),
            all(sbCount >= floor(nComb / length(sbName))),
            all(sbCount <= ceiling(nComb / length(sbName))),
            all(gmCount >= floor(nComb / length(gmName))),
            all(gmCount <= ceiling(nComb / length(gmName))))
}


