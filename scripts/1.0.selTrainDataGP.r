###### 1. Settings ######
##### 1.1. Remove objects #####
rm(list = ls())
gc(reset = TRUE);gc(reset = TRUE)

##### 1.2. Load packages #####
source("scripts/function_scripts/2.0.function.BO.R")

##### 1.3. Set parameters #####
scriptID <- "0.4." ### Training-set file prefix

cropName <- c("Sb", "Gm")

nCV <- 10
nComb <- 520
cv2Seed <- 20260928

nSbCV1 <- 26 ### number of selected Sb in CV1
nGmCV1 <- 20 ### number of selected Gm in CV1

##### 1.4. Define directories #####
dataDir <- "data/"
genoDir <- paste0(dataDir, "genotype/")

saveDir <- paste0(dataDir, "selTrainData_GP/")
if (!dir.exists(saveDir)) {
  dir.create(saveDir, recursive = TRUE)
}

###### 2. Load Data ######
### Selected Individuals
fileNameSelIndLst <- paste0(dataDir, "SelectedIndList_kMedoids_nIndSb_79_nIndGm_66.rds")
selIndLst <- readRDS(file = fileNameSelIndLst)

### selected individuals' name
sbName <- selIndLst[["Sb"]]
gmName <- selIndLst[["Gm"]]

###### 3. Select Training Data ######
##### 3.1. Design-D: all combinations for selected lines #####
fileNameSelCombAll <- paste0(saveDir, scriptID, "trainData_selCombName_allComb.rds")
if (!file.exists(fileNameSelCombAll)) {
  ### Seed for generating Design-D training sets.
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

##### 3.2. Design-B: balanced coverage of all lines with unique pairs #####
fileNameSelCombOnce <- paste0(saveDir, scriptID,
                             "trainData_selCombName_useOnce.rds")
selCombOnceLst <- if (file.exists(fileNameSelCombOnce)) readRDS(fileNameSelCombOnce) else NULL
### Reuse training sets only when marked as balanced unique pairs.
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

###### 4. Check training-set size, uniqueness and balance ######
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

