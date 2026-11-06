###### 0. Settings ######
rm(list = ls())
gc(reset = TRUE);gc(reset = TRUE)

##### 0.1. Load Packages #####
require(gaston)

##### 0.2. Set Parameters #####
scriptID <- "0.1."
mafThre <- 0.1
ldThre <- 0.6 ### LD threshold for cutting markers
### genome average recombination rate (~cM/Mb)
sbMapRate <- 2.2
gmMapRate <- 2.5

##### 0.3. Set Directories #####
dirSave <- "data/genotype/"
dir.create(dirSave, recursive = TRUE, showWarnings = FALSE)




###### 1. Published genotype panels and fixed selected-accession IDs ######
### The published matrices already underwent full-panel MAF/LD filtering and
### removal of markers containing heterozygotes. See README for raw-data provenance.
selectedFile <- "data/SelectedIndList_kMedoids_nIndSb_79_nIndGm_66.rds"
if (!file.exists(selectedFile)) {
  stop("Selected accession IDs are not distributed with this code-only release. ",
    "Supply ", selectedFile, " as described in README.md.", call. = FALSE)
}
### Existing public full-panel files are at the repository root.
fileSb <- paste0(dirSave, "0.1.Sorghum237_geno_Maf_0.1_LD_0.6.rds")
fileGm <- paste0(dirSave, "0.1.Soybean198_geno_Maf_0.1_LD_0.6.rds")
if (!file.exists(fileSb)) fileSb <- basename(fileSb)
if (!file.exists(fileGm)) fileGm <- basename(fileGm)
genoSbMat <- readRDS(fileSb)
genoGmMat <- readRDS(fileGm)
selIndLst <- readRDS(selectedFile)
stopifnot(nrow(genoSbMat) == 237L, nrow(genoGmMat) == 198L,
  length(selIndLst$Sb) == 79L, length(selIndLst$Gm) == 66L)

###### 2. Marker positions and selected panels ######
##### 2.3. Extract marker position information #####
### sorghum
sbFilInfoDf <- data.frame(MrkName = colnames(genoSbMat), 
                          Chr = as.numeric(substr(colnames(genoSbMat), start = 4, stop = 5)), 
                          Pos = as.numeric(substr(colnames(genoSbMat), start = 7, stop = 14)))
sbFilInfoDf$MapDis <- (sbFilInfoDf$Pos / 1e6) * sbMapRate
sbFirstMrk <- match(unique(sbFilInfoDf$Chr), sbFilInfoDf$Chr) ### first marker for each chr
sbFilInfoDf$MapDis[sbFirstMrk] <- 0
### soybean
gmFilInfoDf <- data.frame(MrkName = colnames(genoGmMat), 
                          Chr = as.numeric(substr(colnames(genoGmMat), start = 4, stop = 5)), 
                          Pos = as.numeric(substr(colnames(genoGmMat), start = 7, stop = 14)))
gmFilInfoDf$MapDis <- (gmFilInfoDf$Pos / 1e6) * gmMapRate
gmFirstMrk <- match(unique(gmFilInfoDf$Chr), gmFilInfoDf$Chr) ### first marker for each chr
gmFilInfoDf$MapDis[gmFirstMrk] <- 0


##### 2.4. Filter to Selected individuals #####
### sorghum
genoSelIndSb <- as.matrix(select.snps(x = as.bed.matrix(genoSbMat[selIndLst$Sb, ]), 
                            condition = maf >= mafThre))
sbSelIndInfoDf <- sbFilInfoDf[sbFilInfoDf$MrkName %in% colnames(genoSelIndSb), ]
sbSelIndFirstMrk <- match(unique(sbSelIndInfoDf$Chr), sbSelIndInfoDf$Chr) ### first marker for each chr
sbSelIndInfoDf$MapDis[sbSelIndFirstMrk] <- 0
### soybean
genoSelIndGm <- as.matrix(select.snps(x = as.bed.matrix(genoGmMat[selIndLst$Gm, ]), 
                                      condition = maf >= mafThre))
gmSelIndInfoDf <- gmFilInfoDf[gmFilInfoDf$MrkName %in% colnames(genoSelIndGm), ]
gmSelIndFirstMrk <- match(unique(gmSelIndInfoDf$Chr), gmSelIndInfoDf$Chr) ### first marker for each chr
gmSelIndInfoDf$MapDis[gmSelIndFirstMrk] <- 0




###### 3. Save ######
### sorghum
### genotype 237 lines
fileNameFilGenoSb <- paste0(dirSave, scriptID, 
                            "Sorghum237_geno_Maf_", mafThre, "_LD_", ldThre, 
                            ".rds")
saveRDS(object = genoSbMat, file = fileNameFilGenoSb)
### marker information 237 lines
fileNameFilInfoSb <- paste0(dirSave, scriptID, 
                            "Sorghum237_markerInfo_Maf_", mafThre, "_LD_", ldThre, 
                            ".rds")
saveRDS(object = sbFilInfoDf, file = fileNameFilInfoSb)
### genotype 79 lines
fileNameSelIndGenoSb <- paste0(dirSave, scriptID, 
                            "Sorghum79_geno_Maf_", mafThre, "_LD_", ldThre, 
                            ".rds")
saveRDS(object = genoSelIndSb, file = fileNameSelIndGenoSb)
### marker information 79 lines
fileNameSelIndInfoSb <- paste0(dirSave, scriptID, 
                            "Sorghum79_markerInfo_Maf_", mafThre, "_LD_", ldThre, 
                            ".rds")
saveRDS(object = sbSelIndInfoDf, file = fileNameSelIndInfoSb)

### soybean
### genotype 198 lines
fileNameFilGenoGm <- paste0(dirSave, scriptID, 
                            "Soybean198_geno_Maf_", mafThre, "_LD_", ldThre, 
                            ".rds")
saveRDS(object = genoGmMat, file = fileNameFilGenoGm)
### marker information 198 lines
fileNameFilInfoGm <- paste0(dirSave, scriptID, 
                            "Soybean198_markerInfo_Maf_", mafThre, "_LD_", ldThre, 
                            ".rds")
saveRDS(object = gmFilInfoDf, file = fileNameFilInfoGm)
### genotype 66 lines
fileNameSelIndGenoGm <- paste0(dirSave, scriptID, 
                               "Soybean66_geno_Maf_", mafThre, "_LD_", ldThre, 
                               ".rds")
saveRDS(object = genoSelIndGm, file = fileNameSelIndGenoGm)
### marker information 66 lines
fileNameSelIndInfoGm <- paste0(dirSave, scriptID, 
                               "Soybean66_markerInfo_Maf_", mafThre, "_LD_", ldThre, 
                               ".rds")
saveRDS(object = gmSelIndInfoDf, file = fileNameSelIndInfoGm)