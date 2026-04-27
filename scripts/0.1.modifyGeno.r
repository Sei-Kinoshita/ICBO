###### 0. Settings ######
rm(list = ls())
gc(reset = TRUE);gc(reset = TRUE)

##### 0.1. Load Packages #####
require(data.table)
require(openxlsx)
require(gaston)

##### 0.2. Set Parameters #####
scriptID <- "0.1."
mafThre <- 0.1
ldThre <- 0.6 ### LD threshold for cutting markers
### genome average recombination rate (~cM/Mb)
sbMapRate <- 2.2
gmMapRate <- 2.5

##### 0.3. Set Directories #####
dirRaw <- "rawData/"
dirRawSb <- paste0(dirRaw, "sorghum/")
dirRawGm <- paste0(dirRaw, "soybean/")

dirSave <- "data/genotype/"




###### 1. Load Data ######
### sorghum
fileNameSbGeno <- paste0(dirRawSb, "Sbi454_237_180927.fil.SNP.remHet.MS0.95_bi_MQ20_DP3-1000.MAF0.025.imputed.vcf.gz")
genoSbRaw <- read.vcf(file = fileNameSbGeno)

fileNameSbInfo <- paste0(dirRawSb, "Sbi454_237_180927.fil.SNP.remHet.MS0.95_bi_MQ20_DP3-1000.MAF0.025.imputed_marker.csv")
sbMrkInfo <- fread(fileNameSbInfo)

### soybean (Glycine max)
fileNameGmGeno <- paste0(dirRawGm, "Gm198_HCDB_190207.fil.snp.remHet.MS0.95_bi_MQ20_DP3-1000.MAF0.025.imputed.v2.chrnum.vcf.gz")
genoGmRaw <- read.vcf(fileNameGmGeno)

### selected individuals (by k-medoids)
fileNameSelIndLst <- "data/SelectedIndList_kMedoids_nIndSb_79_nIndGm_66.rds"
selIndLst <- readRDS(file = fileNameSelIndLst)


###### 2. Modify Data ######
##### 2.1. Sorghum #####
### change marker name
genoSbRaw@snps$chr <- sbMrkInfo$chrom
genoSbRaw@snps$id <- paste0("Chr",
                            formatC(sbMrkInfo$chrom, width = 2, flag = "0"),
                            "_",
                            formatC(sbMrkInfo$pos, width = 8, flag = "0"))
### cut by MAF
genoSbMaf <- select.snps(x = genoSbRaw, condition = maf >= mafThre)
### cut by LD
genoSbFil <- as.matrix(LD.thin(genoSbMaf, threshold = ldThre))

### cut hetero markers
isHeteroSb <- apply(X = genoSbFil, MARGIN = 2, 
                    FUN = function(oneMrkVec) any(oneMrkVec == 1))
genoSbMat <- genoSbFil[, !isHeteroSb]


##### 2.2. Soybean #####
### change marker name
genoGmRaw@snps$id <- paste0("Chr",
                            formatC(genoGmRaw@snps$chr, width = 2, flag = "0"),
                            "_",
                            formatC(genoGmRaw@snps$pos, width = 8, flag = "0"))
### cut by MAF
genoGmMaf <- select.snps(x = genoGmRaw, condition = maf >= mafThre)
### cut by LD
genoGmFil <- as.matrix(LD.thin(genoGmMaf, threshold = ldThre))

### cut hetero markers
isHeteroGm <- apply(X = genoGmFil, MARGIN = 2, 
                    FUN = function(oneMrkVec) any(oneMrkVec == 1))
genoGmMat <- genoGmFil[, !isHeteroGm]


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
### marker information 79 lines
fileNameSelIndInfoGm <- paste0(dirSave, scriptID, 
                               "Soybean66_markerInfo_Maf_", mafThre, "_LD_", ldThre, 
                               ".rds")
saveRDS(object = gmSelIndInfoDf, file = fileNameSelIndInfoGm)