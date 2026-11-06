###### 0. Settings ######
### Read saved genotype panels; do not rerun preprocessing or change GP inputs.
require(ggplot2)
require(rrBLUP)

scriptID <- "0.1."
mafThre <- 0.1
dirGeno <- "data/genotype/"
saveDir <- "midstream/"
dir.create(saveDir, recursive = TRUE, showWarnings = FALSE)

cropName <- c("Sb", "Gm")
cropLabel <- c("Sorghum", "Soybean")
nFull <- c(237L, 198L)
nSelected <- c(79L, 66L)
selectedFile <- "data/SelectedIndList_kMedoids_nIndSb_79_nIndGm_66.rds"
if (!file.exists(selectedFile)) {
  stop("Selected accession IDs are not distributed with this code-only release. ",
    "Supply ", selectedFile, " and run scripts/0.1.modifyGeno.r first; see README.md.",
    call. = FALSE)
}
selected <- readRDS(selectedFile)
panelLevels <- c("Full panel", "Selected: same SNPs", "Selected: final SNPs")
panelColors <- c("Full panel" = "#0072B2", "Selected: same SNPs" = "#D55E00",
                 "Selected: final SNPs" = "#009E73")
plotLst <- pcaPlots <- list()

for (cropNo in seq_along(cropName)) {
  crop <- cropName[cropNo]
  label <- cropLabel[cropNo]
  fileFull <- paste0(dirGeno, scriptID, label, nFull[cropNo], "_geno_Maf_0.1_LD_0.6.rds")
  fileFinal <- paste0(dirGeno, scriptID, label, nSelected[cropNo], "_geno_Maf_0.1_LD_0.6.rds")
  full <- readRDS(fileFull)
  final <- readRDS(fileFinal)
  selectedIds <- selected[[crop]]
  stopifnot(nrow(full) == nFull[cropNo], nrow(final) == nSelected[cropNo],
    !anyNA(full), !anyNA(final), !anyDuplicated(rownames(full)),
    !anyDuplicated(colnames(full)), !anyDuplicated(selectedIds),
    all(selectedIds %in% rownames(full)), identical(rownames(final), selectedIds))
  subset <- full[selectedIds, , drop = FALSE]
  stopifnot(identical(subset[, colnames(final), drop = FALSE], final))

  ### Keep every full-panel SNP when measuring the effect of accession selection alone.
  pFull <- colMeans(full) / 2
  pSubset <- colMeans(subset) / 2
  mafFull <- pmin(pFull, 1 - pFull)
  mafSubset <- pmin(pSubset, 1 - pSubset)
  kept <- colnames(full) %in% colnames(final)
  stopifnot(all(mafFull >= mafThre), identical(kept, unname(mafSubset >= mafThre)))
  mafFinal <- mafSubset[kept]
  ### The primary GRM comparison fixes both SNPs and the full-panel reference frequencies.
  gFull <- rrBLUP::A.mat(full - 1, n.core = 1)
  gFixed <- gFull[selectedIds, selectedIds, drop = FALSE]
  ### Separate the additional effects of re-centering and the repeated MAF filter.
  gSubset <- rrBLUP::A.mat(subset - 1, n.core = 1)
  gFinal <- rrBLUP::A.mat(final - 1, n.core = 1)
  offFull <- gFull[lower.tri(gFull)]
  offFixed <- gFixed[lower.tri(gFixed)]
  offSubset <- gSubset[lower.tri(gSubset)]
  offFinal <- gFinal[lower.tri(gFinal)]
  offDf <- rbind(
    data.frame(Panel = panelLevels[1], Basis = "Full-panel reference", G = offFull),
    data.frame(Panel = panelLevels[2], Basis = "Full-panel reference", G = offFixed),
    data.frame(Panel = panelLevels[1], Basis = "Panel-specific reference", G = offFull),
    data.frame(Panel = panelLevels[2], Basis = "Panel-specific reference", G = offSubset),
    data.frame(Panel = panelLevels[3], Basis = "Panel-specific reference", G = offFinal))

  ### Common MAF bins, with each panel normalized by its own SNP count.
  breaks <- seq(0, 0.5, by = 0.025)
  mafHist <- do.call(rbind, lapply(seq_along(panelLevels), function(i) {
    values <- list(mafFull, mafSubset, mafFinal)[[i]]
    h <- hist(values, breaks = breaks, plot = FALSE, right = FALSE, include.lowest = TRUE)
    data.frame(Panel = panelLevels[i], MAF = c(h$breaks, tail(h$breaks, 1)),
      Proportion = c(h$counts / length(values), 0, 0))
  }))
  mafHist$Panel <- factor(mafHist$Panel, levels = panelLevels)
  offDf$Panel <- factor(offDf$Panel, levels = panelLevels)
  pMAF <- ggplot(mafHist, aes(MAF, Proportion, color = Panel)) +
    geom_step(linewidth = 0.8) +
    geom_vline(xintercept = mafThre, linetype = "dotted", color = "grey40") +
    scale_color_manual(values = panelColors, drop = FALSE) +
    scale_x_continuous(limits = c(0, 0.5)) +
    labs(title = paste0(label, ": MAF spectrum"),
      subtitle = paste(nrow(full), "to", nrow(subset), "accessions; bin width = 0.025"),
      x = "Minor allele frequency", y = "Proportion of SNPs", color = NULL) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom", legend.text = element_text(size = 8)) +
    guides(color = guide_legend(nrow = 2))
  pFixed <- ggplot(offDf[offDf$Basis == "Full-panel reference", ], aes(G, color = Panel)) +
    stat_ecdf(geom = "step", linewidth = 0.8) + scale_color_manual(values = panelColors) +
    labs(title = paste0(label, ": GRM, fixed reference"),
      subtitle = "Same SNPs and full-panel allele frequencies",
      x = "Off-diagonal genomic relationship", y = "Cumulative proportion of pairs", color = NULL) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom", legend.text = element_text(size = 8))
  pRecentered <- ggplot(offDf[offDf$Basis == "Panel-specific reference", ], aes(G, color = Panel)) +
    stat_ecdf(geom = "step", linewidth = 0.8) + scale_color_manual(values = panelColors) +
    labs(title = paste0(label, ": GRM, re-estimated reference"),
      subtitle = "Allele frequencies re-estimated within each panel",
      x = "Off-diagonal genomic relationship", y = "Cumulative proportion of pairs", color = NULL) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom", legend.text = element_text(size = 8)) +
    guides(color = guide_legend(nrow = 2))
  plotLst[[crop]] <- list(pMAF, pFixed, pRecentered)
  ### Figure S1: PC coordinates from the full-panel GRM eigendecomposition.
  eig <- eigen(gFull, symmetric = TRUE)
  scores <- sweep(eig$vectors[, 1:2, drop = FALSE], 2, sqrt(pmax(eig$values[1:2], 0)), "*")
  explained <- 100 * pmax(eig$values[1:2], 0) / sum(pmax(eig$values, 0))
  pcaDf <- data.frame(PC1 = scores[, 1], PC2 = scores[, 2], Selected = rownames(gFull) %in% selectedIds)
  pcaPlots[[crop]] <- ggplot(pcaDf, aes(PC1, PC2, color = Selected)) +
    geom_point(size = 2.6) + scale_color_manual(values = c("FALSE" = "grey65", "TRUE" = "#F8766D"),
      labels = c("Not selected", "Selected")) +
    labs(title = label, x = sprintf("PC1 (%.1f%%)", explained[1]), y = sprintf("PC2 (%.1f%%)", explained[2]), color = NULL) +
    theme_bw(base_size = 16) + theme(legend.position = "bottom")
  rm(full, final, subset, gFull, gFixed, gSubset, gFinal)
  gc()
}

###### 3. Figures S1 and S2 ######
cairo_pdf(paste0(saveDir, scriptID, "FigureS1_PCA.pdf"), width = 12, height = 6)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
for (i in 1:2) print(pcaPlots[[i]], newpage = FALSE,
  vp = grid::viewport(layout.pos.row = 1, layout.pos.col = i))
dev.off()
cairo_pdf(paste0(saveDir, scriptID, "FigureS2_SubsetComparison.pdf"), width = 13, height = 13)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(3, 2)))
for (row in 1:3) for (column in 1:2) print(plotLst[[column]][[row]], newpage = FALSE,
  vp = grid::viewport(layout.pos.row = row, layout.pos.col = column))
dev.off()
