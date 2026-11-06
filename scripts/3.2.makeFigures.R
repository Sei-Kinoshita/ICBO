###### Manuscript and supplementary figures ######
rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(rrBLUP)
library(multcompView)
source("scripts/function_scripts/2.1.function.BOEval.R")
source("scripts/function_scripts/1.2.function.GPAccuracy.R")
runMCMCChecks <- FALSE

###### Figure 1: BO workflow (manually drawn schematic) ######

###### Figure 2, S7 and S8: predictive ability ######
scriptID <- "1.2."
saveDir <- "midstream/1.2.GP_Accuracy/"
resDf <- readRDS(paste0(saveDir, scriptID, "GP_PublicationSummary.rds"))
figureGroups <- list(Figure2 = c("Sb", "Gm"), FigureS7 = "Sb", FigureS8 = "Gm")
for (figure in names(figureGroups)) {
  pdf(paste0(saveDir, scriptID, figure, "_PredictiveAbility.pdf"), width = 12, height = 9)
  for (h in c(.3, .7)) {
    dat <- resDf[resDf$Crop %in% figureGroups[[figure]] & resDf$h2 == h, ]
    print(MakePlot(dat, paste0("Predictive ability / h2 = ", h)))
  }
  dev.off()
}

###### Figure S3, S5 and S6: realized and estimated genetic correlations ######
PlotGeneticCorrelationRecovery(attr(resDf, "geneticCorrelationRecovery"), saveDir, scriptID)

###### Figure S4: training designs (manually drawn schematic) ######

###### BO figure inputs ######
scriptID <- "2.1."
saveDir <- "midstream/2.1.BO_Evaluation/"
boIterDf <- readRDS(paste0(saveDir, scriptID, "BO_PublicationSummary.rds"))
parSetDf <- readRDS("data/phenotype/0.3.parmeterSet.rds")
selection <- attr(boIterDf, "selectionMechanism")
indexName <- c("Mu", "PI", "EI", "UCB")
nIndex <- length(indexName)
nRep <- nRepBo <- 10L
nPair <- 5214L
nParUniq <- nrow(unique(parSetDf[c("genCor", "h2", "gamma")]))

###### Figure 3: cycle-wise advantage over Mu ######
PlotCycleAdvantage(attr(boIterDf, "cycleAdvantage"),
  paste0(saveDir, scriptID, "Figure3_CycleAdvantage.pdf"))

###### Figures 4, 5 and S10: candidate selection and discovery ranks ######
mechanismMethods <- indexName
mechanismRho <- c(-.7, .7)
mechanismH2 <- .7
mechanismGamma <- .2
mechanismInitialDf <- selection$initial
mechanismExampleDf <- selection$example
mechanismDiscoveryDf <- selection$discovery
mechanismInitialBlocks <- mechanismInitialDf %>% group_by(genCor, phenotypeRep, method) %>%
  summarise(value = mean(differentFromMu), .groups = "drop") %>% mutate(metric = "differentFromMu")
mechanismInitialSummary <- mechanismInitialBlocks %>% group_by(genCor, method, metric) %>%
  summarise(Mean = mean(value), SE = sd(value) / sqrt(n()), .groups = "drop")
mechanismDiscoveryBlocks <- mechanismDiscoveryDf %>% group_by(genCor, h2, gamma, phenotypeRep, method) %>%
  summarise(percent = 100 * mean(outsideMeanTop56), .groups = "drop")
mechanismDiscoverySummary <- mechanismDiscoveryBlocks %>% group_by(genCor, h2, gamma, method) %>%
  summarise(Mean = mean(percent), SE = sd(percent) / sqrt(n()), .groups = "drop")
##### Publication figures #####
mechanismColors <- setNames(scales::hue_pal()(4), mechanismMethods)
mechanismShapes <- c(Mu = 18, PI = 16, EI = 17, UCB = 15)
mechanismTheme <- theme_bw(base_size = 18, base_family = 'sans') +
  theme(plot.title = element_text(size = 23, face = 'bold', margin = margin(b = 9)),
    plot.subtitle = element_text(size = 16, margin = margin(b = 11)),
    axis.title = element_text(size = 19), axis.text = element_text(size = 17),
    axis.ticks = element_line(linewidth = 0.65), panel.border = element_rect(linewidth = 0.75),
    panel.grid.major = element_line(color = 'grey92', linewidth = 0.25),
    panel.grid.minor = element_blank(), strip.background = element_rect(fill = 'grey96'),
    strip.text = element_text(size = 19), panel.spacing = grid::unit(1.0, 'lines'),
    legend.position = 'bottom', legend.title = element_blank(), legend.text = element_text(size = 18),
    legend.key.width = grid::unit(1.3, 'cm'),
    plot.caption = element_text(size = 12.5, hjust = 0, lineheight = 1.12),
    plot.margin = margin(12, 16, 10, 14))
mechanismExampleDf$method <- factor(mechanismExampleDf$method, levels = mechanismMethods)
mechanismInitialBlocks$method <- factor(mechanismInitialBlocks$method, levels = mechanismMethods)
mechanismInitialSummary$method <- factor(mechanismInitialSummary$method, levels = mechanismMethods)
mechanismExampleX <- range(c(mechanismExampleDf$postMean, mechanismExampleDf$incumbent))
mechanismExampleY <- range(mechanismExampleDf$postSD)
###### Figure 4: selection from a common GP ######
mechanismExamplePlots <- list()
for (rhoNow in mechanismRho) {
  exampleNow <- filter(mechanismExampleDf, genCor == rhoNow)
  mechanismExamplePlots[[length(mechanismExamplePlots) + 1L]] <-
    ggplot(exampleNow, aes(postMean, postSD)) +
    geom_point(color = 'grey77', size = 0.65, alpha = 0.4) +
    geom_vline(data = distinct(exampleNow, method, incumbent), aes(xintercept = incumbent),
      linetype = 'dashed', linewidth = 0.8, color = 'grey35') +
    geom_point(data = filter(exampleNow, selected), aes(color = method, shape = method),
      size = 3.2, show.legend = FALSE) +
    geom_point(data = filter(exampleNow, optimum), shape = 21, fill = 'white', color = 'white', size = 5.6) +
    geom_point(data = filter(exampleNow, optimum), shape = 8, color = 'black', size = 5.0, stroke = 1.05) +
    scale_color_manual(values = mechanismColors) + scale_shape_manual(values = mechanismShapes) +
    scale_x_continuous(limits = mechanismExampleX, breaks = scales::breaks_pretty(n = 5),
      expand = expansion(mult = 0.04)) +
    scale_y_continuous(limits = mechanismExampleY, breaks = scales::breaks_pretty(n = 5),
      expand = expansion(mult = 0.04)) +
    facet_wrap(vars(method), ncol = 2) +
    labs(title = 'Figure 4  Selection from a common GP',
      subtitle = bquote(h^2 == 0.7 ~ ',' ~ gamma == 0.2 ~ ',' ~ rho == .(rhoNow) ~
        '; first model-guided batch'),
      x = 'Posterior mean of the objective', y = 'Posterior SD of the latent objective',
      caption = paste0('Gray: all 5,158 candidates; colored points: 56 selected pairs; black star: the true optimum.\n',
        'Dashed line: observed-pair posterior-mean incumbent (PI/EI). Both scenarios use the same axes.\n',
        'Phenotype replicate 1, initial set 1. The common GP reproduces all four saved selections.\n',
        'SD is conditional on fitted variance components and excludes observation noise.')) +
    mechanismTheme + theme(legend.position = 'none')
}

###### Figure 5b: selected pairs differing from Mu ######
mechanismSummaryPlots <- list()
for (metricNow in 'differentFromMu') {
  blockNow <- filter(mechanismInitialBlocks, metric == metricNow) %>%
    mutate(displayMethod = as.numeric(method) + (phenotypeRep - (nRep + 1) / 2) * 0.035)
  summaryNow <- filter(mechanismInitialSummary, metric == metricNow) %>%
    mutate(displayMethod = as.numeric(method))
  yLabel <- 'Pairs differing from Mu'
  p <- ggplot(summaryNow, aes(displayMethod, Mean, color = method, shape = method)) +
    geom_point(data = blockNow, aes(y = value), size = 2.3, alpha = 0.38, show.legend = FALSE) +
    geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.16, linewidth = 1.05) +
    geom_point(size = 4.6) +
    scale_color_manual(values = mechanismColors, drop = FALSE) +
    scale_shape_manual(values = mechanismShapes, drop = FALSE) +
    scale_x_continuous(breaks = 1:4, labels = mechanismMethods, limits = c(0.65, 4.35)) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 4), expand = expansion(mult = 0.15)) +
    facet_wrap(vars(genCor), nrow = 1,
      labeller = label_bquote(rho == .(as.numeric(as.character(genCor))))) +
    labs(x = NULL, y = yLabel) + mechanismTheme + theme(legend.position = 'none')
  if (metricNow == 'differentFromMu') {
    p <- p + labs(title = 'Figure 5b  Differences in the first selected batch',
      subtitle = expression(h^2 == 0.7 ~ ',' ~ gamma == 0.2),
      caption = paste0('Large points: mean +/- SE of 10 phenotype means, each averaging 10 initial-set runs.\n',
      'Small points: phenotype means; differing pairs = 56 - overlap with Mu.'))
  }
  mechanismSummaryPlots[[length(mechanismSummaryPlots) + 1L]] <- p
}

###### Figure 5a: posterior-mean and acquisition ranks at discovery ######
mechanismRankDf <- mechanismDiscoveryDf %>%
  filter(h2 == mechanismH2, gamma == mechanismGamma, genCor %in% mechanismRho, method != 'Mu') %>%
  mutate(method = factor(method, levels = c('PI', 'EI', 'UCB')),
    genCor = factor(genCor, levels = mechanismRho))
mechanismRankLabels <- mechanismRankDf %>% group_by(genCor, method) %>%
  summarise(label = paste0(sum(outsideMeanTop56), '/100\noutside mean top 56'), .groups = 'drop')
mechanismRankXMax <- max(350, ceiling(max(mechanismRankDf$meanRank) / 100) * 100)
mechanismRankPlot <- ggplot(mechanismRankDf,
  aes(meanRank, acquisitionRank, color = method, shape = method)) +
  annotate('rect', xmin = 56.5, xmax = mechanismRankXMax, ymin = 0.85, ymax = 140, fill = 'grey93') +
  geom_abline(intercept = 0, slope = 1, color = 'grey70', linewidth = 0.65) +
  geom_vline(xintercept = 56.5, linetype = 'dashed', color = 'grey35', linewidth = 0.8) +
  geom_hline(yintercept = 56.5, linetype = 'dashed', color = 'grey35', linewidth = 0.8) +
  geom_point(size = 3.2, alpha = 0.7) +
  geom_text(data = mechanismRankLabels, aes(x = 1, y = 125, label = label),
    inherit.aes = FALSE, hjust = 0, vjust = 1, size = 4.5, lineheight = 0.95, color = 'grey20') +
  scale_color_manual(values = mechanismColors) + scale_shape_manual(values = mechanismShapes) +
  scale_x_log10(breaks = c(1, 10, 56, 300, 1000, 3000), limits = c(0.85, mechanismRankXMax), expand = expansion(mult = 0.01)) +
  scale_y_log10(breaks = c(1, 10, 56), limits = c(0.85, 140), expand = expansion(mult = 0.01)) +
  facet_grid(rows = vars(genCor), cols = vars(method),
    labeller = label_bquote(rows = rho == .(as.numeric(as.character(genCor))))) +
  labs(title = 'Figure 5a  Ranking the optimum at discovery',
    subtitle = expression(h^2 == 0.7 ~ ',' ~ gamma == 0.2 ~ '; ranks under the same GP at discovery'),
    x = 'Rank by posterior mean (log scale)',
    y = 'Rank by acquisition score (log scale)',
    caption = paste0('Each point is one run (100 per panel). Shading: posterior-mean rank > 56.\n',
      'Both ranks use the same GP and unobserved candidate set immediately before discovery.\n',
      'Dashed lines: batch cutoff; gray diagonal: equal ranks. Axes are logarithmic (1 = best).\n',
      'Methods have different histories; the comparison concerns their final selection decisions.\n',
      'Mu has identical ranks and is omitted. Overlapping points are not jittered.')) +
  mechanismTheme + theme(legend.position = 'none')

###### Figure S10: optimum outside the posterior-mean top 56 ######
mechanismScenarioLabels <- mechanismDiscoverySummary %>% distinct(genCor, gamma) %>%
  arrange(genCor, gamma) %>% mutate(label = sprintf('rho == %.1f ~ "," ~ gamma == %.1f', genCor, gamma))
mechanismAllDf <- mechanismDiscoverySummary %>% filter(method != 'Mu') %>%
  mutate(method = factor(method, levels = c('PI', 'EI', 'UCB')),
    scenario = factor(sprintf('rho == %.1f ~ "," ~ gamma == %.1f', genCor, gamma),
      levels = rev(mechanismScenarioLabels$label)))
mechanismAllDodge <- position_dodge(width = 0.9, orientation = 'y', reverse = TRUE)
mechanismRateLimit <- max(40, ceiling(max(mechanismAllDf$Mean + mechanismAllDf$SE) / 10) * 10)
mechanismAllPlot <- ggplot(mechanismAllDf, aes(Mean, scenario, color = method, shape = method, group = method)) +
  geom_vline(xintercept = 0, color = 'grey55', linewidth = 0.5) +
  geom_errorbar(aes(xmin = Mean - SE, xmax = Mean + SE), orientation = 'y',
    position = mechanismAllDodge, width = 0.16, linewidth = 1.05) +
  geom_point(position = mechanismAllDodge, size = 4.4) +
  scale_color_manual(values = mechanismColors, breaks = c('PI', 'EI', 'UCB')) +
  scale_shape_manual(values = mechanismShapes, breaks = c('PI', 'EI', 'UCB')) +
  scale_x_continuous(limits = c(0, mechanismRateLimit), breaks = seq(0, mechanismRateLimit, by = 10),
    labels = function(x) paste0(x, '%'), expand = expansion(mult = c(0.02, 0.03))) +
  scale_y_discrete(labels = function(x) parse(text = x), expand = expansion(add = 0.5)) +
  facet_wrap(vars(h2), nrow = 1, labeller = label_bquote(h^2 == .(h2))) +
  labs(title = 'Figure S10  Optima outside the mean top 56',
    subtitle = 'All 18 scenarios; ranks under the same GP at discovery',
    x = 'Runs with the optimum outside the mean top 56', y = NULL,
    caption = paste0('Points: percentage of 100 runs. Bars: SE of 10 phenotype means (10 initial sets each).\n',
      'The optimum was selected despite a posterior-mean rank > 56 under that method\'s GP.\n',
      'Mu is 0% by construction. Earlier data collection differs among methods.\n',
      'Rows: target DGE-IGE correlation and SMA/GMA variance ratio.')) +
  mechanismTheme + theme(panel.grid.major.y = element_line(color = 'grey92', linewidth = 0.25),
    panel.spacing.x = grid::unit(1.25, 'cm'))


cairo_pdf(paste0(saveDir, scriptID, "Figure4_CommonGPSelection.pdf"), width = 11.5, height = 9.5)
for (p in mechanismExamplePlots) print(p)
dev.off()
cairo_pdf(paste0(saveDir, scriptID, "Figure5_DiscoveryRanks_SelectionDifference.pdf"), width = 11.5, height = 9.5)
print(mechanismRankPlot)
print(mechanismSummaryPlots[[1]])
dev.off()
cairo_pdf(paste0(saveDir, scriptID, "FigureS10_OptimumOutsideMeanTop56.pdf"), width = 11.5, height = 9.5)
print(mechanismAllPlot)
dev.off()

###### Figure S9: conditional random reference and paired cycle savings ######
mainRunDf <- boIterDf %>% transmute(genCor = GenCor, h2 = h2, gamma = Gamma,
  phenotypeRep = parSetDf$rep[ParNo], ParNo = ParNo, method = Index, boRep = BoRepNo,
  nCycles = nIter, nEvaluatedPair = pmin(56L * nIter, nPair))
mainMethodNames <- c("Random", indexName)
stopifnot(nrow(mainRunDf) == 7200L, all(attr(boIterDf, "runTiming")$optimumFound))
mainPhenotypeDf <- mainRunDf %>%
  group_by(genCor, h2, gamma, phenotypeRep, ParNo, method) %>%
  summarise(meanCycles = mean(nCycles),
    nInitialSetRep = n(), .groups = "drop")
stopifnot(all(mainPhenotypeDf$nInitialSetRep == nRepBo))

### Evaluate empirical GP discovery rates and the analytical Random reference at identical batch budgets.
### Balanced run counts give equal weight to the 18 scenarios and to their phenotype replicates.
mainRunCounts <- mainRunDf %>% count(genCor, h2, gamma, method)
stopifnot(nrow(mainRunCounts) == nParUniq * nIndex,
  all(mainRunCounts$n == nRep * nRepBo))
mainBudgetVec <- sort(unique(c(0L, seq(56L, nPair, by = 56L), nPair)))
mainDiscoveryDf <- bind_rows(lapply(indexName, function(methodNow) {
  completedAt <- mainRunDf$nEvaluatedPair[mainRunDf$method == methodNow]
  data.frame(method = methodNow, evaluatedPairs = mainBudgetVec,
    discoveryPercent = 100 * vapply(mainBudgetVec, function(budget) {
      mean(completedAt <= budget)
    }, numeric(1)))
}))
### The unique optimum is excluded from the initial 56 pairs. Its rank among the
### remaining pairs is uniform under sampling without replacement. Include the final partial batch.
mainInitialPairs <- mainBatchSize <- 56L
mainRemainingPairs <- nPair - mainInitialPairs
stopifnot(mainRemainingPairs > 0)
mainRandomExpectedCycles <- 1 + mean(ceiling(seq_len(mainRemainingPairs) / mainBatchSize))
mainRandomDiscoveryDf <- data.frame(method = "Random", evaluatedPairs = mainBudgetVec,
  discoveryPercent = 100 * pmin(pmax(mainBudgetVec - mainInitialPairs, 0) / mainRemainingPairs, 1))
mainDiscoveryDf <- bind_rows(mainRandomDiscoveryDf, mainDiscoveryDf)
mainDiscoveryDf$method <- factor(mainDiscoveryDf$method, levels = mainMethodNames)

### Match methods within phenotype replicate before estimating mean cycle savings and intervals.
mainMuDf <- mainPhenotypeDf %>% filter(method == "Mu") %>%
  select(genCor, h2, gamma, phenotypeRep, ParNo, muMeanCycles = meanCycles)
mainPairedDf <- mainPhenotypeDf %>% filter(method %in% setdiff(indexName, "Mu")) %>%
  inner_join(mainMuDf, by = c("genCor", "h2", "gamma", "phenotypeRep", "ParNo")) %>%
  mutate(cyclesSaved = muMeanCycles - meanCycles)
mainContrastDf <- mainPairedDf %>%
  group_by(genCor, h2, gamma, method) %>%
  summarise(Mean = mean(cyclesSaved), SE = sd(cyclesSaved) / sqrt(n()),
    nPhenotypeRep = n(), .groups = "drop") %>%
  mutate(CILower = ifelse(SE > 0, Mean - qt(0.975, nPhenotypeRep - 1) * SE, NA_real_),
    CIUpper = ifelse(SE > 0, Mean + qt(0.975, nPhenotypeRep - 1) * SE, NA_real_))
stopifnot(nrow(mainContrastDf) == nParUniq * (nIndex - 1L),
  all(mainContrastDf$nPhenotypeRep == nRep))

mainColors <- c(Random = "grey40", setNames(scales::hue_pal()(nIndex), indexName))
mainLineTypes <- c(Random = "longdash", Mu = "solid", PI = "dashed", EI = "dotted", UCB = "dotdash")
mainMethodLabels <- c(Random = "Random (analytical)", setNames(indexName, indexName))
mainScenarioLabels <- mainContrastDf %>% distinct(genCor, gamma) %>% arrange(genCor, gamma) %>%
  mutate(label = sprintf("rho == %.1f ~ ',' ~ gamma == %.1f", genCor, gamma))
mainContrastDf <- mainContrastDf %>%
  mutate(scenario = factor(sprintf("rho == %.1f ~ ',' ~ gamma == %.1f", genCor, gamma),
                          levels = rev(mainScenarioLabels$label)),
    h2 = factor(h2, levels = c(0.3, 0.7)), method = factor(method, levels = setdiff(indexName, "Mu")))
mainTheme <- theme_bw(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 19),
    plot.subtitle = element_text(size = 14, margin = margin(b = 10)),
    axis.title = element_text(size = 18), axis.text = element_text(size = 16),
    axis.ticks = element_line(linewidth = 0.6), panel.border = element_rect(linewidth = 0.7),
    panel.grid.major = element_line(linewidth = 0.3, color = "grey90"),
    panel.grid.minor = element_blank(), legend.text = element_text(size = 16),
    legend.title = element_text(size = 16), legend.key.width = grid::unit(1.8, "cm"),
    plot.caption = element_text(size = 11, hjust = 0, lineheight = 1.15),
    plot.margin = margin(12, 18, 10, 12))

mainPlotABase <- ggplot(mainDiscoveryDf,
  aes(evaluatedPairs, discoveryPercent, color = method, linetype = method, group = method)) +
  geom_step(linewidth = 1, direction = "hv") +
  scale_color_manual(values = mainColors, breaks = mainMethodNames, labels = mainMethodLabels, drop = FALSE) +
  scale_linetype_manual(values = mainLineTypes, breaks = mainMethodNames, labels = mainMethodLabels, drop = FALSE) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20),
    labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0.015, 0.025))) +
  labs(x = "Evaluated pairs", y = "Optimum discovery rate", color = "Method", linetype = "Method") +
  mainTheme + theme(legend.position = "bottom")
mainPlotAFull <- mainPlotABase +
  scale_x_continuous(breaks = c(0, 1000, 2000, 3000, 4000, nPair),
    labels = scales::label_comma(), expand = expansion(mult = c(0, 0.015))) +
  labs(title = "A  Optimum discovery at a common evaluation budget",
    subtitle = "Full range: empirical GP discovery rates and an analytical Random reference") +
  theme(legend.position = "bottom")
mainDodge <- position_dodge(width = 0.65, orientation = "y", reverse = TRUE)
mainPlotB <- ggplot(mainContrastDf, aes(x = Mean, y = scenario, color = method, group = method)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.75) +
  geom_errorbar(aes(xmin = CILower, xmax = CIUpper), orientation = "y",
    position = mainDodge, width = 0.16, linewidth = 0.8, na.rm = TRUE) +
  geom_point(aes(shape = method), position = mainDodge, size = 3.4) +
  scale_color_manual(values = mainColors, breaks = setdiff(indexName, "Mu")) +
  scale_shape_manual(values = c(PI = 16, EI = 17, UCB = 15)) +
  scale_y_discrete(labels = function(x) parse(text = x), expand = expansion(add = 0.6)) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 6)) +
  facet_wrap(vars(h2), nrow = 1,
    labeller = label_bquote(h^2 == .(as.numeric(as.character(h2))))) +
  labs(title = "Figure S9b  Cycles saved relative to Mu",
    subtitle = "Positive values favor the acquisition method; all 18 scenarios are shown",
    x = "Cycles saved relative to Mu", y = NULL, color = "Method", shape = "Method",
    caption = paste0("Paired mean differences and pointwise 95% t intervals across 10 phenotype replicates; 10 initial-set runs averaged within replicate.\n",
      "Initial batch = cycle 1; 56 pairs per cycle. Scenario labels: target DGE-IGE correlation and SMA/GMA variance ratio.\n",
      "Intervals are conditional on the initial sets. The true optimum is a simulation evaluation endpoint.")) +
  mainTheme + theme(legend.position = "bottom", strip.text = element_text(size = 17),
    panel.spacing.x = grid::unit(1.2, "lines"), panel.grid.major.y = element_line(color = "grey92"))


cairo_pdf(paste0(saveDir, scriptID, "FigureS9_Discovery_CycleSavings.pdf"), width = 13, height = 10)
print(mainPlotAFull + labs(title = "Figure S9a  Optimum discovery"))
print(mainPlotB)
dev.off()


###### Figure S1 and S2: PCA, MAF and genomic relationships ######
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
  stop("Missing selected accession IDs: ", selectedFile, call. = FALSE)
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

###### Optional reviewer-response figures: MCMC running means ######
if (runMCMCChecks) {
  parameters <- readRDS("data/phenotype/0.3.parmeterSet.rds")
  modelDirs <- c("Model1_DGE", "Model2_DGE_IGE", "Model3_GMA_SMA",
    "Model4_Multi_DGE_IGE", "Model5_Multi_GMA_SMA")
  gp <- EvaluateGpRunningMeans(parameters, modelDirs, "midstream/1.1.GP_CV/")
  gp$plotData <- filter(gp$plotData, Model == "Model3", trainingDesign == "All")
  PlotMcmcRunningMeans(gp, "midstream/ResponseFigure1_GP_RunningMeans.pdf", "GP")
  bo <- EvaluateBoRunningMeans(parameters, "midstream/2.0.BO/Model3_GMA_SMA/", cycles = c(1L, 3L))
  bo$plotData <- filter(bo$plotData, method == "EI", h2 == .7, gamma == .2)
  PlotMcmcRunningMeans(bo, "midstream/ResponseFigure2_BO_RunningMeans.pdf", "BO")
}
