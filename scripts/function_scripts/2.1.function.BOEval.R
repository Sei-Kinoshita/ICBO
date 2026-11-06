###### Running means from saved parameter draws (no model refitting) ######
ReadSavedMcmcTrace <- function(path, nColumn = 1L, nIter = 10000L, thin = 20L) {
  values <- scan(path, quiet = TRUE)
  if (length(values) != (nIter / thin) * nColumn || any(!is.finite(values)))
    stop("Incomplete or non-finite MCMC trace: ", path)
  matrix(values, ncol = nColumn, byrow = TRUE)
}

RevisionPlotTheme <- function() {
  theme_bw(base_size = 18) + theme(plot.title = element_text(size = 22, face = "bold"),
    plot.subtitle = element_text(size = 15), axis.title = element_text(size = 19),
    axis.text = element_text(size = 16), strip.text = element_text(size = 17),
    legend.position = "bottom", legend.text = element_text(size = 17),
    plot.caption = element_text(size = 12, hjust = 0),
    panel.spacing = grid::unit(1, "lines"), panel.grid.minor = element_blank(),
    plot.margin = margin(t = 12, r = 16, b = 12, l = 18, unit = "pt"))
}

McmcExampleScenarios <- function(parSetDf) {
  ### Fixed scientific examples, not selected for apparent chain stability.
  selected <- which(parSetDf$genCor == -.7 & parSetDf$rep == 1L &
    ((parSetDf$h2 == .3 & parSetDf$gamma == .4) |
     (parSetDf$h2 == .7 & parSetDf$gamma == .2)))
  stopifnot(length(selected) == 2L)
  selected[order(parSetDf$h2[selected])]
}

ReadMcmcRunningMeans <- function(manifest, nIter = 10000L, burnIn = 2000L, thin = 20L) {
  stopifnot(nIter > burnIn, nIter %% thin == 0L, burnIn %% thin == 0L)
  rows <- summaries <- vector("list", nrow(manifest))
  for (i in seq_len(nrow(manifest))) {
    spec <- manifest[i, ]
    values <- ReadSavedMcmcTrace(spec$sourceFile, spec$nColumn, nIter, thin)
    if (spec$traceColumn == 0L) {
      stopifnot(all(values[, 1] > 0), all(values[, 3] > 0),
        all(values[, 1] * values[, 3] > values[, 2]^2))
      value <- values[, 2] / sqrt(values[, 1] * values[, 3])
    } else value <- values[, spec$traceColumn]
    iteration <- seq_along(value) * thin
    value <- value[iteration > burnIn]
    iteration <- iteration[iteration > burnIn]
    runningMean <- cumsum(value) / seq_along(value)
    stopifnot(abs(tail(runningMean, 1) - mean(value)) < 1e-12)
    if ("expectedMean" %in% names(spec) && is.finite(spec$expectedMean))
      stopifnot(abs(mean(value) - spec$expectedMean) < 1e-6)
    rows[[i]] <- cbind(spec[rep(1L, length(value)), , drop = FALSE],
      data.frame(Iteration = iteration, RunningMean = runningMean))
    rownames(rows[[i]]) <- NULL
    half <- length(value) / 2L
    quarter <- 3L * length(value) / 4L
    summaries[[i]] <- cbind(spec, data.frame(nPosteriorDraws = length(value),
      posteriorMean = mean(value), posteriorSD = sd(value),
      firstHalfMean = mean(head(value, half)), lastHalfMean = mean(tail(value, half)),
      meanAt75Percent = runningMean[quarter],
      lastQuarterMeanChange = tail(runningMean, 1) - runningMean[quarter]))
  }
  list(settings = list(nChains = 1L, nIter = nIter, burnIn = burnIn, thin = thin,
    scope = "Fixed scenarios: rho=-0.7, h2=0.3/gamma=0.4 and h2=0.7/gamma=0.2; phenotype and design replicate 1",
    mean = "Cumulative arithmetic mean of retained draws from each individual chain; original parameter scale",
    interpretation = "Visual stability of parameter means only; not proof of convergence or an audit of every fit"),
    summary = bind_rows(summaries), plotData = bind_rows(rows))
}

PlotMcmcRunningMeans <- function(result, fileName, analysis = c("GP", "BO")) {
  analysis <- match.arg(analysis)
  dat <- result$plotData
  crops <- c(Sb = "Sorghum", Gm = "Soybean")
  components <- c("DGE variance", "IGE variance", "SMA variance", "Residual variance",
    "DGE-IGE correlation")
  panelLevels <- unlist(lapply(components, function(x) paste(x, crops, sep = " | ")))
  dat$Panel <- factor(paste(dat$component, crops[dat$Crop], sep = " | "), levels = panelLevels)
  colors <- c(Mu = "#F8766D", PI = "#7CAE00", EI = "#00BFC4", UCB = "#C77CFF")
  cairo_pdf(fileName, width = 12.5, height = 12, family = "sans", onefile = TRUE)
  on.exit(dev.off())
  pageGroups <- if (analysis == "GP") c("Model", "trainingDesign") else c("method", "ParNo")
  pages <- unique(dat[pageGroups])
  for (i in seq_len(nrow(pages))) {
    keep <- rep(TRUE, nrow(dat))
    for (key in pageGroups) keep <- keep & dat[[key]] == pages[[key]][i]
    sub <- dat[keep, ]
    if (analysis == "GP") {
      sub$Series <- factor(paste0("h\u00b2 = ", sub$h2, ", \u03b3 = ", sub$gamma),
        levels = c("h\u00b2 = 0.3, \u03b3 = 0.4", "h\u00b2 = 0.7, \u03b3 = 0.2"))
      design <- c(All = "Design-D", Once = "Design-B")[pages$trainingDesign[i]]
      title <- paste0("GP prediction: ", sub("Model", "Model ", pages$Model[i]), " / ", design)
      subtitle <- "\u03c1 = -0.7; phenotype replicate 1; training-set replicate 1"
      p <- ggplot(sub, aes(Iteration, RunningMean, color = Series, linetype = Series)) +
        geom_line(linewidth = 1.15) +
        scale_color_manual(values = c("#0072B2", "#D55E00")) +
        scale_linetype_manual(values = c("solid", "longdash"))
      note <- if (any(sub$component == "DGE-IGE correlation"))
        "Correlation panels identify the genotype panel; variance panels identify the response crop.\n" else ""
    } else {
      sub$Series <- factor(paste0("Cycle ", sub$cycle), levels = paste0("Cycle ", c(1, 3, 5)))
      title <- paste0("BO: ", pages$method[i], " / Model 3")
      subtitle <- paste0("\u03c1 = -0.7; h\u00b2 = ", sub$h2[1], "; \u03b3 = ", sub$gamma[1],
        "; phenotype / initial-set replicate 1")
      p <- ggplot(sub, aes(Iteration, RunningMean, linetype = Series)) +
        geom_line(color = colors[pages$method[i]], linewidth = 1.2) +
        scale_linetype_manual(values = c("Cycle 1" = "solid", "Cycle 3" = "longdash", "Cycle 5" = "dotted"))
      note <- "Cycles use different training data and have different posterior means; compare stability within each line.\n"
    }
    print(p + facet_wrap(~Panel, ncol = 2, scales = "free_y") +
      scale_x_continuous(breaks = c(2000, 4000, 6000, 8000, 10000), labels = scales::label_comma()) +
      scale_y_continuous(breaks = scales::breaks_pretty(n = 4)) +
      labs(title = title, subtitle = subtitle, x = "MCMC iteration", y = "Cumulative posterior mean",
        color = NULL, linetype = NULL,
        caption = paste0("Mean of draws after burn-in (2,000); thinning = 20; 400 retained draws per chain.\n",
          note, "Fixed illustrative fits; original parameter scales, free y axes. Mean stability alone does not establish convergence.")) +
      RevisionPlotTheme() + theme(strip.text = element_text(size = 16),
        legend.key.width = grid::unit(2, "cm")))
  }
  invisible(nrow(pages))
}

EvaluateBoRunningMeans <- function(parSetDf, boDir, boRep = 1L, cycles = c(1L, 3L, 5L),
    nIter = 10000L, burnIn = 2000L, thin = 20L) {
  rows <- list()
  skipped <- list()
  for (method in c("Mu", "PI", "EI", "UCB")) {
    for (parNo in McmcExampleScenarios(parSetDf)) {
      root <- paste0(boDir, "parNo_", parNo, "/", method, "/")
      selected <- readRDS(paste0(root, "BO_result/2.0.selectedPairNo_allIteration_BORepNo_", boRep, ".rds"))
      for (cycle in cycles) {
        ### The terminal batch has no following fit; do not invent or reuse its trace.
        if (cycle >= length(selected)) {
          skipped[[length(skipped) + 1L]] <- data.frame(method = method, ParNo = parNo, cycle = cycle)
          next
        }
        fit <- readRDS(paste0(root, "BGLR_result/2.0.fitting_result_BORepNo_", boRep,
          "_iterNo_", cycle, ".rds"))
        for (crop in c("Sb", "Gm")) {
          settings <- unlist(fit[[crop]]$mcmcSettings)
          stopifnot(all(settings[c("nIter", "burnIn", "thin")] == c(nIter, burnIn, thin)))
          prefix <- paste0(root, "BGLR_parameter/2.0.", crop, "_BORepNo_", boRep, "_iterNo_", cycle, "_")
          kernels <- if (crop == "Sb") c("Sb", "Gm", "GbyG") else c("Gm", "Sb", "GbyG")
          rows[[length(rows) + 1L]] <- data.frame(method = method, ParNo = parNo,
            h2 = parSetDf$h2[parNo], gamma = parSetDf$gamma[parNo], phenotypeRep = parSetDf$rep[parNo],
            boRep = boRep, cycle = cycle, Crop = crop,
            component = c("DGE variance", "IGE variance", "SMA variance", "Residual variance"),
            sourceFile = paste0(prefix, c(paste0("ETA_", kernels, "_varU.dat"), "varE.dat")),
            nColumn = 1L, traceColumn = 1L)
        }
      }
    }
  }
  result <- ReadMcmcRunningMeans(bind_rows(rows), nIter, burnIn, thin)
  result$settings$skippedAfterStopping <- bind_rows(skipped)
  result
}

###### Equal-budget paired advantage curves ######
EvaluateCycleAdvantage <- function(boResLst, parSetDf, nRepBo = 10L,
    nBootstrap = 10000L, bootstrapSeed = 20261101L) {
  perRun <- bind_rows(lapply(boResLst, function(dat) {
    long <- tidyr::pivot_longer(dat, starts_with("Iter_"), names_to = "cycle", values_to = "bestTrueValue") %>%
      mutate(cycle = as.integer(sub("Iter_", "", cycle)))
    mu <- long %>% filter(Index == "Mu") %>% select(ParNo, BoRepNo, cycle, muBest = bestTrueValue)
    long %>% filter(Index != "Mu") %>% inner_join(mu, by = c("ParNo", "BoRepNo", "cycle")) %>%
      mutate(Advantage = bestTrueValue - muBest) %>% rename(method = Index, boRep = BoRepNo) %>%
      left_join(data.frame(ParNo = seq_len(nrow(parSetDf)), parSetDf), by = "ParNo") %>%
      rename(phenotypeRep = rep)
  }))
  stopifnot(all(is.finite(perRun$Advantage)), max(abs(perRun$Advantage[perRun$cycle == 1L])) < 1e-10)
  blocks <- perRun %>% group_by(genCor, h2, gamma, ParNo, phenotypeRep, method, cycle) %>%
    summarise(Advantage = mean(Advantage), nInitialRep = n(), .groups = "drop")
  stopifnot(all(blocks$nInitialRep == nRepBo))
  grouped <- perRun %>% group_by(genCor, h2, gamma, method, cycle)
  groups <- group_split(grouped)
  summary <- group_keys(grouped)
  nPhenotype <- n_distinct(perRun$phenotypeRep)
  stopifnot(nPhenotype == 10L, nRepBo == 10L)
  values <- vapply(groups, function(dat) {
    dat <- arrange(dat, boRep, phenotypeRep)
    stopifnot(nrow(dat) == nPhenotype * nRepBo,
      n_distinct(dat$phenotypeRep, dat$boRep) == nPhenotype * nRepBo,
      identical(sort(unique(dat$phenotypeRep)), seq_len(nPhenotype)),
      identical(sort(unique(dat$boRep)), seq_len(nRepBo)))
    dat$Advantage
  }, numeric(nPhenotype * nRepBo))
  ### Initial-set IDs recur across phenotypes: bootstrap both crossed factors.
  ### The same weights preserve method pairing and complete cycle trajectories.
  hadSeed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (hadSeed) oldSeed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (hadSeed) assign(".Random.seed", oldSeed, envir = .GlobalEnv) else
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv))
  set.seed(bootstrapSeed)
  phenotypeWeights <- rmultinom(nBootstrap, nPhenotype, rep(1 / nPhenotype, nPhenotype))
  initialWeights <- rmultinom(nBootstrap, nRepBo, rep(1 / nRepBo, nRepBo))
  weights <- t(vapply(seq_len(nBootstrap), function(i)
    as.vector(outer(phenotypeWeights[, i], initialWeights[, i])) /
      (nPhenotype * nRepBo), numeric(nPhenotype * nRepBo)))
  stopifnot(max(abs(rowSums(weights) - 1)) < 1e-12)
  bootstrapMeans <- weights %*% values
  summary$Mean <- colMeans(values)
  summary$SE <- apply(bootstrapMeans, 2, sd)
  summary$CILower <- apply(bootstrapMeans, 2, quantile, probs = .025, names = FALSE)
  summary$CIUpper <- apply(bootstrapMeans, 2, quantile, probs = .975, names = FALSE)
  summary$nRuns <- nPhenotype * nRepBo
  summary$nPhenotypeRep <- nPhenotype
  summary$nInitialRep <- nRepBo
  summary$budgetPairCount <- summary$cycle * 56L
  list(settings = list(metric = "best true objective among evaluated pairs (method) minus the matched Mu value",
    objective = "Sum of the two standardized true genomic values; evaluation only",
    budget = "Cycle 1 is the common initial 56 pairs; each further cycle adds up to 56 pairs; equal field-evaluation budget, not equal CPU time",
    stopping = "Carry the attained optimum forward after stopping; retain all 100 runs per scenario at every cycle",
    nBootstrap = nBootstrap, bootstrapSeed = bootstrapSeed,
    uncertainty = "Pointwise 95% percentile intervals from a crossed-factor (pigeonhole) bootstrap: independently resample 10 phenotype IDs and 10 initial-set IDs, then average their 100 intersections; preserve matched method contrasts and trajectories. Approximate, potentially conservative, not multiplicity-adjusted",
    reference = "Owen and Eckles (2012), Bootstrapping data arrays of arbitrary order, doi:10.1214/12-AOAS547"),
    perRun = perRun, phenotypeMeans = blocks, summary = summary)
}

PlotCycleAdvantage <- function(result, fileName) {
  dat <- result$summary %>% mutate(
    scenario = paste0("rho == ", genCor, "*','~~gamma == ", gamma),
    method = factor(method, levels = c("PI", "EI", "UCB")))
  scenarioOrder <- expand.grid(gamma = c(0, .2, .4), genCor = c(-.7, -.1, .7))
  dat$scenario <- factor(dat$scenario, levels = with(scenarioOrder,
    paste0("rho == ", genCor, "*','~~gamma == ", gamma)))
  colors <- c(PI = "#7CAE00", EI = "#00BFC4", UCB = "#C77CFF")
  cairo_pdf(fileName, width = 11.5, height = 9.5, family = "sans", onefile = TRUE)
  on.exit(dev.off())
  for (h in c(.3, .7)) {
    p <- ggplot(dat[dat$h2 == h, ], aes(cycle, Mean, color = method, fill = method, linetype = method)) +
      geom_hline(yintercept = 0, color = "grey35", linewidth = .65) +
      geom_ribbon(aes(ymin = CILower, ymax = CIUpper), alpha = .12, color = NA, show.legend = FALSE) +
      geom_line(linewidth = 1.15) +
      geom_point(data = function(x) x[x$cycle %in% c(1, 3, 6, 9, 12, 18, 24, 30, 40), ],
        aes(shape = method), size = 2.7) +
      facet_wrap(~scenario, ncol = 3, scales = "free", labeller = label_parsed) +
      scale_color_manual(values = colors) + scale_fill_manual(values = colors) +
      scale_linetype_manual(values = c(PI = "dashed", EI = "dotted", UCB = "dotdash")) +
      scale_shape_manual(values = c(PI = 16, EI = 17, UCB = 15)) +
      scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
      scale_y_continuous(breaks = scales::breaks_pretty(n = 4)) +
      labs(title = bquote("Best-value advantage relative to Mu"~~(h^2 == .(h))),
        x = "Evaluation cycle",
        y = "Best true objective: method minus Mu", color = NULL, fill = NULL, linetype = NULL, shape = NULL,
        caption = paste0("100 paired runs (10 phenotypes x 10 initial sets); shading: pointwise 95% crossed-bootstrap CI.\n",
          "56 pairs per cycle; initial batch = cycle 1. Positive values favor PI/EI/UCB. Axes vary by cell.")) +
      RevisionPlotTheme() + theme(
        strip.text = element_blank(), strip.background = element_blank(),
        panel.spacing = grid::unit(1.2, "mm"),
        axis.text = element_text(size = 17), axis.title = element_text(size = 19),
        axis.text.x = element_text(margin = margin(t = 2)),
        axis.text.y = element_text(margin = margin(r = 2)),
        plot.title = element_text(size = 21, face = "bold", margin = margin(b = 8)),
        plot.caption = element_text(size = 11.5, hjust = 0, margin = margin(t = 4)),
        legend.text = element_text(size = 18), legend.key.width = grid::unit(1.2, "cm"),
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = grid::unit(1, "mm"),
        plot.margin = margin(t = 8, r = 8, b = 8, l = 8))
    ### Shared row/column headers with independently scaled panels.
    ### Keeping facet_wrap avoids forcing a common y range within each rho row.
    grid::grid.newpage()
    panelPlot <- ggplotGrob(p)
    panels <- panelPlot$layout[grepl("^panel-", panelPlot$layout$name), ]
    stopifnot(nrow(panels) == 9L)
    panelPlot <- gtable::gtable_add_rows(panelPlot, grid::unit(c(9, 8), "mm"),
      pos = min(panels$t) - 1L)
    panels <- panelPlot$layout[grepl("^panel-", panelPlot$layout$name), ]
    panelPlot <- gtable::gtable_add_cols(panelPlot, grid::unit(c(8, 10), "mm"),
      pos = max(panels$r))
    panels <- panelPlot$layout[grepl("^panel-", panelPlot$layout$name), ]
    rowPositions <- sort(unique(panels$t))
    colPositions <- sort(unique(panels$l))
    header <- function(label, rotation = 0, size = 18) {
      grid::grobTree(grid::rectGrob(gp = grid::gpar(fill = "grey86", col = "grey35", lwd = .8)),
        grid::textGrob(label, rot = rotation,
          gp = grid::gpar(fontsize = size, fontface = "bold", fontfamily = "sans")))
    }
    panelPlot <- gtable::gtable_add_grob(panelPlot,
      header("Variance ratio of SMA to (DGE + IGE) (\u03b3)", size = 19),
      t = min(rowPositions) - 2L, l = min(colPositions), r = max(panels$r),
      name = "gamma-title")
    for (j in seq_len(3L)) {
      panelPlot <- gtable::gtable_add_grob(panelPlot, header(c("0", "0.2", "0.4")[j]),
        t = min(rowPositions) - 1L, l = colPositions[j],
        r = max(panels$r[panels$l == colPositions[j]]), name = paste0("gamma-value-", j))
      panelPlot <- gtable::gtable_add_grob(panelPlot, header(c("-0.7", "-0.1", "0.7")[j], rotation = 270),
        t = rowPositions[j], b = max(panels$b[panels$t == rowPositions[j]]),
        l = max(panels$r) + 1L, name = paste0("rho-value-", j))
    }
    panelPlot <- gtable::gtable_add_grob(panelPlot,
      header("Genetic correlation between DGE and IGE (\u03c1)", rotation = 270, size = 19),
      t = min(rowPositions), b = max(panels$b), l = max(panels$r) + 2L, name = "rho-title")
    grid::grid.draw(panelPlot)
  }
}

###### Exact paired comparisons for the ten independent phenotype replicates ######
ExactPairedWilcoxon <- function(x, y) {
  stopifnot(length(x) == length(y), all(is.finite(x)), all(is.finite(y)))
  ### Remove numerical round-off before identifying zeros and tied absolute differences.
  difference <- round(x - y, 12)
  difference <- difference[difference != 0]
  if (!length(difference)) return(1)
  if (length(difference) > 20L) stop("Exact enumeration requires at most 20 paired blocks.")
  ranks <- rank(abs(difference), ties.method = "average")
  observed <- sum(ranks[difference > 0])
  ### Enumerate all equally likely signs under the symmetric paired-difference null.
  nullSum <- 0
  for (rankNow in ranks) nullSum <- c(nullSum, nullSum + rankNow)
  center <- sum(ranks) / 2
  mean(abs(nullSum - center) >= abs(observed - center) - 1e-12)
}

###### Summary tables for revision diagnostics ######
ReadPublicationSelection <- function(parSetDf, boDir, phenoDir, nRep = 10L, nRepBo = 10L,
    kappa = round(qnorm(.95), 2)) {
  ### Reconstruct Figures 4, 5 and S10 from the original saved fits.
  ### Reconstruct diagnostics from saved fits and selections; no model fitting or new observations.
  mechanismMethods <- c('Mu', 'PI', 'EI', 'UCB')
  mechanismRho <- c(-0.7, 0.7)
  mechanismH2 <- 0.7
  mechanismGamma <- 0.2
  mechanismExampleRep <- mechanismExampleBoRep <- 1L
  mechanismInitialRows <- mechanismExampleRows <- mechanismDiscoveryRows <- list()
  mechanismScenarios <- unique(parSetDf[c('genCor', 'h2', 'gamma')])

  ##### 3.8.1. Saved-fit diagnostics #####
  for (scenarioNo in seq_len(nrow(mechanismScenarios))) {
    scenarioNow <- mechanismScenarios[scenarioNo, ]
    mechanismFocus <- scenarioNow$h2 == mechanismH2 && scenarioNow$gamma == mechanismGamma &&
      scenarioNow$genCor %in% mechanismRho
    phenoNow <- readRDS(paste0(phenoDir, '0.3.allSampledRes_phenotype_genCor_',
      scenarioNow$genCor, '_h2_', scenarioNow$h2, '_gamma_', scenarioNow$gamma, '_nRep_', nRep, '.rds'))
    parNoNow <- which(parSetDf$genCor == scenarioNow$genCor &
      parSetDf$h2 == scenarioNow$h2 & parSetDf$gamma == scenarioNow$gamma)
    for (parNo in parNoNow) {
      phenotypeRep <- parSetDf$rep[parNo]
      phenoRepNow <- phenoNow[[phenotypeRep]]
      mechanismTruth <- as.numeric(scale(phenoRepNow$dA + phenoRepNow$sBtoA + phenoRepNow$vA) +
        scale(phenoRepNow$dB + phenoRepNow$sAtoB + phenoRepNow$vB))
      mechanismOptimum <- which.max(mechanismTruth)
      stopifnot(sum(mechanismTruth == max(mechanismTruth)) == 1L)
      for (boRep in seq_len(nRepBo)) {
        mechanismBase <- paste0(boDir, 'parNo_', parNo, '/')
        mechanismSelections <- setNames(lapply(mechanismMethods, function(methodNow) {
          readRDS(paste0(mechanismBase, methodNow,
            '/BO_result/2.0.selectedPairNo_allIteration_BORepNo_', boRep, '.rds'))
        }), mechanismMethods)

        if (mechanismFocus) {
          ### All four first selections are reconstructed from the same initial GP.
          initialFit <- readRDS(paste0(mechanismBase,
            'EI/BGLR_result/2.0.fitting_result_BORepNo_', boRep, '_iterNo_1.rds'))
          initialObserved <- attr(initialFit, 'fitContext')$trainingPairNo
          initialCandidates <- setdiff(seq_along(mechanismTruth), initialObserved)
          initialMu <- initialFit$Sb$mean_A_B_AB + initialFit$Gm$mean_A_B_AB
          initialSd <- sqrt(initialFit$Sb$sd_A_B_AB^2 + initialFit$Gm$sd_A_B_AB^2)
          initialBest <- max(initialFit$Sb$meanObserved_A_B_AB + initialFit$Gm$meanObserved_A_B_AB)
          initialAF <- cal.AF.BGLR(initialBest, initialCandidates, initialMu, initialSd, kappa)
          initialMuTop <- head(order(initialMu, decreasing = TRUE), 56L)
          for (methodNow in mechanismMethods) {
            initialPick <- head(order(initialAF[[methodNow]], decreasing = TRUE), 56L)
            selectedNow <- mechanismSelections[[methodNow]]
            stopifnot(identical(as.integer(initialObserved), as.integer(selectedNow[[1]])),
              identical(as.integer(initialCandidates[initialPick]), as.integer(selectedNow[[2]])))
            mechanismInitialRows[[length(mechanismInitialRows) + 1L]] <- data.frame(
              genCor = scenarioNow$genCor, phenotypeRep = phenotypeRep, boRep = boRep, method = methodNow,
              differentFromMu = 56L - length(intersect(initialPick, initialMuTop)))
            if (phenotypeRep == mechanismExampleRep && boRep == mechanismExampleBoRep) {
              mechanismExampleRows[[length(mechanismExampleRows) + 1L]] <- data.frame(
                genCor = scenarioNow$genCor, method = methodNow, pairNo = initialCandidates,
                postMean = initialMu, postSD = initialSd, incumbent = initialBest,
                selected = seq_along(initialCandidates) %in% initialPick,
                optimum = initialCandidates == mechanismOptimum)
            }
          }
        }

        ### At discovery, compare the acquisition rank with the mean rank under that very same fit.
        for (methodNow in mechanismMethods) {
          selectedNow <- mechanismSelections[[methodNow]]
          discoveryCycle <- length(selectedNow)
          stopifnot(discoveryCycle >= 2L, mechanismOptimum %in% selectedNow[[discoveryCycle]],
            !mechanismOptimum %in% unlist(selectedNow[seq_len(discoveryCycle - 1L)]))
          discoveryFit <- readRDS(paste0(mechanismBase, methodNow,
            '/BGLR_result/2.0.fitting_result_BORepNo_', boRep,
            '_iterNo_', discoveryCycle - 1L, '.rds'))
          discoveryObserved <- attr(discoveryFit, 'fitContext')$trainingPairNo
          stopifnot(identical(as.integer(discoveryObserved),
            as.integer(unlist(selectedNow[seq_len(discoveryCycle - 1L)]))))
          discoveryCandidates <- setdiff(seq_along(mechanismTruth), discoveryObserved)
          discoveryMu <- discoveryFit$Sb$mean_A_B_AB + discoveryFit$Gm$mean_A_B_AB
          discoverySd <- sqrt(discoveryFit$Sb$sd_A_B_AB^2 + discoveryFit$Gm$sd_A_B_AB^2)
          discoveryBest <- NA_real_
          if (methodNow %in% c('PI', 'EI')) {
            stopifnot(identical(attr(discoveryFit, 'fitContext')$incumbentDefinition,
              'max_observed_posterior_mean'))
            discoveryBest <- max(discoveryFit$Sb$meanObserved_A_B_AB + discoveryFit$Gm$meanObserved_A_B_AB)
          }
          discoveryAF <- cal.AF.BGLR(discoveryBest, discoveryCandidates, discoveryMu, discoverySd, kappa)
          discoveryOrder <- order(discoveryAF[[methodNow]], decreasing = TRUE)
          discoveryMuOrder <- order(discoveryMu, decreasing = TRUE)
          discoveryOptimumIndex <- match(mechanismOptimum, discoveryCandidates)
          stopifnot(identical(as.integer(discoveryCandidates[head(discoveryOrder, 56L)]),
            as.integer(selectedNow[[discoveryCycle]])))
          mechanismDiscoveryRows[[length(mechanismDiscoveryRows) + 1L]] <- data.frame(
            genCor = scenarioNow$genCor, h2 = scenarioNow$h2, gamma = scenarioNow$gamma,
            phenotypeRep = phenotypeRep, boRep = boRep, method = methodNow,
            discoveryCycle = discoveryCycle,
            meanRank = match(discoveryOptimumIndex, discoveryMuOrder),
            acquisitionRank = match(discoveryOptimumIndex, discoveryOrder),
            outsideMeanTop56 = !discoveryOptimumIndex %in% head(discoveryMuOrder, 56L))
        }
      }
    }
    cat('Publication figures: checked scenario ', scenarioNo, '/', nrow(mechanismScenarios), '\n', sep = '')
    flush.console()
  }
  mechanismInitialDf <- bind_rows(mechanismInitialRows)
  mechanismExampleDf <- bind_rows(mechanismExampleRows)
  mechanismDiscoveryDf <- bind_rows(mechanismDiscoveryRows)
  stopifnot(nrow(mechanismInitialDf) == 2L * nRep * nRepBo * 4L,
    nrow(mechanismDiscoveryDf) == 18L * nRep * nRepBo * 4L,
    all(mechanismDiscoveryDf$acquisitionRank <= 56L),
    all(mechanismDiscoveryDf$outsideMeanTop56 == (mechanismDiscoveryDf$meanRank > 56L)),
    !any(mechanismDiscoveryDf$outsideMeanTop56[mechanismDiscoveryDf$method == 'Mu']))
    list(initial = mechanismInitialDf, example = mechanismExampleDf, discovery = mechanismDiscoveryDf)
}
ReadPublicationBoResults <- function(parSetDf, boDir, nRepBo = 10L) {
  methods <- c("Mu", "PI", "EI", "UCB")
  curves <- counts <- timings <- list()
  scenarios <- unique(parSetDf[c("genCor", "h2", "gamma")])
  for (scenarioNo in seq_len(nrow(scenarios))) {
    key <- scenarios[scenarioNo, ]
    parNos <- which(parSetDf$genCor == key$genCor & parSetDf$h2 == key$h2 & parSetDf$gamma == key$gamma)
    trajectories <- keys <- list()
    for (parNo in parNos) for (boRep in seq_len(nRepBo)) {
      initial <- NULL
      for (method in methods) {
        root <- paste0(boDir, "parNo_", parNo, "/", method, "/BO_result/2.0.")
        best <- readRDS(paste0(root, "g_best_allIteration_BORepNo_", boRep, ".rds"))
        selected <- readRDS(paste0(root, "selectedPairNo_allIteration_BORepNo_", boRep, ".rds"))
        if (is.null(initial)) initial <- selected[[1]]
        stopifnot(identical(initial, selected[[1]]), length(best) == length(selected),
          all(diff(best) >= -1e-10), isTRUE(attr(selected, "runTiming")$optimumFound))
        trajectories[[length(trajectories) + 1L]] <- best
        keys[[length(keys) + 1L]] <- data.frame(Index = method, BoRepNo = boRep, ParNo = parNo)
        counts[[length(counts) + 1L]] <- data.frame(nIter = length(best), Index = method,
          GenCor = key$genCor, h2 = key$h2, Gamma = key$gamma, BoRepNo = boRep, ParNo = parNo)
        timings[[length(timings) + 1L]] <- attr(selected, "runTiming")
      }
    }
    lastCycle <- max(lengths(trajectories))
    padded <- t(vapply(trajectories, function(x) c(x, rep(tail(x, 1), lastCycle - length(x))), numeric(lastCycle)))
    colnames(padded) <- paste0("Iter_", seq_len(lastCycle))
    curves[[scenarioNo]] <- cbind(as.data.frame(padded), bind_rows(keys))
  }
  attr(curves, "cycleConvention") <- "initial_batch_is_cycle_1"
  countDf <- bind_rows(counts)
  attr(countDf, "runTiming") <- bind_rows(timings)
  list(curves = curves, counts = countDf)
}
