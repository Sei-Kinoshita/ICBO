###### 1. Model comparisons ######
#' Create significance letters for model comparisons
#'
#' Average Fisher-z correlations over CV replicates (and crops for pooled plots)
#' within each phenotype replicate. Compare models using phenotype replicates
#' as paired blocks, with Friedman and Holm-adjusted exact paired Wilcoxon tests.
#' Inference is conditional on the training sets used in the simulation.
#'
#' @param dat A data frame containing `genCor`, `h2`, `gamma`, `CVMethod`,
#'   `rep`, `CV`, `Crop`, `Model`, and `Accuracy`.
#'
#' @return Significance letters and y-positions, with test results in the
#'   `modelComparison` attribute. Holm correction is applied within each panel.
MakeLettersDf <- function(dat) {
  groupLst <- dat %>%
    group_by(genCor, h2, gamma, CVMethod) %>%
    group_split()
  comparisonLst <- lapply(groupLst, function(oneDf) {
    oneDf$Model <- factor(oneDf$Model)
    if (any(!is.finite(oneDf$Accuracy)) || any(abs(oneDf$Accuracy) > 1)) {
      stop("GP model comparisons require finite correlations in [-1, 1].")
    }
    cellCount <- with(oneDf, table(factor(rep), factor(CV), factor(Crop), Model))
    if (any(cellCount != 1L)) {
      stop("GP model comparisons require matched phenotype, CV and crop observations for every model.")
    }
    ### The stored correlations are rounded; keep exact +/-1 finite for Fisher z.
    corBound <- 1 - sqrt(.Machine$double.eps)
    oneDf$FisherZ <- atanh(pmin(pmax(oneDf$Accuracy, -corBound), corBound))
    blockDf <- oneDf %>%
      group_by(rep, Model) %>%
      summarise(FisherZ = mean(FisherZ), .groups = "drop")
    scoreMat <- as.matrix(xtabs(FisherZ ~ rep + Model, data = blockDf))
    if (nrow(scoreMat) < 2L || ncol(scoreMat) < 2L) {
      stop("GP model comparisons require at least two phenotype replicates and two models.")
    }
    allTied <- all(apply(scoreMat, 1, function(x) all(x == x[1])))
    friedmanP <- if (allTied) 1 else friedman.test(scoreMat)$p.value
    modelPair <- combn(colnames(scoreMat), 2)
    pValue <- apply(modelPair, 2, function(pair) {
      x <- scoreMat[, pair[1]]
      y <- scoreMat[, pair[2]]
      if (all(x == y)) return(1)
      ExactPairedWilcoxon(x, y)
    })
    pAdjusted <- p.adjust(pValue, method = "holm")
    names(pAdjusted) <- paste(modelPair[1, ], modelPair[2, ], sep = "-")
    let <- multcompView::multcompLetters(pAdjusted)
    keyDf <- as.data.frame(oneDf[1, c("genCor", "h2", "gamma", "CVMethod")])
    cldTbl <- data.frame(Model = names(let$Letters), label = unname(let$Letters))
    ymaxTbl <- oneDf %>%
      group_by(Model) %>%
      summarise(y = max(Accuracy) + 0.02, .groups = "drop")
    lettersDf <- cbind(keyDf[rep(1, nrow(cldTbl)), ],
                       left_join(cldTbl, ymaxTbl, by = "Model"))
    testDf <- cbind(keyDf[rep(1, ncol(modelPair)), ],
      data.frame(Crop = paste(sort(unique(as.character(oneDf$Crop))), collapse = "+"),
        Model1 = modelPair[1, ], Model2 = modelPair[2, ],
        nPhenotypeRep = nrow(scoreMat), nTrainingRep = length(unique(oneDf$CV)),
        friedmanP = friedmanP, pValue = unname(pValue), pAdjusted = unname(pAdjusted),
        adjustment = "holm", test = "exact paired Wilcoxon",
        comparisonFamily = "10 model contrasts within scenario, CV design and crop panel",
        meanFisherZDifference = apply(modelPair, 2, function(pair) {
          mean(scoreMat[, pair[1]] - scoreMat[, pair[2]])
        }),
        block = "phenotype replicate; mean Fisher z over training replicates and included crops"))
    list(letters = lettersDf, tests = testDf)
  })
  lettersDf <- bind_rows(lapply(comparisonLst, function(x) x$letters))
  attr(lettersDf, "modelComparison") <- bind_rows(lapply(comparisonLst, function(x) x$tests))
  lettersDf
}

###### 2. Predictive-ability plots ######
#' Create a faceted boxplot of accuracy
#'
#' Generate a boxplot of `Accuracy` by `CVMethod`, grouped by `Model`,
#' faceted by `genCor` and `gamma`, and annotated with significance letters.
#'
#' @param dat A data frame containing `genCor`, `gamma`, `CVMethod`,
#'   `Model`, and `Accuracy`.
#' @param plotTitle A character string for the plot title. Default is `NULL`.
#'
#' @return A ggplot object.
MakePlot <- function(dat, plotTitle = NULL) {
  
  dat$CVMethod <- factor(dat$CVMethod, 
                         levels = c("All", "Once"), labels = c("CV1", "CV2"))
  dat$Model <- factor(dat$Model)
  dat$genCor <- factor(dat$genCor)
  dat$gamma <- factor(dat$gamma)
  
  lettersDf <- MakeLettersDf(dat)
  
  p <- ggplot(dat, aes(x = CVMethod, y = Accuracy, fill = Model)) +
    geom_boxplot(position = position_dodge(width = 0.8),
                 width = 0.7, outlier.shape = 16, outlier.size = 1.5) +
    geom_text(data = lettersDf, aes(x = CVMethod, y = y, label = label, group = Model),
      position = position_dodge(width = 0.8),
      vjust = 0, size = 5, show.legend = FALSE) +
    facet_grid(genCor ~ gamma) +
    scale_x_discrete(labels = c(CV1 = "Design-D", CV2 = "Design-B")) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = plotTitle, x = "Training design", y = "Predictive ability", fill = "Model",
         caption = "Exact paired Wilcoxon; 10 phenotype blocks; Holm correction within each scenario and training design.") +
    theme_bw(base_size = 14) +
    theme(strip.text = element_text(face = "bold", size = 14),
          axis.text = element_text(size = 17),
          axis.title = element_text(size = 14),
          plot.title = element_text(size = 15, face = "bold"),
          legend.position = "none")
  attr(p, "modelComparison") <- attr(lettersDf, "modelComparison")
  p
}

###### 3. MCMC running means ######
EvaluateGpRunningMeans <- function(parSetDf, modelDirs, dirGpPred, trainingRep = 1L,
    nIter = 10000L, burnIn = 2000L, thin = 20L) {
  rows <- list()
  for (modelNo in seq_along(modelDirs)) {
    for (design in c("All", "Once")) {
      root <- paste0(dirGpPred, modelDirs[modelNo], "/", design, "/CV", trainingRep, "/")
      for (parNo in McmcExampleScenarios(parSetDf)) {
        for (crop in if (modelNo <= 3L) c("Sb", "Gm") else "Multitrait") {
          prefix <- paste0(root, "modelPar/1.1.", crop, "_parNo_", parNo, "_", design, "_CV_", trainingRep, "_")
          fit <- readRDS(paste0(root, "allRes/1.1.", crop, "_fitting_allRes_parNo_", parNo,
            "_", design, "_CV_", trainingRep, ".rds"))
          addTrace <- function(suffix, nColumn, columns, crops, components, expected) {
            rows[[length(rows) + 1L]] <<- data.frame(trainingDesign = design, Model = paste0("Model", modelNo),
              ParNo = parNo, h2 = parSetDf$h2[parNo], gamma = parSetDf$gamma[parNo],
              phenotypeRep = parSetDf$rep[parNo], trainingRep = trainingRep,
              Crop = crops, component = components, sourceFile = paste0(prefix, suffix),
              nColumn = nColumn, traceColumn = columns, expectedMean = expected)
          }
          if (modelNo <= 3L) {
            stopifnot(fit$nIter == nIter, fit$burnIn == burnIn, fit$thin == thin)
            addTrace("varE.dat", 1L, 1L, crop, "Residual variance", fit$varE)
            if (modelNo == 1L) {
              addTrace("ETA_GenoRKHS_varU.dat", 1L, 1L, crop, "DGE variance", fit$ETA$GenoRKHS$varU)
            } else {
              addTrace("ETA_GenoSb_varU.dat", 1L, 1L, crop,
                if (crop == "Sb") "DGE variance" else "IGE variance", fit$ETA$GenoSb$varU)
              addTrace("ETA_GenoGm_varU.dat", 1L, 1L, crop,
                if (crop == "Gm") "DGE variance" else "IGE variance", fit$ETA$GenoGm$varU)
              if (modelNo == 3L) addTrace("ETA_GenoGbyG_varU.dat", 1L, 1L, crop,
                "SMA variance", fit$ETA$GenoGbyG$varU)
            }
          } else {
            addTrace("R.dat", 3L, c(1L, 3L), c("Sb", "Gm"), "Residual variance", diag(fit$resCov$R))
            for (kernelNo in 1:2) {
              kernel <- c("GenoSb", "GenoGm")[kernelNo]
              components <- if (kernelNo == 1L) c("DGE variance", "IGE variance") else
                c("IGE variance", "DGE variance")
              addTrace(paste0("Omega_", kernelNo, ".dat"), 3L, c(1L, 3L, 0L),
                c("Sb", "Gm", c("Sb", "Gm")[kernelNo]), c(components, "DGE-IGE correlation"),
                c(diag(fit$ETA[[kernel]]$Cov$Omega), NA_real_))
            }
            if (modelNo == 5L) addTrace("Omega_3.dat", 3L, c(1L, 3L), c("Sb", "Gm"),
              "SMA variance", diag(fit$ETA$GenoGbyG$Cov$Omega))
          }
        }
      }
    }
  }
  ReadMcmcRunningMeans(bind_rows(rows), nIter, burnIn, thin)
}

#' Check simulated genotype-level DGE-IGE correlations and multivariate covariance recovery.
#' Read existing covariance traces only; do not refit models or retain MCMC draws.
###### 4. Genetic-correlation recovery ######
EvaluateGeneticCorrelationRecovery <- function(parSetDf, pairNameDf, dirPheno,
    dirGpPred, saveDir, scriptID, cvCombLst,
    nIter = 10000L, burnIn = 2000L, thin = 20L) {
  stopifnot(nIter > burnIn, thin > 0, nIter %% thin == 0)
  cropNames <- c(Sb = "Sorghum", Gm = "Soybean")
  modelDirs <- c(Model4 = "Model4_Multi_DGE_IGE", Model5 = "Model5_Multi_GMA_SMA")
  designNames <- c(All = "Design-D", Once = "Design-B")
  meanCorrelation <- function(x) tanh(mean(atanh(x)))
  retained <- seq_len(nIter / thin) * thin > burnIn

  ### Read MCMC settings without executing the prediction script.
  for (script in "scripts/1.1.GPPred.r") {
    expressions <- parse(script)
    for (parameter in c("nIter", "burnIn", "thin")) {
      assigned <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
        identical(x[[2]], as.name(parameter)), as.list(expressions))
      stopifnot(length(assigned) == 1L,
        identical(as.numeric(eval(assigned[[1]][[3]])), as.numeric(get(parameter))))
    }
  }

  ### E[sample covariance of centered genotype effects] = kernelFactor * Omega.
  ### The simulation uses fixed marker-effect realizations, so covariance/kernelFactor is
  ### a moment-matched finite-panel benchmark, not a known generating kernel parameter.
  kernelFactors <- sapply(c("Sb", "Gm"), function(crop) {
    path <- if (crop == "Sb") "data/genotype/0.3.Sb_GRM_AMat_nInd_79.rds" else
      "data/genotype/0.3.Gm_GRM_AMat_nInd_66.rds"
    grm <- readRDS(path)
    stopifnot(setequal(rownames(grm), pairNameDf[[crop]]))
    (sum(diag(grm)) - sum(grm) / nrow(grm)) / (nrow(grm) - 1)
  })
  stopifnot(all(kernelFactors > 0))
  scenarioDf <- unique(parSetDf[, c("genCor", "h2", "gamma")])
  realizedRows <- list()
  for (scenarioNo in seq_len(nrow(scenarioDf))) {
    key <- scenarioDf[scenarioNo, ]
    phenoFile <- paste0(dirPheno, "0.3.allSampledRes_phenotype_genCor_", key$genCor,
      "_h2_", key$h2, "_gamma_", key$gamma, "_nRep_", length(unique(parSetDf$rep)), ".rds")
    pheno <- readRDS(phenoFile)
    for (repNo in seq_along(pheno)) {
      p <- pheno[[repNo]]
      parNo <- which(parSetDf$genCor == key$genCor & parSetDf$h2 == key$h2 &
        parSetDf$gamma == key$gamma & parSetDf$rep == repNo)
      stopifnot(length(parNo) == 1L)
      scales <- c(Sb = sd(p$yA), Gm = sd(p$yB))
      stopifnot(max(abs(as.numeric(scale(p$yA)) - as.numeric(p$scaled_yA))) < 1e-10,
        max(abs(as.numeric(scale(p$yB)) - as.numeric(p$scaled_yB))) < 1e-10)
      for (crop in names(cropNames)) {
        genotype <- unique(pairNameDf[[crop]])
        first <- match(genotype, pairNameDf[[crop]])
        d <- if (crop == "Sb") p$dA else p$dB
        s <- if (crop == "Sb") p$sAtoB else p$sBtoA
        otherCrop <- if (crop == "Sb") "Gm" else "Sb"
        stopifnot(identical(as.numeric(d), as.numeric(d[first][match(pairNameDf[[crop]], genotype)])),
          identical(as.numeric(s), as.numeric(s[first][match(pairNameDf[[crop]], genotype)])))
        r <- cor(d[first], s[first])
        savedR <- p$diagnostics$realizedDgeIgeCor[p$diagnostics$Crop == crop]
        stopifnot(length(savedR) == 1L, abs(r - savedR) < 1e-10)
        u <- cbind(DGE = d[first] / scales[crop], IGE = s[first] / scales[otherCrop])
        empiricalCov <- cov(u)
        benchmark <- empiricalCov / kernelFactors[crop]
        realizedRows[[length(realizedRows) + 1L]] <- data.frame(
          ParNo = parNo, genCor = key$genCor, h2 = key$h2, gamma = key$gamma,
          phenotypeRep = repNo, Crop = crop, nGenotype = length(genotype),
          realizedCorrelation = r, deviationFromTarget = r - key$genCor,
          dgePhenotypeSD = unname(scales[crop]), igePhenotypeSD = unname(scales[otherCrop]),
          realizedDGEVariance = empiricalCov[1, 1], realizedIGEVariance = empiricalCov[2, 2],
          realizedDGE_IGECovariance = empiricalCov[1, 2], kernelFactor = unname(kernelFactors[crop]),
          benchmarkOmegaDGE = benchmark[1, 1], benchmarkOmegaIGE = benchmark[2, 2],
          benchmarkOmegaDGE_IGE = benchmark[1, 2],
          benchmarkOmegaCovarianceAtTarget = key$genCor * sqrt(benchmark[1, 1] * benchmark[2, 2]))
      }
    }
  }
  realizedDf <- bind_rows(realizedRows)
  stopifnot(nrow(realizedDf) == 2L * nrow(parSetDf), !anyDuplicated(realizedDf[, c("ParNo", "Crop")]))

  ### Each saved row is vech(Omega) = (trait1 variance, covariance, trait2 variance).
  ### GenoSb has DGE/IGE in trait order; GenoGm has IGE/DGE, which must be reversed.
  posteriorRows <- list()
  for (model in names(modelDirs)) {
    for (design in names(designNames)) {
      trainingReps <- seq_along(cvCombLst[[design]])
      for (trainingRep in trainingReps) {
        for (parNo in seq_len(nrow(parSetDf))) {
            root <- paste0(dirGpPred, modelDirs[model], "/", design, "/CV", trainingRep, "/")
            prefix <- paste0(root, "modelPar/1.1.Multitrait_parNo_", parNo, "_", design, "_CV_", trainingRep, "_")
            fitFile <- paste0(root, "allRes/1.1.Multitrait_fitting_allRes_parNo_", parNo, "_", design, "_CV_", trainingRep, ".rds")
            nTrainingPair <- nrow(unique(cvCombLst[[design]][[trainingRep]]))
          if (!file.exists(fitFile)) stop("Missing completed multivariate fit: ", fitFile)
          ### Validate burn-in, thinning, trait order and trace precision against one saved fit per design/model.
          checkFit <- parNo == 1L && trainingRep == trainingReps[1]
          if (checkFit) referenceFit <- readRDS(fitFile)
          for (cropNo in seq_along(cropNames)) {
            crop <- names(cropNames)[cropNo]
            path <- paste0(prefix, "Omega_", cropNo, ".dat")
            if (!file.exists(path)) stop("Missing covariance trace: ", path)
            values <- scan(path, what = double(), quiet = TRUE)
            stopifnot(length(values) == 3L * length(retained), all(is.finite(values)))
            draws <- matrix(values, ncol = 3L, byrow = TRUE)[retained, , drop = FALSE]
            stopifnot(all(draws[, 1] > 0), all(draws[, 3] > 0),
              all(draws[, 1] * draws[, 3] > draws[, 2]^2))
            if (checkFit) {
              omega <- referenceFit$ETA[[paste0("Geno", crop)]]$Cov$Omega
              stopifnot(max(abs(colMeans(draws) - omega[lower.tri(omega, diag = TRUE)])) < 1e-7)
            }
            if (crop == "Gm") draws <- draws[, c(3, 2, 1), drop = FALSE]
            colnames(draws) <- c("omegaDGE", "omegaDGE_IGE", "omegaIGE")
            rho <- draws[, 2] / sqrt(draws[, 1] * draws[, 3])
            means <- colMeans(draws)
            limits <- apply(draws, 2, quantile, probs = c(0.025, 0.975), names = FALSE)
            truth <- realizedDf[realizedDf$ParNo == parNo & realizedDf$Crop == crop, ]
            ci <- quantile(rho, c(0.025, 0.975), names = FALSE)
            posteriorRows[[length(posteriorRows) + 1L]] <- cbind(truth,
              data.frame(Model = model, trainingDesign = design,
                trainingRep = trainingRep,
                nTrainingPair = nTrainingPair, nPosteriorDraws = length(rho),
                posteriorMeanCorrelation = mean(rho), posteriorMedianCorrelation = median(rho),
                correlationCrILower = ci[1], correlationCrIUpper = ci[2],
                correlationFromMeanOmega = means[2] / sqrt(means[1] * means[3]),
                omegaDGE = means[1], omegaDGE_IGE = means[2], omegaIGE = means[3],
                omegaDGELower = limits[1, 1], omegaDGEUpper = limits[2, 1],
                omegaDGE_IGELower = limits[1, 2], omegaDGE_IGEUpper = limits[2, 2],
                omegaIGELower = limits[1, 3], omegaIGEUpper = limits[2, 3],
                deviationFromRealized = mean(rho) - truth$realizedCorrelation,
                posteriorDeviationFromTarget = mean(rho) - truth$genCor,
                realizedWithinCrI = truth$realizedCorrelation >= ci[1] && truth$realizedCorrelation <= ci[2],
                targetWithinCrI = truth$genCor >= ci[1] && truth$genCor <= ci[2], row.names = NULL))
          }
          if (checkFit) rm(referenceFit)
        }
        cat("Covariance recovery:", model, design, "training replicate", trainingRep, "complete\n")
      }
    }
  }
  posteriorDf <- bind_rows(posteriorRows)
  expectedFits <- nrow(parSetDf) * length(modelDirs) * sum(lengths(cvCombLst))
  stopifnot(nrow(posteriorDf) == 2L * expectedFits,
    !anyDuplicated(posteriorDf[, c("ParNo", "Crop", "Model", "trainingDesign", "trainingRep")]))

  ### Average correlations on the Fisher-z scale; do not pool posterior samples from different fits.
  phenotypeDf <- posteriorDf %>%
    group_by(genCor, h2, gamma, phenotypeRep, ParNo, Crop, Model, trainingDesign) %>%
    summarise(estimatedCorrelation = meanCorrelation(posteriorMeanCorrelation),
      realizedCorrelation = first(realizedCorrelation),
      trainingRepSD = if (n() > 1L) sd(posteriorMeanCorrelation) else NA_real_,
      nTrainingRep = n(), .groups = "drop")
  recoverySummary <- phenotypeDf %>%
    group_by(genCor, h2, gamma, Crop, Model, trainingDesign) %>%
    summarise(nPhenotypeRep = n(), nTrainingRepPerPhenotype = first(nTrainingRep),
      estimatedCorrelationSD = sd(estimatedCorrelation),
      realizedCorrelationSD = sd(realizedCorrelation),
      estimatedCorrelation = meanCorrelation(estimatedCorrelation),
      realizedCorrelation = meanCorrelation(realizedCorrelation),
      meanTrainingRepSD = mean(trainingRepSD), .groups = "drop")
  fitSummary <- posteriorDf %>% group_by(genCor, h2, gamma, Crop, Model, trainingDesign) %>%
    summarise(nFits = n(), nTrainingPair = first(nTrainingPair), nPosteriorDraws = first(nPosteriorDraws),
      kernelFactor = first(kernelFactor),
      meanDeviationFromTarget = mean(posteriorDeviationFromTarget),
      meanDeviationFromRealized = mean(deviationFromRealized),
      rmseToRealized = sqrt(mean(deviationFromRealized^2)),
      realizedCrIInclusionRate = mean(realizedWithinCrI), targetCrIInclusionRate = mean(targetWithinCrI),
      meanCorrelationCrIWidth = mean(correlationCrIUpper - correlationCrILower),
      across(c(omegaDGE, omegaDGE_IGE, omegaIGE), mean),
      across(starts_with("benchmarkOmega"), mean), .groups = "drop")
  recoverySummary <- left_join(recoverySummary, fitSummary,
    by = c("genCor", "h2", "gamma", "Crop", "Model", "trainingDesign"))
  recoverySummary$Design <- unname(designNames[recoverySummary$trainingDesign])
  realizedSummary <- realizedDf %>% group_by(genCor, h2, gamma, Crop) %>%
    summarise(nPhenotypeRep = n(),
      realizedCorrelationSD = sd(realizedCorrelation),
      minimumCorrelation = min(realizedCorrelation), maximumCorrelation = max(realizedCorrelation),
      realizedCorrelation = meanCorrelation(realizedCorrelation),
      meanDeviationFromTarget = mean(deviationFromTarget), .groups = "drop")
  stopifnot(all(recoverySummary$nPhenotypeRep == length(unique(parSetDf$rep))),
    all(realizedSummary$nPhenotypeRep == length(unique(parSetDf$rep))))

  list(settings = list(nIter = nIter, burnIn = burnIn, thin = thin,
      correlation = "Posterior mean of Omega12 / sqrt(Omega11 * Omega22); not correlation of posterior mean effects",
      traitMapping = c(GenoSb = "trait1=DGE, trait2=IGE", GenoGm = "trait1=IGE, trait2=DGE; reordered to DGE, IGE"),
      covarianceScale = "Standardized phenotype units; DGE scaled by own-crop phenotype SD, IGE by partner-crop phenotype SD",
      covarianceBenchmark = "Realized centered genotype covariance divided by (trace(G)-sum(G)/n)/(n-1); finite-panel moment matching, not a known generating Omega",
      targetCovarianceBenchmark = "Target rho times the geometric mean of the two realized benchmark variances",
      interval = "Equal-tailed 95% posterior credible interval per fit after burn-in; reference inclusion rates are descriptive, not a calibration claim",
      averaging = "Correlations averaged on Fisher-z scale across training sets, then phenotype replicates; SD on the correlation scale. Matrix elements and signed errors use arithmetic means",
      validation = "Saved posterior covariance matrices agree with retained trace means within 1e-7 for one fit per model/design; all traces checked for length, finiteness and positive definiteness",
      source = "https://github.com/gdlc/BGLR-R/blob/master/R/Multitrait.R"),
    realized = realizedDf, realizedSummary = realizedSummary, posterior = posteriorDf,
    phenotypeMeans = phenotypeDf, summary = recoverySummary)
}

###### 5. Genetic-correlation plots ######
PlotGeneticCorrelationRecovery <- function(result, saveDir, scriptID = "1.2.") {
  cropNames <- c(Sb = "Sorghum", Gm = "Soybean")
  designNames <- c(All = "Design-D", Once = "Design-B")
  realizedDf <- result$realized
  realizedSummary <- result$realizedSummary
  phenotypeDf <- result$phenotypeMeans
  ### All matrix elements use standardized phenotype units and the DGE, IGE order for both crops.
  fileName <- paste0(saveDir, scriptID, "FigureS3_RealizedCorrelation.pdf")
  pdf(fileName, width = 12, height = 9, useDingbats = FALSE)
  plotTheme <- theme_bw(base_size = 15) + theme(legend.position = "bottom",
    strip.text = element_text(size = 15), axis.text = element_text(size = 13),
    plot.title = element_text(face = "bold"), plot.caption = element_text(size = 10, hjust = 0))
  for (crop in names(cropNames)) {
    dat <- realizedDf[realizedDf$Crop == crop, ]
    means <- realizedSummary[realizedSummary$Crop == crop, ]
    dat$displayTarget <- dat$genCor + (dat$phenotypeRep - (max(dat$phenotypeRep) + 1) / 2) * 0.007
    p <- ggplot(dat, aes(displayTarget, realizedCorrelation)) +
      geom_abline(slope = 1, intercept = 0, color = "grey45", linetype = "dashed", linewidth = 0.6) +
      geom_point(color = "#0072B2", size = 1.8, alpha = 0.7) +
      geom_errorbar(data = means, aes(x = genCor, ymin = realizedCorrelation - realizedCorrelationSD,
        ymax = realizedCorrelation + realizedCorrelationSD), width = 0.045, linewidth = 0.6) +
      geom_point(data = means, aes(x = genCor), shape = 18, size = 3) +
      facet_grid(h2 ~ gamma, labeller = label_both) +
      scale_x_continuous(breaks = c(-0.7, -0.1, 0.7), limits = c(-0.85, 0.85)) +
      scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
      labs(title = paste0(cropNames[crop], ": simulated DGE-IGE correlation"),
        x = expression("Target correlation (" * rho * ")"), y = "Realized genotype-level correlation",
        caption = "Blue points: 10 phenotype replicates per scenario (small horizontal offsets for visibility).\nBlack diamonds: Fisher-z averaged correlations; bars: +/- one SD across replicates. Dashed line: target = realized.\nDGE and IGE refer to the same genotype: its direct effect and its effect on the partner crop.") + plotTheme
    print(p)
  }
  dev.off()
  for (design in c("All", "Once")) {
    figure <- if (design == "All") "FigureS5" else "FigureS6"
    pdf(paste0(saveDir, scriptID, figure, "_CorrelationRecovery.pdf"), width = 12, height = 9, useDingbats = FALSE)
    for (crop in names(cropNames)) {
      dat <- phenotypeDf[phenotypeDf$trainingDesign == design & phenotypeDf$Crop == crop, ]
      p <- ggplot(dat, aes(realizedCorrelation, estimatedCorrelation, color = Model,
                           shape = factor(genCor))) +
        geom_abline(slope = 1, intercept = 0, color = "grey45", linetype = "dashed", linewidth = 0.6) +
        geom_point(size = 2.2, alpha = 0.85) +
        facet_grid(h2 ~ gamma, labeller = label_both) +
        scale_color_manual(values = c(Model4 = "#00BFC4", Model5 = "#C77CFF")) +
        scale_shape_manual(values = c(16, 17, 15)) +
        scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
        scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
        labs(title = paste0(cropNames[crop], ": covariance recovery / ", designNames[design]),
          subtitle = "520 observed pairs per fit",
          x = "Realized genotype-level correlation", y = "Estimated DGE-IGE correlation",
          color = "Model", shape = expression("Target " * rho),
          caption = paste0("Each point: one phenotype replicate; estimates are posterior means of correlation from retained covariance draws.\n",
            "Ten training-set estimates averaged on the Fisher-z scale within each phenotype replicate. ",
            "Dashed line: estimated = realized.\nPosterior 95% credible intervals and full 2 x 2 covariance summaries are retained in the existing accuracy RDS.")) + plotTheme
      print(p)
    }
    dev.off()
  }
}
