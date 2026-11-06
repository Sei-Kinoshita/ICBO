###### 1. Function to Calculate posterior mean & var ######
##### 1.1. For Model 3 #####
###' @param yObs [vector] nObsx1 observed pair phenotype
###' @param obs_idx [vector] mx1 observed pair index
###' @param mu_hat [scalar] estimated fixed effect, fit$mu from BGLR()
###' @param K_A [matrix] nPair x nPair matrix, Z^TGZ=K, for Sb
###' @param K_B [matrix] nPair x nPair matrix, Z^TGZ=K, for Gm
###' @param K_AB [matrix] nPair x nPair matrix, Z^TGZ=K, for GbyG
###' @param sigma2_A [scalar] estimated BLUP variance, fit$ETA$Sb$varU from BGLR()
###' @param sigma2_B [scalar] estimated BLUP variance, fit$ETA$Gm$varU from BGLR()
###' @param sigma2_AB [scalar] estimated BLUP variance, fit$ETA$GbyG$varU from BGLR()
###' @param sigma2_e [scalar] error variance, fit$varE from BGLR()
###' @param predictObserved [logical] also return observed pairs' latent genetic means
###' @return posterior mean and variance for BO acquisition function
CalcPostMeanVar <- function(yObs, obs_idx, mu_hat,
                            K_A, K_B, K_AB,
                            sigma2_A, sigma2_B, sigma2_AB, sigma2_e,
                            predictObserved = FALSE) {
  
  nObs <- length(obs_idx)
  ### error
  r <- yObs - mu_hat
  ### subset K matrix
  ### K matrix: observed × unobserved
  K_A_Obs <- K_A[obs_idx, obs_idx]
  K_B_Obs <- K_B[obs_idx, obs_idx]
  K_AB_Obs <- K_AB[obs_idx, obs_idx]
  ### K matrix: unobserved × observed
  K_A_New_Obs  <- K_A[-obs_idx, obs_idx]
  K_B_New_Obs  <- K_B[-obs_idx, obs_idx]
  K_AB_New_Obs <- K_AB[-obs_idx, obs_idx]
  ### K matrix: diagnal element
  K_A_diag  <- diag(K_A)[-obs_idx]
  K_B_diag  <- diag(K_B)[-obs_idx]
  K_AB_diag <- diag(K_AB)[-obs_idx]
  ### sum 3 kernels
  K_tot <- sigma2_A * K_A_Obs + sigma2_B * K_B_Obs + sigma2_AB * K_AB_Obs
  ### variance of observed y: Var(y) = K_tot + σ_e^2 I
  C <- K_tot + sigma2_e * diag(nObs)
  ### inverse matrix for the following analysis
  C_inv <- solve(C)
  
  ### posterior mean for unobserved pair (for each genetic effect)
  mean_A <- as.vector(sigma2_A * K_A_New_Obs %*% C_inv %*% r)
  mean_B <- as.vector(sigma2_B * K_B_New_Obs %*% C_inv %*% r)
  mean_AB <- as.vector(sigma2_AB * K_AB_New_Obs %*% C_inv %*% r)
  mean_A_B <- mean_A + mean_B
  mean_A_B_AB <- mean_A + mean_B + mean_AB
  
  ### posterior variance for unobserved pair (for each genetic effect)
  H_A <- K_A_New_Obs %*% C_inv
  H_B <- K_B_New_Obs %*% C_inv
  H_AB <- K_AB_New_Obs %*% C_inv
  
  var_A <- sigma2_A * K_A_diag - sigma2_A^2 * rowSums(H_A * K_A_New_Obs)
  var_B <- sigma2_B * K_B_diag - sigma2_B^2 * rowSums(H_B * K_B_New_Obs)
  var_AB <- sigma2_AB * K_AB_diag - sigma2_AB^2 * rowSums(H_AB * K_AB_New_Obs)
  ### posterior covariance
  cov_A_B <- sigma2_A * sigma2_B * rowSums(H_A * K_B_New_Obs)
  cov_A_AB <- sigma2_A * sigma2_AB * rowSums(H_A * K_AB_New_Obs)
  cov_B_AB <- sigma2_B * sigma2_AB * rowSums(H_B * K_AB_New_Obs)
  ### posterior variance for u^A+u^B
  var_A_B <- var_A + var_B - 2 * cov_A_B
  ### posterior variance for u^A+u^B+u^AB
  var_A_B_AB <- var_A + var_B + var_AB - 2 * (cov_A_B + cov_A_AB + cov_B_AB)
  
  postRes <- list(mean_A = mean_A, sd_A = sqrt(var_A),
       mean_B = mean_B, sd_B = sqrt(var_B),
       mean_AB = mean_AB, sd_AB = sqrt(var_AB),
       mean_A_B = mean_A_B, sd_A_B = sqrt(var_A_B), ### u^A+u^B
       mean_A_B_AB = mean_A_B_AB, sd_A_B_AB  = sqrt(var_A_B_AB)) ### u^A+u^B+u^AB
  if (predictObserved) {
    ### Same fit and centered latent genetic scale as the candidate means, in obs_idx order.
    postRes$meanObserved_A_B_AB <- as.vector(K_tot %*% (C_inv %*% r))
  }
  return(postRes)
}

###### 2. Function to Calculate acquisition function ######
###' @param candIdx [vector] name of candidate individual
###' @param gBest [scalar] maximum current posterior objective mean over observed pairs (EI/PI)
###' @param muCand [vector] nCand x 1 vector for candidate pairs' mean
###' @param sdCand [vector] nCand x 1 vector for candidate pairs' SD
###' @param kappa [scalar] parameter for UCB
cal.AF.BGLR <- function(gBest, candIdx, muCand, sdCand, kappa, ...) {
  
  ### Calculate acquisition function for candidate individuals
  z <- (muCand - gBest) / sdCand ### normalized distribution
  ### PI
  pi <- pnorm(q = z, mean = 0, sd = 1)
  ### EI
  ei <- (muCand - gBest) * pnorm(q = z, mean = 0, sd = 1) + 
    sdCand * dnorm(x = z, mean = 0, sd = 1)
  ### UCB
  ucb <- muCand + kappa * sdCand
  
  data.frame(ID = candIdx, 
             Mu = muCand, 
             PI = pi, 
             EI = ei, 
             UCB = ucb)
}

###### 3. Batch diagnostics (not used for selection) ######
BatchDiversity <- function(pairNo, pairDf, grm) {
  batchGrm <- grm[pairNo, pairNo, drop = FALSE]
  distance <- sqrt(pmax(outer(diag(batchGrm), diag(batchGrm), "+") -
                          2 * batchGrm, 0))
  diag(distance) <- Inf
  data.frame(nSorghumGenotype = length(unique(pairDf$Sb[pairNo])),
             nSoybeanGenotype = length(unique(pairDf$Gm[pairNo])),
             minimumPairwiseDistance = if (length(pairNo) > 1) min(distance) else NA_real_,
             meanNearestNeighborDistance = if (length(pairNo) > 1) {
               mean(apply(distance, 1, min))
             } else NA_real_)
}

###### 4. Diagnostics for the existing single chain ######
SingleChainDiagnostics <- function(saveAt, crop, burnIn, thin, essThreshold) {
  kernelName <- c("Sb", "Gm", "GbyG")
  componentName <- if (crop == "Sb") c("DGE", "IGE", "SMA") else c("IGE", "DGE", "SMA")
  do.call(rbind, lapply(seq_along(kernelName), function(i) {
    sampleFile <- paste0(saveAt, "ETA_", kernelName[i], "_varU.dat")
    draws <- scan(sampleFile, quiet = TRUE)
    draws <- draws[seq_along(draws) * thin > burnIn]
    ess <- if (length(draws) > 1 && all(is.finite(draws)) && var(draws) > 0) {
      as.numeric(coda::effectiveSize(coda::mcmc(draws)))
    } else NA_real_
    data.frame(Crop = crop, component = componentName[i], nChains = 1L,
               ESS = ess, essThreshold = essThreshold,
               convergencePass = is.finite(ess) && round(ess, 8) >= essThreshold)
  }))
}

###### 5. CV2-style balanced training pairs ######
SelectBroadCoverage <- function(pairDf, trainingSize, excludedPairNo = integer(0)) {
  sbName <- unique(as.character(pairDf$Sb))
  gmName <- unique(as.character(pairDf$Gm))
  stopifnot(trainingSize <= nrow(pairDf) - length(excludedPairNo))
  pairKey <- paste(pairDf$Sb, pairDf$Gm)
  for (attempt in 1:1000) {
    ### Reuse CV2's balanced marginal counts, while enforcing unique pairs.
    sbCount <- tabulate(rep(sample(seq_along(sbName)), length.out = trainingSize),
                        nbins = length(sbName))
    gmRemain <- tabulate(rep(sample(seq_along(gmName)), length.out = trainingSize),
                         nbins = length(gmName))
    sbOrder <- sample(seq_along(sbName))
    sbOrder <- sbOrder[order(sbCount[sbOrder], decreasing = TRUE)]
    selected <- integer(0)
    feasible <- TRUE
    for (sb in sbOrder) {
      pairNo <- match(paste(sbName[sb], gmName), pairKey)
      available <- which(gmRemain > 0 & !is.na(pairNo) & !(pairNo %in% excludedPairNo))
      available <- available[sample.int(length(available))]
      available <- available[order(gmRemain[available], decreasing = TRUE)]
      if (length(available) < sbCount[sb]) {
        feasible <- FALSE
        break
      }
      chosen <- head(available, sbCount[sb])
      selected <- c(selected, pairNo[chosen])
      gmRemain[chosen] <- gmRemain[chosen] - 1L
    }
    if (feasible && all(gmRemain == 0)) {
      return(selected[sample.int(length(selected))])
    }
  }
  stop("Could not construct unique training pairs with balanced marginal counts.")
}
