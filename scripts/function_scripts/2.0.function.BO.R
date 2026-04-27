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
###' @return posterior mean and variance for BO acquisition function
CalcPostMeanVar <- function(yObs, obs_idx, mu_hat,
                            K_A, K_B, K_AB,
                            sigma2_A, sigma2_B, sigma2_AB, sigma2_e) {
  
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
  
  list(mean_A = mean_A, sd_A = sqrt(var_A),
       mean_B = mean_B, sd_B = sqrt(var_B),
       mean_AB = mean_AB, sd_AB = sqrt(var_AB),
       mean_A_B = mean_A_B, sd_A_B = sqrt(var_A_B), ### u^A+u^B
       mean_A_B_AB = mean_A_B_AB, sd_A_B_AB  = sqrt(var_A_B_AB)) ### u^A+u^B+u^AB
}




##### 1.2. For Model 5 #####
###' @param YObs [matrix] nObsx2 observed pair phenotype（column for different crops）
###' @param obs_idx [vector] mx1 observed pair index
###' @param mu_hat  [vector] 2x1 estimated fixed effect, fit$mu from Multitrait()
###' @param K_A,K_B,K_AB [matrix] nPair x nPair kernels（K_A for Sb, K_B for Gm, K_AB for GbyG
###' @param G1,G2,G3 [matrix] 2x2 var-cov matrix between DGE & IGE, ETA[[r]] from Multitrait
###' @param R0 [matrix] 2x2 residual var-cov matrix resCov$R from Multitrait()
###' @return mean_u1, mean_u2, mean_u3（nCandx2 length vector）
###' @return cov_u1, cov_u2, cov_u3（2 x 2 x nCand array）
###' @return mean_uSum（nCand x 1 vector）
###' @return cov_uSum（nCand x 1 vector）
CalcPostMeanVar_MT <- function(yObs, obs_idx, mu_hat,
                               K_A, K_B, K_AB,
                               G1, G2, G3,
                               R0) {
  nPair <- nrow(K_A)
  nObs <- length(obs_idx)
  nCand <- nPair - nObs
  ### observed error
  rMat <- sweep(yObs, 2, mu_hat, "-")
  rVec <- as.vector(rMat)  ### 2*nObs vector
  
  ### Z^TGZ (=K) for each block of K
  K1_oo <- K_A[obs_idx, obs_idx, drop = FALSE]
  K2_oo <- K_B[obs_idx, obs_idx, drop = FALSE]
  K3_oo <- K_AB[obs_idx, obs_idx, drop = FALSE]
  
  K1_no <- K_A[-obs_idx, obs_idx, drop = FALSE]
  K2_no <- K_B[-obs_idx, obs_idx, drop = FALSE]
  K3_no <- K_AB[-obs_idx, obs_idx, drop = FALSE]
  
  K1_nn <- K_A[-obs_idx, -obs_idx, drop = FALSE]
  K2_nn <- K_B[-obs_idx, -obs_idx, drop = FALSE]
  K3_nn <- K_AB[-obs_idx, -obs_idx, drop = FALSE]
  
  ### observed var-cov matrix considering error
  Voo <- kronecker(K1_oo, G1) + kronecker(K2_oo, G2) + kronecker(K3_oo, G3) +
    kronecker(diag(nObs), R0)
  Voo_inv <- solve(Voo)
  
  ### Vno_u = Σ_r (K_r_no ⊗ G_r) s
  Vno_u <- kronecker(K1_no, G1) + kronecker(K2_no, G2) + kronecker(K3_no, G3)
  
  ### posterior mean of new u
  mean_vec_u <- Vno_u %*% (Voo_inv %*% rVec)
  mean_u <- matrix(mean_vec_u, nrow = nCand, ncol = 2)
  
  ### prior var-cov of u
  Vnn_u <- kronecker(K1_nn, G1) + kronecker(K2_nn, G2) + kronecker(K3_nn, G3)
  
  ### posterior var-cov of u
  Vpost_u <- Vnn_u - Vno_u %*% Voo_inv %*% t(Vno_u)
  
  Cov_u <- array(NA_real_, dim = c(2, 2, nCand))
  for (i in 1:nCand) {
    idx <- (2 * i - 1):(2 * i)
    Cov_u[, , i] <- Vpost_u[idx, idx, drop = FALSE]
  }
  
  ### uSum = u^A + u^B
  w <- c(1, 1)
  mean_uSum <- as.numeric(mean_u %*% w)
  var_uSum <- numeric(nCand)
  for (i in 1:nCand) {
    var_uSum[i] <- as.numeric(t(w) %*% Cov_u[, , i] %*% w)
  }
  
  ### posterior for each term
  Vno1 <- kronecker(K1_no, G1); Vnn1 <- kronecker(K1_nn, G1)
  Vno2 <- kronecker(K2_no, G2); Vnn2 <- kronecker(K2_nn, G2)
  Vno3 <- kronecker(K3_no, G3); Vnn3 <- kronecker(K3_nn, G3)
  
  mean1 <- matrix(Vno1 %*% (Voo_inv %*% rVec), nrow = nCand, ncol = 2)
  mean2 <- matrix(Vno2 %*% (Voo_inv %*% rVec), nrow = nCand, ncol = 2)
  mean3 <- matrix(Vno3 %*% (Voo_inv %*% rVec), nrow = nCand, ncol = 2)
  
  Vpost1 <- Vnn1 - Vno1 %*% Voo_inv %*% t(Vno1)
  Vpost2 <- Vnn2 - Vno2 %*% Voo_inv %*% t(Vno2)
  Vpost3 <- Vnn3 - Vno3 %*% Voo_inv %*% t(Vno3)
  
  cov1 <- array(NA_real_, dim = c(2, 2, nCand))
  cov2 <- array(NA_real_, dim = c(2, 2, nCand))
  cov3 <- array(NA_real_, dim = c(2, 2, nCand))
  for (i in 1:nCand) {
    idx <- (2 * i - 1):(2 * i)
    cov1[,,i] <- Vpost1[idx, idx, drop=FALSE]
    cov2[,,i] <- Vpost2[idx, idx, drop=FALSE]
    cov3[,,i] <- Vpost3[idx, idx, drop=FALSE]
  }
  
  list(
    mean_u1 = mean1, cov_u1 = cov1,
    mean_u2 = mean2, cov_u2 = cov2,
    mean_u3 = mean3, cov_u3 = cov3,
    mean_uSum = mean_uSum, cov_uSum = var_uSum
  )
}




###### 2. Function to Calculate acquisition function ######
###' @param candIdx [vector] name of candidate individual
###' @param gBest [scalar] max value of true genotypic value
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
