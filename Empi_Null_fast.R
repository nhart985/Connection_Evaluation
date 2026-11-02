library(cobs)
library(MASS)   # for rlm()
library(fmsb)
library(mvtnorm)
library(condMVNorm)

#─── Overall empirical‐null ────────────────────────────────────────────────────
# Calculate the negative log-likelihood assuming the Z-scores follow a homogeneous empirical null distribution
# arg = c(mu0, sd0)
negloglik_overall_fast <- function(arg, z, trunc_lower, trunc_upper, p0) {
  mu0  <- arg[1]
  sd0  <- arg[2]
  
  in_trunc <- (z >= trunc_lower & z <= trunc_upper)
  z0       <- z[in_trunc]
  N0       <- sum(in_trunc)
  N1       <- length(z) - N0
  
  # truncate‐probability (constant here)
  Q   <- pnorm((trunc_upper - mu0)/sd0) - pnorm((trunc_lower - mu0)/sd0)
  
  # sum of squared deviations
  ssq <- sum((z0 - mu0)^2)
  
  # negative log-likelihood
  - ( N0*log(p0)
      - N0*log(sd0)
      - ssq/(2 * sd0^2)
      + N1 * log(1 - p0 * Q)
  )
}

#─── Individualized empirical‐null ────────────────────────────────────────────
# Calculate the negative log-likelihood assuming the Z-scores variance is a function of the effective sample size
negloglik_indiv_fast <- function(phi, z, n_tilde, trunc_lower, trunc_upper, p0) {
  in_trunc <- (z >= trunc_lower & z <= trunc_upper)
  z0       <- z[in_trunc]; n0 <- n_tilde[in_trunc]
  z1       <- z[!in_trunc]; n1 <- n_tilde[!in_trunc]
  N0       <- length(z0)
  
  # variance‐vectors
  var0 <- 1 + n0 * phi
  var1 <- 1 + n1 * phi
  
  # truncated‐probabilities (vector if bounds vary; scalar otherwise)
  Q    <- pnorm(trunc_upper[!in_trunc] / sqrt(var1)) - pnorm(trunc_lower[!in_trunc] / sqrt(var1))
  
  ll0  <- - sum(
    log(var0)/2 +
      (z0^2) / (2 * var0)
  )
  ll1  <-  N0*log(p0) + sum(log(1 - p0 * Q))
  
  -(ll0 + ll1)
}

#─── Confounder‐adjusted empirical‐null ────────────────────────────────────────
# Calculate the negative log-likelihood assuming the Z-scores' mean is a function of the mean of covariates, and variance is a function of the effective sample size
# a function of the mean of covariates and the effective sample size
# arg = c(coefs, phi)
negloglik_confounder_fast <- function(arg, z, x_bar, n_tilde, trunc_lower, trunc_upper, p0) {
  # arg = c(coefs, phi)
  betas <- arg[-length(arg)]
  phi   <- arg[length(arg)]
  
  in_trunc <- (z >= trunc_lower & z <= trunc_upper)
  z0       <- z[in_trunc]
  n0       <- n_tilde[in_trunc]
  z1       <- z[!in_trunc]
  n1       <- n_tilde[!in_trunc]
  
  X0 <- as.matrix(x_bar)[in_trunc, , drop=FALSE]
  X1 <- as.matrix(x_bar)[!in_trunc, , drop=FALSE]
  
  N0 <- length(z0)
  
  var0 <- 1 + n0 * phi
  var1 <- 1 + n1 * phi
  
  # CORRECT: Matrix multiplication to get a VECTOR of individual means
  mu0_vec <- (X0 * sqrt(n0)) %*% betas
  mu1_vec <- (X1 * sqrt(n1)) %*% betas
  
  # CORRECT: Safely handle scalar vs vector bounds
  t_up <- if(length(trunc_upper) > 1) trunc_upper[!in_trunc] else trunc_upper
  t_low <- if(length(trunc_lower) > 1) trunc_lower[!in_trunc] else trunc_lower
  
  # truncated-probabilities
  upper1 <- (t_up - mu1_vec) / sqrt(var1)
  lower1 <- (t_low - mu1_vec) / sqrt(var1)
  Q      <- pnorm(upper1) - pnorm(lower1)
  
  ll0 <- - sum(
    log(var0)/2 +
      (z0 - mu0_vec)^2 / (2 * var0)
  )
  ll1 <-  N0*log(p0) + sum(log(1 - p0 * Q))
  
  return( -(ll0 + ll1) )
}

#─── Grid‐search version for confounder model ─────────────────────────────────
negloglik_confounder_gs_fast <- function(arg, z, x_bar, n_tilde, trunc_lower, trunc_upper, p0, phi) {
  # arg = coefficient vector; phi passed separately
  in_trunc <- (z >= trunc_lower & z <= trunc_upper)
  z0       <- z[in_trunc];   n0 <- n_tilde[in_trunc]
  z1       <- z[!in_trunc];  n1 <- n_tilde[!in_trunc]
  
  X0 <- as.matrix(x_bar)[in_trunc, , drop=FALSE]
  X1 <- as.matrix(x_bar)[!in_trunc, , drop=FALSE]
  
  N0 <- length(z0)
  
  var0 <- 1 + n0 * phi
  var1 <- 1 + n1 * phi
  
  mu0_vec <- (X0 * sqrt(n0)) %*% arg
  mu1_vec <- (X1 * sqrt(n1)) %*% arg
  # truncated‐probabilities (no mean shift in GS)
  Q <- pnorm((trunc_upper[!in_trunc] - mu1_vec) / sqrt(var1)) - pnorm((trunc_lower[!in_trunc] - mu1_vec)/ sqrt(var1))
  
  ll0 <- - sum(
    log(var0)/2 +
      (z0 - mu0_vec)^2 / (2 * var0)
  )
  ll1 <-  N0*log(p0) + sum(log(1 - p0 * Q))
  
  -(ll0 + ll1)
}

# ─── Confounder-adjust model with Poisson Z-scores ─────────────
negloglik_confounder_Pois <- function(arg, z, x_bar, E, E2, trunc_lower, trunc_upper, p0) {
  # arg = coefficient vector; phi passed separately
  in_trunc <- (z >= trunc_lower & z <= trunc_upper)
  z0       <- z[in_trunc];   n0 <- n[in_trunc]
  z1       <- z[!in_trunc];  n1 <- n[!in_trunc]
  
  P = ifelse(is.null(dim(x_bar)),1,dim(x_bar)[2])
  get_pdr=function(x_bar,arg) {
    return((as.matrix(x_bar)%*%arg[1:P])[,1])
  }
  
  temp=exp(get_pdr(x_bar,arg)+arg[P+1]/2+arg[P+2]/2)
  m=sqrt(E)*(temp-1)
  v=temp*(1+temp*(exp(arg[P+1])-1)*E+temp*exp(arg[P+1])*(exp(arg[P+2])-1)*E2/E)
  
  var0 <- v[in_trunc]
  var1 <- v[!in_trunc]
  
  mu0_vec <- m[in_trunc]
  
  # truncated‐probabilities (no mean shift in GS)
  Q <- pnorm(trunc_upper / sqrt(var1)) - pnorm(trunc_lower / sqrt(var1))
  
  ll0 <- - sum(
    log(var0)/2 +
      (z0 - mu0_vec)^2 / (2 * var0)
  )
  ll1 <-  N0*log(p0) + sum(log(1 - p0 * Q))
  
  -(ll0 + ll1)
}

#─── CRE empirical‐null ─────────────────────────────────────
negloglik_confounder_CRE_fast <- function(phi, z, n, trunc_lower, trunc_upper, p0) {
  in_trunc <- (z >= trunc_lower & z <= trunc_upper)
  z0       <- z[in_trunc];   n0 <- n[in_trunc]
  z1       <- z[!in_trunc];  n1 <- n[!in_trunc]
  N0       <- length(z0)
  
  var0 <- n0 * phi
  var1 <- n1 * phi
  
  Q <- pnorm(trunc_upper[!in_trunc] / sqrt(var1)) - pnorm(trunc_lower[!in_trunc] / sqrt(var1))
  
  ll0 <- - sum(
    log(var0)/2 +
      (z0^2) / (2 * var0)
  )
  ll1 <-  N0*log(p0) + sum(log(1 - p0 * Q))
  
  -(ll0 + ll1)
}

#' Overall Empirical Null
#' 
#' @param z a vector of Z-scores to be standardized
#' @param p_grid a grid of the assumed probability that a cluster is 'null',
#' with default from 0.8 to 0.999 with step 0.005
#' @param xlim a two-dimensional vector specify the center and the devaition of
#' the truncated interval, with default set to NULL

empirical_null_overall <- function(z, p_grid = seq(0.8,0.999,0.005), xlim = NULL) {
  
  N = length(z)
  
  if (is.null(xlim)) {
    b <- ifelse(N>500000, 1, 4.3 * exp(-0.26*log(N,10)))
    xlim <- c(median(z),b*IQR(z)/(2*qnorm(0.75)))
  }
  trunc_lower <- xlim[1]-xlim[2]
  trunc_upper <- xlim[1]+xlim[2]
  
  z0 <- z[which(z>=trunc_lower & z<=trunc_upper)]
  
  eval_p_grid <- sapply(p_grid, function(p){
    optim(par=c(mean(z0), sd(z0)), fn=negloglik_overall,
          z = z, trunc_lower = trunc_lower, 
          trunc_upper = trunc_upper, p0 = p,
          hessian = FALSE, method="Nelder-Mead")$value
  })
  
  p0_selected <- p_grid[which.min(eval_p_grid)]
  res_optim <- optim(par=c(mean(z0), sd(z0)), fn=negloglik_overall,
                     z = z, trunc_lower = trunc_lower, 
                     trunc_upper = trunc_upper, p0 = p0_selected,
                     hessian = FALSE, method="Nelder-Mead")
  return(list(est=c(res_optim$par, p0)))
}

#' Individualized Empirical Null account for cluster-level random quantity
#' 
#' @param z a vector of Z-scores to be standardized
#' @param p_grid a grid of the assumed probability that a cluster is 'null',
#' with default from 0.8 to 0.999 with step 0.005
#' @param phi_grid a grid of the parameter phi, with default from 
#' 0.01 to 5 with step 0.01
#' @param est_method estimation methods for the phi parameter, with options 
#' "optim" and "grid search". The "optim" use optim function, while the 
#' "grid search" just applied the grid search method on the given phi_grid.
#' The default is set to "optim"
#' @param cutoff a cutoff value to define the critical value of the truncated interval
#' @param sigma_est an estimate for the standard deviation to define the truncated interval
#' @param range a range parameter to define the truncated interval based on effective sample
#' size n. If range is missing and sigma_est is provided, the truncated interval will define
#' based on sigma_est.

empirical_null_ind_RQ <- function(z, n, p_grid=seq(0.8,0.999,0.005), 
                                  phi_grid=seq(0.01,5,0.01), # seq(0,5,0.1)
                                  est_method = "optim", # may change here to specific methods in optim
                                  cutoff=0.95, sigma_est, range) {
  N <- length(z)
  c_val <- qnorm((1-cutoff)/2, lower.tail=F)
  
  if(missing(sigma_est) & missing(range)) {
    rlm_est <- rlm(z~1)
    range=(rlm_est$s^2-1)/median(n)
    xlim <- c_val*sqrt(1+n*range)
  } else if(!missing(sigma_est) & missing(range)) {
    xlim <- c_val*sigma_est
  } else{
    xlim <- c_val*sqrt(n*range)
  }
  
  trunc_lower <- -xlim
  trunc_upper <- xlim
  
  eval_res <- sapply(p_grid, function(p){
    eval_phi_grid <- sapply(phi_grid, function(g){
      negloglik_indiv_fast(phi = g, z, n, trunc_lower, trunc_upper, p0 = p)
    })
    return(c(p = p, phi = phi_grid[which.min(eval_phi_grid)], neglik = min(eval_phi_grid)))
  })
  
  phi_selected <- eval_res["phi", which.min(eval_res["neglik",])]
  p_selected <- eval_res["p", which.min(eval_res["neglik",])]
  
  return(list(phi=sqrt(phi_selected), p = p_selected))
}

empirical_null_ind_CF_CRE <- function(z, n, p_grid = seq(0.8, 0.999, 0.005),
                                      phi_grid = seq(0.01, 5, 0.01),
                                      cutoff = 0.95, sigma_est, range) {
  N <- length(z)
  c_val <- qnorm(cutoff)
  
  if (missing(sigma_est) & missing(range)) {
    rlm_est <- MASS::rlm(z ~ 1)
    range <- (rlm_est$s^2 - 1) / median(n)
    xlim <- c_val * sqrt(n * range)
  } else if (!missing(sigma_est) & missing(range)) {
    xlim <- c_val * sigma_est
  } else {
    xlim <- c_val * sqrt(n * range)
  }
  
  trunc_lower <- -xlim
  trunc_upper <- xlim
  
  eval_res <- sapply(phi_grid, function(g) {
    eval_p_grid <- sapply(p_grid, function(p) {
      negloglik_confounder_CRE_fast(phi = g, z, n, trunc_lower, trunc_upper, p0 = p)
    })
    return(c(p = p_grid[which.min(eval_p_grid)], neglik = min(eval_p_grid)))
  })
  
  phi_selected <- phi_grid[which.min(eval_res["neglik", ])]
  p_selected <- eval_res["p", which.min(eval_res["neglik", ])]
  
  return(list(phi = sqrt(phi_selected), p = p_selected))
}

# empirical_null_ind_CF_CRE <- function(z, n, p_grid=seq(0.8,0.999,0.005), 
#                                       phi_grid=seq(0.01,5,0.01), est_method = "optim", # may change here to specific methods in optim
#                                       cutoff=0.95, sigma_est, range) {
#   N <- length(z)
#   c_val <- qnorm(cutoff)
#   
#   if(missing(sigma_est) & missing(range)) {
#     rlm_est <- rlm(z~1)
#     range=(rlm_est$s^2-1)/median(n)
#     xlim <- c_val*sqrt(n*range)
#   } else if(!missing(sigma_est) & missing(range)) {
#     xlim <- c_val*sigma_est
#   } else{
#     xlim <- c_val*sqrt(n*range)
#   }
#   
#   trunc_lower <- -xlim
#   trunc_upper <- xlim
#   
#   eval_p_grid <- sapply(p_grid, function(p){
#     optim(par=c(range), fn=negloglik_confounder_CRE_fast,
#           z = z, n = n, trunc_lower = trunc_lower, 
#           trunc_upper = trunc_upper, p0 = p,
#           hessian = FALSE, method="Nelder-Mead")$value
#   })
#   
#   # eval_res <- sapply(phi_grid, function(g){
#   #   eval_p_grid <- sapply(p_grid, function(p){
#   #     negloglik_confounder_CRE_fast(phi = g, z, n, trunc_lower, trunc_upper, p0 = p)
#   #   })
#   #   return(c(p = p_grid[which.min(eval_p_grid)], neglik = min(eval_p_grid)))
#   # })
#   
#   p_selected <- p_grid[which.min(eval_p_grid)]
#   phi_selected <- optim(par=c(range), fn=negloglik_confounder_CRE_fast,
#                 z = z, n = n, trunc_lower = trunc_lower, 
#                 trunc_upper = trunc_upper, p0 = p_selected,
#                 hessian = FALSE, method="Nelder-Mead")$par
#   
#   return(list(phi=phi_selected, p = p_selected))
# }



# empirical_null_ind_CF <- function(z, n, x_bar, p_grid=seq(0.8,0.999,0.005), 
#                                   phi_grid=seq(0.01,5,0.01), est_method = "optim", # may change here to specific methods in optim
#                                   cutoff=0.95, sigma_est, range) {
#   N <- length(z)
#   c_val <- qnorm(cutoff)
#   
#   rlm_est <- rlm(z~x_bar:sqrt(n)-1)
#   
#   if(missing(sigma_est) & missing(range)) {
#     range=(rlm_est$s^2-1)/median(n)
#     xlim <- c_val*sqrt(1+n*range)
#   } else if(!missing(sigma_est) & missing(range)) {
#     xlim <- c_val*sigma_est
#   } else{
#     xlim <- c_val*sqrt(1+n*range)
#   }
#   
#   trunc_lower <- -xlim + as.matrix(sqrt(n)*x_bar) %*% as.matrix(rlm_est$coefficients)
#   trunc_upper <- xlim + as.matrix(sqrt(n)*x_bar) %*% as.matrix(rlm_est$coefficients)
#   
#   # eval_res <- sapply(phi_grid, function(g){
#   #   eval_p_grid <- sapply(p_grid, function(p){
#   #     ifelse(is.null(dim(x_bar)), 
#   #            optimize(f=negloglik_confounder_gs_fast, interval = c(-10, 10),
#   #                     z = z, x_bar = x_bar, n_tilde = n, trunc_lower = trunc_lower, 
#   #                     trunc_upper = trunc_upper, phi = g, p0 = p)$objective,
#   #            optim(par=c(rlm_est$coefficients), fn=negloglik_confounder_gs_fast,
#   #                  z = z, x_bar = x_bar, n_tilde = n, trunc_lower = trunc_lower, 
#   #                  trunc_upper = trunc_upper, phi = g, p0 = p,
#   #                  hessian = FALSE, method="Nelder-Mead")$value)
#   #   })
#   #   return(c(p = p_grid[which.min(eval_p_grid)], neglik = min(eval_p_grid)))
#   # })
#   # 
#   # phi_selected <- phi_grid[which.min(eval_res["neglik",])]
#   # p_selected <- eval_res["p", which.min(eval_res["neglik",])]
#   # zeta <- ifelse(is.null(dim(x_bar)), 
#   #                optimize(f=negloglik_confounder_gs_fast, interval = c(-10, 10),
#   #                         z = z, x_bar = x_bar, n_tilde = n, trunc_lower = trunc_lower, 
#   #                         trunc_upper = trunc_upper, phi = phi_selected, p0 = p_selected)$minimum,
#   #                optim(par=c(rlm_est$coefficients), fn=negloglik_confounder_gs_fast,
#   #                      z = z, x_bar = x_bar, n_tilde = n, trunc_lower = trunc_lower, 
#   #                      trunc_upper = trunc_upper, phi = phi_selected, p0 = p_selected,
#   #                      hessian = FALSE, method="Nelder-Mead")$par)
#   
#   eval_p_grid <- sapply(p_grid, function(p){
#     optim(par=c(rlm_est$coefficients,range), fn=negloglik_confounder_fast,
#           z = z, x_bar = x_bar, n_tilde = n, trunc_lower = trunc_lower, 
#           trunc_upper = trunc_upper, p0 = p,
#           hessian = FALSE, method="Nelder-Mead")$value
#   })
#   
#   p_selected <- p_grid[which.min(eval_p_grid)]
#   pars <- optim(par=c(rlm_est$coefficients,range), fn=negloglik_confounder_fast,
#                 z = z, x_bar = x_bar, n_tilde = n, trunc_lower = trunc_lower, 
#                 trunc_upper = trunc_upper, p0 = p_selected,
#                 hessian = FALSE, method="Nelder-Mead")$par
#   return(list(zeta = pars[-length(pars)],
#               phi = pars[length(pars)],
#               p0 = p_selected))
# }

empirical_null_ind_CF <- function(z, n, x_bar, p_grid = seq(0.8, 0.999, 0.005), 
                                       phi_bounds = c(0.01, 5), # Replaced grid with bounds
                                       cutoff = 0.95, sigma_est, range) {
  
  N <- length(z)
  c_val <- qnorm(cutoff)
  
  rlm_est <- MASS::rlm(z ~ x_bar:sqrt(n) - 1)
  
  # 1. Fixed the missing 'range' calculation bug
  if (missing(sigma_est) & missing(range)) {
    range <- (rlm_est$s^2 - 1) / median(n)
    xlim <- c_val * sqrt(1 + n * range)
  } else if (!missing(sigma_est) & missing(range)) {
    xlim <- c_val * sigma_est
    range <- (rlm_est$s^2 - 1) / median(n) # Fixed!
  } else {
    xlim <- c_val * sqrt(1 + n * range)
  }
  
  # 2. Make sure the starting guess for phi is strictly within our bounds
  init_phi <- max(phi_bounds[1], min(phi_bounds[2], range))
  init_pars <- c(rlm_est$coefficients, init_phi)
  
  trunc_lower <- -xlim + as.matrix(sqrt(n) * x_bar) %*% as.matrix(rlm_est$coefficients)
  trunc_upper <- xlim + as.matrix(sqrt(n) * x_bar) %*% as.matrix(rlm_est$coefficients)
  
  # 3. Setup lower and upper bounds for L-BFGS-B
  # Coefficients (zeta) are unbounded (-Inf to Inf), but phi is strictly bounded.
  num_coefs <- length(rlm_est$coefficients)
  lower_bounds <- c(rep(-Inf, num_coefs), phi_bounds[1])
  upper_bounds <- c(rep(Inf, num_coefs), phi_bounds[2])
  
  # 4. Joint Optimization over the p_grid
  eval_p_grid <- sapply(p_grid, function(p) {
    optim(par = init_pars, 
          fn = negloglik_confounder_fast,
          z = z, x_bar = x_bar, n_tilde = n, 
          trunc_lower = trunc_lower, trunc_upper = trunc_upper, p0 = p,
          method = "L-BFGS-B",         # Changed to Bounded Optimization
          lower = lower_bounds,        # Applied lower bounds
          upper = upper_bounds,        # Applied upper bounds
          hessian = FALSE)$value
  })
  
  # 5. Extract the best p0 and re-run to get the final parameters
  p_selected <- p_grid[which.min(eval_p_grid)]
  
  best_fit <- optim(par = init_pars, 
                    fn = negloglik_confounder_fast,
                    z = z, x_bar = x_bar, n_tilde = n, 
                    trunc_lower = trunc_lower, trunc_upper = trunc_upper, p0 = p_selected,
                    method = "L-BFGS-B", 
                    lower = lower_bounds, 
                    upper = upper_bounds,
                    hessian = FALSE)
  
  pars <- best_fit$par
  
  return(list(zeta = pars[1:num_coefs],
              phi = pars[num_coefs + 1],
              p0 = p_selected))
}

empirical_null_ind_CF_pois <- function(z, E, E_2, x_bar, p_grid=seq(0.8,0.999,0.005), 
                                       phi_grid=seq(0.01,5,0.01), est_method = "optim", # may change here to specific methods in optim
                                       cutoff=0.95, sigma_est, range) {
  N <- length(z)
  c_val <- qnorm(cutoff)
  
  P = ifelse(is.null(dim(x_bar)),1,dim(x_bar)[2])
  get_pdr=function(x_bar,arg) {
    return((as.matrix(x_bar)%*%arg[1:P])[,1])
  }
  
  temp=exp(get_pdr(x_bar,arg)+arg[P+1]/2+arg[P+2]/2)
  m=sqrt(E)*(temp-1)
  v=temp*(1+temp*(exp(arg[P+1])-1)*E+temp*exp(arg[P+1])*(exp(arg[P+2])-1)*E2/E)
  
  rlm_est <- rlm(z~x_bar-1)
  
  if(missing(sigma_est) & missing(range)) {
    range=(rlm_est$s^2-1)/median(n)
    xlim <- c_val*sqrt(1+n*range)
  } else if(!missing(sigma_est) & missing(range)) {
    xlim <- c_val*sigma_est
  } else{
    xlim <- c_val*sqrt(1+n*range)
  }
  
  trunc_lower <- -xlim
  trunc_upper <- xlim
  
  eval_res <- sapply(phi_grid, function(g){
    eval_p_grid <- sapply(p_grid, function(p){
      ifelse(is.null(dim(x_bar)),
             optimize(f=negloglik_confounder_gs, interval = c(-10, 10),
                      z = z, x_bar = x_bar, n = n, trunc_lower = trunc_lower,
                      trunc_upper = trunc_upper, phi = g, p0 = p)$objective,
             optim(par=c(rlm_est$coefficients), fn=negloglik_confounder_gs,
                   z = z, x_bar = x_bar, n = n, trunc_lower = trunc_lower,
                   trunc_upper = trunc_upper, phi = g, p0 = p,
                   hessian = FALSE, method="Nelder-Mead")$value)
    })
    return(c(p = p_grid[which.min(eval_p_grid)], neglik = min(eval_p_grid)))
  })
  
  phi_selected <- phi_grid[which.min(eval_res["neglik",])]
  p_selected <- eval_res["p", which.min(eval_res["neglik",])]
  zeta <- ifelse(is.null(dim(x_bar)), 
                 optimize(f=negloglik_confounder_gs, interval = c(-10, 10),
                          z = z, x_bar = x_bar, n = n, trunc_lower = trunc_lower, 
                          trunc_upper = trunc_upper, phi = phi_selected, p0 = p_selected)$minimum,
                 optim(par=c(rlm_est$coefficients), fn=negloglik_confounder_gs,
                       z = z, x_bar = x_bar, n = n, trunc_lower = trunc_lower, 
                       trunc_upper = trunc_upper, phi = phi_selected, p0 = p_selected,
                       hessian = FALSE, method="Nelder-Mead")$par)
  return(list(zeta = zeta,
              phi = phi_selected,
              p0 = p_selected))
}


#-------------------------------------------------------
# Define the negative log‐likelihood for Poisson Z‐scores
#-------------------------------------------------------
negloglik_pois2 <- function(theta,
                            Z,        # vector of observed Z_i
                            Xbar,     # matrix (n_clusters × p) of cluster‑means
                            E,        # vector E_i = Σ_j exp(X_ij^T β)
                            E2,       # vector E_i^(2)=Σ_j exp(X_ij^T β)^2
                            trunc_l,  # scalar or vector lower trunc bounds
                            trunc_u,  # scalar or vector upper trunc bounds
                            p0) {     # probability of 'null' cluster
  p     <- ncol(Xbar)
  zeta     <- theta[1:p]
  sa2   <- theta[p + 1]   # σ_α^2
  se2   <- theta[p + 2]   # σ_ε^2
  
  # linear predictor
  eta     <- as.numeric(Xbar %*% zeta)
  
  # 1) mean
  shift   <- sa2/2 + se2/2
  mu      <- exp(eta + shift) * sqrt(E)
  
  # 2) variance
  T1      <- exp(eta + shift)
  T2      <- exp(2*eta + sa2 + se2)
  sigma2  <- T1 + T2 * (
    exp(sa2)*(exp(se2)-1)*(E2/E) +
      (exp(sa2)-1)*E
  )
  sigma   <- sqrt(sigma2)
  
  # 3) truncated‑normal mixture
  in_tr   <- (Z >= trunc_l & Z <= trunc_u)
  Q       <- pnorm((trunc_u - mu)/sigma) - pnorm((trunc_l - mu)/sigma)
  
  # null‑component log‑densities for those IN the truncation window
  ll0     <- sum(in_tr * ( log(p0)
                           - log(sigma)
                           - (Z - mu)^2 / (2*sigma2)
  ))
  
  # alt‑component for those OUTSIDE the window
  ll1     <- sum((!in_tr) * log(1 - p0*Q))
  
  # return negative log‑lik
  -(ll0 + ll1)
}

init_trunc_pois <- function(Z, Xbar, E, E2, cutoff = 0.95) {
  n      <- length(Z)
  c_val  <- qnorm((1 + cutoff)/2)
  
  fit0   <- rlm(log(Z / sqrt(E)) ~ Xbar - 1)
  zeta0  <- fit0$coefficients
  s0     <- fit0$s
  
  eta     <- as.numeric(Xbar %*% zeta0)
  shift   <- s0^2 / 2
  
  mu_hat  <- exp(eta + shift) * sqrt(E)
  T1      <- exp(eta + shift)
  T2      <- exp(2*eta + 2*shift)
  var_hat <- T1 + T2 * ((exp(s0^2)-1)*(E2/E) + (exp(s0^2)-1)*E)
  
  halfw   <- c_val * sqrt(var_hat)
  lower   <- mu_hat - halfw
  upper   <- mu_hat + halfw
  
  list(
    zeta_init   = zeta0,
    sigma2_init = s0^2,
    trunc_lower = lower,
    trunc_upper = upper
  )
}

estimate_pois <- function(Z, Xbar, E, E2,
                          p_grid = seq(0.8,0.99,0.005),
                          init     = NULL,
                          cutoff   = 0.95) {
  p <- ncol(Xbar)
  
  # 1) compute initial truncation & a rough sigma^2
  init_info <- init_trunc_pois(Z, Xbar, E, E2, cutoff)
  trunc_l   <- init_info$trunc_lower
  trunc_u   <- init_info$trunc_upper
  
  # 2) build a default 'init' vector if none supplied:
  #    [ zeta (p), sigmaAlpha2, sigmaEps2 ]
  if (is.null(init)) {
    # split the total variance equally as a starting guess
    half_sigma2 <- init_info$sigma2_init / 2
    init <- c(
      init_info$zeta_init,       # length p
      sigmaAlpha2 = half_sigma2,
      sigmaEps2   = half_sigma2
    )
  }
  
  # 3) grid‐search on p0, optimizing (zeta, sigmaAlpha2, sigmaEps2)
  results <- lapply(p_grid, function(p0_val) {
    fit <- optim(
      par     = init,
      fn      = negloglik_pois2,
      Z       = Z,
      Xbar    = Xbar,
      E       = E,
      E2      = E2,
      trunc_l = trunc_l,
      trunc_u = trunc_u,
      p0      = p0_val,
      method  = "L-BFGS-B",
      lower   = c(rep(-Inf,p), 1e-6, 1e-6),  # variances > 0
      upper   = c(rep( Inf,p), Inf,   Inf)
    )
    c(p0 = p0_val, negloglik = fit$value, fit$par)
  })
  
  # 4) assemble & pick best
  df <- do.call(rbind, results)
  colnames(df) <- c("p0", "negloglik",
                    sprintf("zeta%02d", 1:p),
                    "sigmaAlpha2", "sigmaEps2")
  best <- df[ which.min(df[,"negloglik"]), ]
  as.list(best)
}


negloglik_indiv_fast <- function(phi, z, n_tilde, trunc_lower, trunc_upper, p0) {
  
  in_trunc <- (z >= trunc_lower & z <= trunc_upper)
  z0       <- z[in_trunc]; n0 <- n_tilde[in_trunc]
  z1       <- z[!in_trunc]; n1 <- n_tilde[!in_trunc]
  N0       <- length(z0)

  # variance‐vectors
  var0 <- 1 + n0 * phi
  var1 <- 1 + n1 * phi
  
  # truncated‐probabilities (vector if bounds vary; scalar otherwise)
  Q    <- pnorm(trunc_upper[!in_trunc] / sqrt(var1)) - pnorm(trunc_lower[!in_trunc] / sqrt(var1))
  
  ll0  <- - sum(
    log(var0)/2 +
      (z0^2) / (2 * var0)
  )
  ll1  <-  N0*log(p0) + sum(log(1 - p0 * Q))
  
  -(ll0 + ll1)
}

empirical_null_ind_RQ_corr <- function(z, n, a_id, g_id, p_grid=seq(0.8,0.999,0.005),
                                  est_method = "optim", # may change here to specific methods in optim
                                  cutoff=0.95, sigma_est, range) {
  N <- length(z)
  c_val <- qnorm((1-cutoff)/2, lower.tail=F)
  
  if(missing(sigma_est) & missing(range)) {
    rlm_est <- rlm(z~1)
    range=(rlm_est$s^2-1)/median(n)
    xlim <- c_val*sqrt(1+n*range)
  } else if(!missing(sigma_est) & missing(range)) {
    xlim <- c_val*sigma_est
  } else{
    xlim <- c_val*sqrt(n*range)
  }
  
  trunc_lower <- -xlim
  trunc_upper <- xlim
  
  parms_out=matrix(NA,10,3)
  neg_log_L=vector("numeric")
  for(i in 1:10) {
    fit=optim(par=c(0,0,0),fn=negloglik_indiv_fast_corr_optim, lower=c(0,0,0), upper=c(0.1,0.1,0.1), method="L-BFGS-B",z=z, n=n, a_id=a_id, g_id=g_id, trunc_lower=trunc_lower, trunc_upper=trunc_upper, p0=1/i)
    parms_out[i,]=fit$par
    neg_log_L[i]=fit$value
  }
  
  parms_selected=parms_out[which.min(neg_log_L),]
  return(list(vars=parms_selected,phi=sqrt(sum(parms_selected))))
}

negloglik_indiv_fast_corr_optim <- function(par, z, n, a_id, g_id, trunc_lower, trunc_upper, p0) {
  return(negloglik_indiv_fast_corr(par[1],par[2],par[3], z, n, a_id, g_id, trunc_lower, trunc_upper, p0))
}

negloglik_indiv_fast_corr <- function(s2a, s2g, s2u, z, n_tilde, a_id, g_id, trunc_lower, trunc_upper, p0) {
  
  in_trunc <- (z >= trunc_lower & z <= trunc_upper)
  z0 <- z[in_trunc]; n0 <- n_tilde[in_trunc]
  z1 <- z[!in_trunc]; n1 <- n_tilde[!in_trunc]
  N <- length(z)
  N0 <- length(z0)
  N1 <- N-N0
  a_id0 <- a_id[in_trunc]
  g_id0 <- g_id[in_trunc]
  a_id1 <- a_id[!in_trunc]
  g_id1 <- g_id[!in_trunc]

  get_indices <- function(id) {
    dups <- (duplicated(id) | duplicated(id,fromLast=T))
    ind_split <- split(which(dups), id[dups])
    f <- function(x) {
      return(t(combn(x,m=2)))
    }
    ind_pairs <- do.call(rbind,lapply(ind_split,f))
    ind_pairs <- rbind(ind_pairs,ind_pairs[,c(2,1)])
  }
  
  # variance‐vectors
  var0 <- 1 + n0*(s2a+s2g+s2u)
  var1 <- 1 + n1*(s2a+s2g+s2u)
  sigma0 <- diag(var0)
  sigma1 <- diag(var1)
  a_ind0 <- get_indices(a_id0)
  g_ind0 <- get_indices(g_id0)
  sigma0[a_ind0] <- s2a*sqrt(n0[a_ind0[,1]]*n0[a_ind0[,2]])
  sigma0[g_ind0] <- s2g*sqrt(n0[g_ind0[,1]]*n0[g_ind0[,2]])
  a_ind1 <- get_indices(a_id1)
  g_ind1 <- get_indices(g_id1)
  sigma1[a_ind1] <- s2a*sqrt(n1[a_ind1[,1]]*n1[a_ind1[,2]])
  sigma1[g_ind1] <- s2g*sqrt(n1[g_ind1[,1]]*n1[g_ind1[,2]])

  mu <- rep(0,length(z1))
  Q_first  <- pnorm(trunc_upper[!in_trunc][1]/sqrt(var1[1])) - pnorm(trunc_lower[!in_trunc][1]/sqrt(var1[1]))
  ll1 <- log(1-p0*Q_first)
  for(i in 2:N1) {
    cond_pars <- condMVN(mean=mu,sigma=sigma1,dependent.ind=i,given.ind=1:(i-1),X.given=z1[1:(i-1)])
    Q_i  <- pnorm(trunc_upper[!in_trunc][i],cond_pars$condMean,sqrt(cond_pars$condVar))-pnorm(trunc_lower[!in_trunc][i],cond_pars$condMean,sqrt(cond_pars$condVar))
    ll1 <- ll1+log(1-p0*Q_i)
  }
  ll0  <- N0*log(p0)+dmvnorm(z0,mean=rep(0,length(z0)),sigma=sigma0,log=T)
  
  out <- -(ll0 + ll1)
  if(is.infinite(out) & out > 0) {
    out <- .Machine$double.xmax
  }
  return(out)
}



