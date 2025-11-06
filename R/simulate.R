# R/simulate.R
#
# Simulate target and source data for the transfer learning Cox model.

simulate_once_all <- function(n_target = 500, n_source = 2000,
                              beta2_source = 0.5, kappa_source = 1.5,
                              c1 = 4.55, covariate_shift = TRUE, L = 5) {
  if (!(L %in% c(5, 6))) stop("Only L = 5 or L = 6 is supported.")
  
  # --- Target data ---
  if (L == 5) {
    X1_t <- runif(n_target)
    X2_t <- rbinom(n_target, 1, 0.5)
    X3_t <- runif(n_target)
    X4_t <- runif(n_target)
    X5_t <- runif(n_target)
    x_t <- cbind(X1_t, X2_t, X3_t, X4_t, X5_t)
    beta_t <- c(-0.5, 0.5, 0.2, 0.1, 0.1)
  } else {
    Sigma <- outer(1:L, 1:L, function(i, j) 0.5^abs(i - j))
    x_t <- MASS::mvrnorm(n_target, mu = rep(0, L), Sigma = Sigma)
    beta_t <- c(rep(0.5, 5), rep(0, L - 5))
  }
  linpred_t <- x_t %*% beta_t
  scale_t <- sqrt(2 / 2) * exp(-linpred_t / 2)
  etime_t <- rweibull(n_target, shape = 2, scale = scale_t)
  ctime_t <- runif(n_target, 0, c1)
  stime_t <- pmin(etime_t, ctime_t)
  type_t <- as.numeric(etime_t <= ctime_t)
  
  target <- list(x = x_t, stime = stime_t, type = type_t)
  
  # --- Source data ---
  if (L == 5) {
    X1_s <- runif(n_source)
    X2_s <- rbinom(n_source, 1, 0.5)
    X3_s <- if (covariate_shift) rbeta(n_source, 1, 2) else runif(n_source)
    X4_s <- runif(n_source)
    X5_s <- runif(n_source)
    X_s <- cbind(X1_s, X2_s, X3_s, X4_s, X5_s)
    beta_s <- c(-0.5, beta2_source, 0.2, 0.1, 0.1)
  } else {
    Sigma <- outer(1:L, 1:L, function(i, j) 0.5^abs(i - j))
    if (covariate_shift) {
      eps <- matrix(rnorm(L, sd = 0.3))
      Sigma <- Sigma + eps %*% t(eps)
    }
    X_s <- MASS::mvrnorm(n_source, mu = rep(0, L), Sigma = Sigma)
    beta_s <- c(0.5, beta2_source, 0.5, 0.5, 0.5, rep(0, L - 5))
  }
  linpred_s <- X_s %*% beta_s
  scale_s <- sqrt(2 / kappa_source) * exp(-linpred_s / 2)
  etime_s <- rweibull(n_source, shape = 2, scale = scale_s)
  ctime_s <- runif(n_source, 0, c1)
  stime_s <- pmin(etime_s, ctime_s)
  type_s <- as.numeric(etime_s <= ctime_s)
  
  source <- list(list(x = X_s, stime = stime_s, type = type_s))
  
  return(list(target = target, source = source))
}


