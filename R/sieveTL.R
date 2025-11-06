# R/sieveTL.R
#
# Transfer learning Cox model with spline basis.
# Provides:
#   - transfer_loss_fast / transfer_grad_fast
#   - debias_loss_fast / debias_grad_fast
#   - sievetl (estimates beta, gamma)
#   - make_baseline_from_gamma (baseline hazard and cumulative baseline hazard)

library(pracma)
library(splines)

# ---------------- Transfer loss ----------------

transfer_loss_fast <- function(theta, p, L,
                               all.x, all.stime, all.type,
                               G_mat, time_grid, g_time_mat) {
  # theta = [gamma(1:p), beta(1:L)]
  gamma <- theta[1:p]
  beta  <- theta[(p + 1):(L + p)]
  
  xb <- as.vector(all.x %*% beta)
  if (p == 1) {
    g_gamma <- as.vector(G_mat * gamma)
  } else {
    g_gamma <- as.vector(G_mat %*% gamma)
  }
  term1 <- all.type * (xb + g_gamma)
  
  exp_gamma_g <- exp(g_time_mat %*% gamma)
  int_vec <- vapply(all.stime, function(Ti) {
    idx <- which(time_grid <= Ti)
    if (length(idx) == 0) return(0)
    trapz(time_grid[idx], exp_gamma_g[idx])
  }, numeric(1))
  
  term2 <- exp(xb) * int_vec
  loss <- -mean(term1 - term2, na.rm = TRUE)
  return(loss)
}

# ---------------- Smoothed L1 helpers (pseudo-Huber) ----------------

.smooth_abs      <- function(u, eps) sqrt(u*u + eps*eps) - eps
.smooth_abs_grad <- function(u, eps) u / sqrt(u*u + eps*eps)

# ---------------- Transfer gradient (no penalty) ----------------

transfer_grad_fast <- function(theta, p, L,
                               all.x, all.stime, all.type,
                               G_mat, time_grid, g_time_mat) {
  gamma <- theta[1:p]
  beta  <- theta[(p + 1):(L + p)]
  
  xb      <- as.vector(all.x %*% beta)
  exp_xb  <- exp(xb)
  eg_time <- as.vector(g_time_mat %*% gamma)
  exp_eg  <- exp(eg_time)
  
  # integral of exp(gamma' g(t)) from 0 to T_i
  int_vec <- vapply(all.stime, function(Ti){
    idx <- which(time_grid <= Ti)
    if (length(idx) < 2) return(0)
    pracma::trapz(time_grid[idx], exp_eg[idx])
  }, numeric(1))
  
  # for each i, j: integral of g_j(t) * exp(gamma' g(t)) from 0 to T_i
  term_mat <- matrix(0, nrow = length(all.stime), ncol = p)
  for (i in seq_along(all.stime)) {
    idx <- which(time_grid <= all.stime[i])
    if (length(idx) < 2) next
    for (j in 1:p) {
      term_mat[i, j] <- pracma::trapz(
        time_grid[idx],
        g_time_mat[idx, j] * exp_eg[idx]
      )
    }
  }
  
  grad_gamma <- -colMeans(G_mat * all.type, na.rm = TRUE) +
    colMeans(exp_xb * term_mat, na.rm = TRUE)
  grad_beta  <- -colMeans(all.x * all.type,  na.rm = TRUE) +
    colMeans(all.x * (exp_xb * int_vec),     na.rm = TRUE)
  
  c(grad_gamma, grad_beta)
}

# ---------------- Debias loss (with smoothed L1) ----------------

debias_loss_fast <- function(delta, lambda_zeta, lambda_eta, L, p,
                             target, G_mat_target, time_grid_target,
                             g_time_mat_target, gamma_hat_A, beta_hat_A,
                             eps = 1e-4) {
  # delta = [zeta(1:p), eta(1:L)]
  zeta <- delta[1:p]
  eta  <- delta[(p + 1):(L + p)]
  
  gamma <- gamma_hat_A + zeta
  beta  <- beta_hat_A  + eta
  
  xb      <- as.vector(target$x %*% beta)
  g_gamma <- if (p == 1)
    as.vector(G_mat_target * gamma)
  else
    as.vector(G_mat_target %*% gamma)
  term1   <- target$type * (xb + g_gamma)
  
  eg_time <- as.vector(g_time_mat_target %*% gamma)
  exp_eg  <- exp(eg_time)
  
  int_vec <- vapply(target$stime, function(Ti){
    idx <- which(time_grid_target <= Ti)
    if (length(idx) < 2) return(0)
    pracma::trapz(time_grid_target[idx], exp_eg[idx])
  }, numeric(1))
  
  term2 <- exp(xb) * int_vec
  nll   <- -mean(term1 - term2, na.rm = TRUE)
  
  # smoothed L1 penalties
  pen <- lambda_zeta * sum(.smooth_abs(zeta, eps)) +
    lambda_eta * sum(.smooth_abs(eta, eps))
  
  nll + pen
}

# ---------------- Debias gradient (with smoothed L1) ----------------

debias_grad_fast <- function(delta, lambda_zeta, lambda_eta, L, p,
                             target, G_mat_target, time_grid_target,
                             g_time_mat_target, gamma_hat_A, beta_hat_A,
                             eps = 1e-4) {
  zeta <- delta[1:p]
  eta  <- delta[(p + 1):(L + p)]
  
  gamma <- gamma_hat_A + zeta
  beta  <- beta_hat_A  + eta
  
  xb      <- as.vector(target$x %*% beta)
  exp_xb  <- exp(xb)
  eg_time <- as.vector(g_time_mat_target %*% gamma)
  exp_eg  <- exp(eg_time)
  
  int_vec <- vapply(target$stime, function(Ti){
    idx <- which(time_grid_target <= Ti)
    if (length(idx) < 2) return(0)
    pracma::trapz(time_grid_target[idx], exp_eg[idx])
  }, numeric(1))
  
  term_mat <- matrix(0, nrow = length(target$stime), ncol = p)
  for (i in seq_along(target$stime)) {
    idx <- which(time_grid_target <= target$stime[i])
    if (length(idx) < 2) next
    for (j in 1:p) {
      term_mat[i, j] <- pracma::trapz(
        time_grid_target[idx],
        g_time_mat_target[idx, j] * exp_eg[idx]
      )
    }
  }
  
  grad_zeta <- -colMeans(G_mat_target * target$type, na.rm = TRUE) +
    colMeans(exp_xb * term_mat, na.rm = TRUE)
  grad_eta  <- -colMeans(target$x * target$type,  na.rm = TRUE) +
    colMeans(target$x * (exp_xb * int_vec),      na.rm = TRUE)
  
  grad_zeta <- grad_zeta + lambda_zeta * .smooth_abs_grad(zeta, eps)
  grad_eta  <- grad_eta  + lambda_eta  * .smooth_abs_grad(eta,  eps)
  
  c(grad_zeta, grad_eta)
}

# ---------------- Main wrapper: sievetl (spline only) ----------------

sievetl <- function(target, source = NULL,
                    p = 5, L = 5,
                    lambda_zeta = 0.05, lambda_eta = 0.05) {
  
  transfer.source.id <- if (is.null(source)) integer(0) else seq_along(source)
  
  # stack data (target + sources)
  all.x <- as.matrix(do.call(
    rbind,
    lapply(c(0, transfer.source.id), function(k){
      if (k == 0) target$x else source[[k]]$x
    })
  ))
  all.stime <- unlist(lapply(c(0, transfer.source.id), function(k){
    if (k == 0) target$stime else source[[k]]$stime
  }))
  all.type <- unlist(lapply(c(0, transfer.source.id), function(k){
    if (k == 0) target$type else source[[k]]$type
  }))
  all_time <- all.stime
  
  # spline basis constructed from combined times
  B0 <- splines::bs(all_time, df = p, degree = 3, intercept = TRUE)
  at    <- attributes(B0)
  knots <- at$knots
  bkn   <- at$Boundary.knots
  deg   <- at$degree
  
  g_funcs <- lapply(seq_len(ncol(B0)), function(j){
    function(t){
      Bt <- splines::bs(t, degree = deg, knots = knots,
                        Boundary.knots = bkn, intercept = TRUE)
      as.numeric(Bt[, j])
    }
  })
  
  # G matrix over all samples
  if (p == 1) {
    G_mat <- matrix(sapply(all.stime, function(ti) g_funcs[[1]](ti)), ncol = 1)
  } else {
    G_mat <- t(sapply(all.stime,
                      function(ti) sapply(g_funcs, function(gj) gj(ti))))
  }
  
  # integration grid
  time_grid  <- seq(min(all_time), max(all_time), length.out = 500)
  g_time_mat <- if (p == 1) {
    matrix(g_funcs[[1]](time_grid), ncol = 1)
  } else {
    tmp <- sapply(g_funcs, function(gj) gj(time_grid))
    if (is.vector(tmp)) matrix(tmp, ncol = p) else tmp
  }
  
  # transfer step (unpenalized)
  theta0 <- rep(0, L + p)  # [gamma(1:p), beta(1:L)]
  loss_wrap <- function(th)
    transfer_loss_fast(th, p, L, all.x, all.stime, all.type,
                       G_mat, time_grid, g_time_mat)
  grad_wrap <- function(th)
    transfer_grad_fast(th, p, L, all.x, all.stime, all.type,
                       G_mat, time_grid, g_time_mat)
  
  fit <- nlminb(start = theta0,
                objective = loss_wrap,
                gradient  = grad_wrap,
                control   = list(iter.max = 5000, eval.max = 5000))
  
  theta_hat_A <- fit$par
  gamma_hat_A <- theta_hat_A[1:p]
  beta_hat_A  <- theta_hat_A[(p + 1):(L + p)]
  
  # debias step (smoothed L1)
  if (p == 1) {
    G_mat_target <- matrix(g_funcs[[1]](target$stime), ncol = 1)
  } else {
    G_mat_target <- t(sapply(target$stime,
                             function(ti) sapply(g_funcs, function(gj) gj(ti))))
  }
  time_grid_target  <- seq(min(target$stime), max(target$stime), length.out = 500)
  g_time_mat_target <- if (p == 1) {
    matrix(g_funcs[[1]](time_grid_target), ncol = 1)
  } else {
    tmp <- sapply(g_funcs, function(gj) gj(time_grid_target))
    if (is.vector(tmp)) matrix(tmp, ncol = p) else tmp
  }
  
  delta0 <- rep(0, L + p)  # [zeta(1:p), eta(1:L)]
  deb_loss <- function(d) debias_loss_fast(d, lambda_zeta, lambda_eta, L, p,
                                           target, G_mat_target, time_grid_target,
                                           g_time_mat_target, gamma_hat_A, beta_hat_A)
  deb_grad <- function(d) debias_grad_fast(d, lambda_zeta, lambda_eta, L, p,
                                           target, G_mat_target, time_grid_target,
                                           g_time_mat_target, gamma_hat_A, beta_hat_A)
  
  fit_debias <- optim(par = delta0, fn = deb_loss, gr = deb_grad,
                      method = "L-BFGS-B",
                      control = list(maxit = 2000, factr = 1e-10, pgtol = 1e-6))
  
  delta_hat <- fit_debias$par
  zeta_hat  <- delta_hat[1:p]
  eta_hat   <- delta_hat[(p + 1):(L + p)]
  
  list(
    gamma_hat_A = gamma_hat_A,
    beta_hat_A  = beta_hat_A,
    zeta_hat    = zeta_hat,
    eta_hat     = eta_hat,
    gamma.hat   = gamma_hat_A + zeta_hat,
    beta.hat    = beta_hat_A  + eta_hat,
    g_funcs     = g_funcs
  )
}

# ---------------- Baseline hazard and cumulative hazard ----------------

make_baseline_from_gamma <- function(gamma.hat, g_funcs, t0) {
  stopifnot(is.numeric(gamma.hat), is.list(g_funcs))
  if (length(gamma.hat) != length(g_funcs)) {
    stop("Length of gamma.hat must match length of g_funcs.")
  }
  if (!is.numeric(t0)) stop("t0 must be a numeric vector.")
  
  t0 <- sort(unique(as.numeric(t0)))
  if (length(t0) < 2L) {
    stop("t0 must contain at least two distinct time points for numerical integration.")
  }
  
  G <- sapply(g_funcs, function(gf) {
    val <- gf(t0)
    if (!is.numeric(val)) stop("Each element of g_funcs must return a numeric vector.")
    as.numeric(val)
  })
  if (is.null(dim(G))) G <- matrix(G, ncol = 1)
  if (ncol(G) != length(gamma.hat)) {
    stop("Number of columns of G must match length of gamma.hat.")
  }
  
  logh0 <- drop(G %*% gamma.hat)
  h0    <- exp(logh0)
  
  dt <- diff(t0)
  if (any(dt <= 0)) stop("t0 must be strictly increasing.")
  H0 <- c(0, cumsum(dt * (head(h0, -1) + tail(h0, -1)) / 2))
  
  data.frame(
    time  = t0,
    logh0 = logh0,
    h0    = h0,
    H0    = H0,
    row.names = NULL
  )
}
