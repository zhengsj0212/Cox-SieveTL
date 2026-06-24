# Hoffman2 job-array driver for sieveTL spline simulation
# This script runs ONE (or a small batch) replicate(s) based on SGE_TASK_ID.


# --- Transfer Loss Function ---
transfer_loss_fast <- function(theta, p, L,
                               all.x, all.stime, all.type,
                               G_mat, time_grid, g_time_mat) {
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

# --- Transfer Gradient Function ---
# Smooth L1 (pseudo-Huber) so the objective/gradient stay stable
# --- helpers for smoothed L1 (pseudo-Huber) ---
.smooth_abs       <- function(u, eps) sqrt(u*u + eps*eps) - eps
.smooth_abs_grad  <- function(u, eps) u / sqrt(u*u + eps*eps)

# ---------------- transfer gradient (no penalty) ----------------
transfer_grad_fast <- function(theta, p, L,
                               all.x, all.stime, all.type,
                               G_mat, time_grid, g_time_mat) {
  # theta = [gamma(1:p), beta(1:L)]
  gamma <- theta[1:p]
  beta  <- theta[(p + 1):(L + p)]
  
  xb      <- as.vector(all.x %*% beta)
  exp_xb  <- exp(xb)
  eg_time <- as.vector(g_time_mat %*% gamma)   # length(time_grid)
  exp_eg  <- exp(eg_time)
  
  # ∫_0^{T_i} exp(gamma' g(t)) dt
  int_vec <- vapply(all.stime, function(Ti){
    idx <- which(time_grid <= Ti)
    if (length(idx) < 2) return(0)
    pracma::trapz(time_grid[idx], exp_eg[idx])
  }, numeric(1))
  
  # For each i, j: ∫ g_j(t) exp(gamma' g(t)) dt
  term_mat <- matrix(0, nrow = length(all.stime), ncol = p)
  for (i in seq_along(all.stime)) {
    idx <- which(time_grid <= all.stime[i])
    if (length(idx) < 2) next
    for (j in 1:p) {
      term_mat[i, j] <- pracma::trapz(time_grid[idx], g_time_mat[idx, j] * exp_eg[idx])
    }
  }
  
  grad_gamma <- -colMeans(G_mat * all.type, na.rm = TRUE) +
    colMeans(exp_xb * term_mat, na.rm = TRUE)
  grad_beta  <- -colMeans(all.x * all.type,  na.rm = TRUE) +
    colMeans(all.x * (exp_xb * int_vec), na.rm = TRUE)
  
  c(grad_gamma, grad_beta)
}

# ---------------- debias loss (with smoothed L1) ----------------
debias_loss_fast <- function(delta, lambda_zeta, lambda_eta, L, p,
                             target, G_mat_target, time_grid_target,
                             g_time_mat_target, gamma_hat_A, beta_hat_A,
                             eps = 1e-4) {
  # delta = [zeta(1:p), eta(1:L)]
  zeta <- delta[1:p]
  eta  <- delta[(p + 1):(L + p)]
  
  gamma <- gamma_hat_A + zeta
  beta  <- beta_hat_A  + eta
  
  n       <- nrow(target$x)
  xb      <- as.vector(target$x %*% beta)
  g_gamma <- if (p == 1) as.vector(G_mat_target * gamma) else as.vector(G_mat_target %*% gamma)
  term1   <- target$type * (xb + g_gamma)
  
  eg_time <- as.vector(g_time_mat_target %*% gamma)
  exp_eg  <- exp(eg_time)
  
  int_vec <- vapply(target$stime, function(Ti){
    idx <- which(time_grid_target <= Ti)
    if (length(idx) < 2) return(0)
    pracma::trapz(time_grid_target[idx], exp_eg[idx])
  }, numeric(1))
  
  term2 <- exp(xb) * int_vec
  nll   <- -mean(term1 - term2, na.rm = TRUE)             # averaged NLL
  
  # smoothed L1 penalties, scaled by n to match averaged NLL
  pen  <- (lambda_zeta ) * sum(.smooth_abs(zeta, eps)) +
    (lambda_eta  ) * sum(.smooth_abs(eta,  eps))
  
  nll + pen
}

# ---------------- debias gradient (with smoothed L1) ----------------
debias_grad_fast <- function(delta, lambda_zeta, lambda_eta, L, p,
                             target, G_mat_target, time_grid_target,
                             g_time_mat_target, gamma_hat_A, beta_hat_A,
                             eps = 1e-4) {
  zeta <- delta[1:p]
  eta  <- delta[(p + 1):(L + p)]
  
  gamma <- gamma_hat_A + zeta
  beta  <- beta_hat_A  + eta
  
  n       <- nrow(target$x)
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
      term_mat[i, j] <- pracma::trapz(time_grid_target[idx], g_time_mat_target[idx, j] * exp_eg[idx])
    }
  }
  
  grad_zeta <- -colMeans(G_mat_target * target$type, na.rm = TRUE) +
    colMeans(exp_xb * term_mat, na.rm = TRUE)
  grad_eta  <- -colMeans(target$x * target$type,  na.rm = TRUE) +
    colMeans(target$x * (exp_xb * int_vec), na.rm = TRUE)
  
  # add smoothed L1 gradients (scaled by n)
  grad_zeta <- grad_zeta + (lambda_zeta ) * .smooth_abs_grad(zeta, eps)
  grad_eta  <- grad_eta  + (lambda_eta  ) * .smooth_abs_grad(eta,  eps)
  
  c(grad_zeta, grad_eta)
}

soft_thresh <- function(z, a) {
  sign(z) * pmax(abs(z) - a, 0)
}

run_one_sparse_iteration <- function(delta_old, theta_hat_A,
                                     lambda_zeta, lambda_eta,
                                     L, p,
                                     target, g_funcs,
                                     G_mat_target, time_grid_target, g_time_mat_target,
                                     gamma_hat_A, beta_hat_A,
                                     iter_update_method = "diag_soft_threshold",
                                     ridge_iter = 1e-6) {
  gamma_old <- gamma_hat_A + delta_old[1:p]
  beta_old <- beta_hat_A + delta_old[(p + 1):(p + L)]

  grad <- debias_grad_fast(
    delta = delta_old,
    lambda_zeta = 0,
    lambda_eta = 0,
    L = L,
    p = p,
    target = target,
    G_mat_target = G_mat_target,
    time_grid_target = time_grid_target,
    g_time_mat_target = g_time_mat_target,
    gamma_hat_A = gamma_hat_A,
    beta_hat_A = beta_hat_A
  )
  H <- compute_target_hessian_theta(
    data_k = target,
    gamma_hat = gamma_old,
    beta_hat = beta_old,
    g_funcs = g_funcs,
    time_grid = time_grid_target
  )
  d <- length(delta_old)
  H_reg <- H + ridge_iter * diag(d)
  lambda_vec <- c(rep(lambda_zeta, p), rep(lambda_eta, L))
  Hdiag <- pmax(diag(H_reg), 1e-12)
  z <- delta_old - grad / Hdiag
  delta_iter <- soft_thresh(z, lambda_vec / Hdiag)

  theta_iter <- theta_hat_A + delta_iter
  list(
    grad = grad,
    H = H,
    H_reg = H_reg,
    delta_iter = delta_iter,
    theta_iter = theta_iter,
    gamma_hat_iter = theta_iter[1:p],
    beta_hat_iter = theta_iter[(p + 1):(p + L)]
  )
}

default_penalty_reference_d <- 13
default_penalty_reference_n0 <- 100

penalty_scale_factor <- function(d, n0) {
  if (!is.finite(d) || d < 1) stop("d must be positive for penalty scaling.")
  if (!is.finite(n0) || n0 < 1) stop("n0 must be positive for penalty scaling.")
  sqrt(log(d) / n0)
}

lambda_to_penalty_c <- function(lambda, d_ref = default_penalty_reference_d,
                                n0_ref = default_penalty_reference_n0) {
  lambda / penalty_scale_factor(d_ref, n0_ref)
}

# ---------------- main wrapper ----------------
# ============== sievetl_approx（g_funcs 按 target$stime 定义）================
sievetl_approx <- function(target, source = NULL,
                            L = 2,
                            lambda_zeta = 0.05, lambda_eta = 0.05,
                            c_zeta = NULL, c_eta = NULL,
                            n0_penalty = NULL,
                            c = 1, c_sieve = NULL,
                            transfer_start = NULL,
                            debias_start = NULL,
                            iter_update_method = "diag_soft_threshold",
                            ridge_iter = 1e-6) {
  if (!is.null(c_sieve)) c <- c_sieve
  transfer.source.id <- seq_along(source)
  
  # ---------- stack data ----------
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
  
  # ---------- spline basis (pooled knots) ----------
  # N distinct pooled EVENT times, J = floor(c * N^(1/3))
  # Internal knots at pooled time quantiles, boundary knots at min/max pooled time
  # p_hat = J + 4 for cubic (degree=3) spline with intercept
  pooled_stime <- all.stime
  pooled_event_times <- pooled_stime[all.type == 1]
  N <- length(unique(pooled_event_times))
  J <- floor(c * N^(1/3))
  if (J < 1) J <- 0
  knots <- if (J > 0) {
    stats::quantile(pooled_stime, probs = (1:J) / (J + 1),
                    type = 7, names = FALSE)
  } else {
    NULL
  }
  bkn <- range(pooled_stime)
  deg <- 3
  
  B0 <- splines::bs(pooled_stime, degree = deg, knots = knots,
                    Boundary.knots = bkn, intercept = TRUE)
  p <- ncol(B0)
  if (is.null(n0_penalty)) n0_penalty <- nrow(target$x)
  d_penalty <- p + L
  lambda_scale <- penalty_scale_factor(d_penalty, n0_penalty)
  if (!is.null(c_zeta)) lambda_zeta <- c_zeta * lambda_scale
  if (!is.null(c_eta)) lambda_eta <- c_eta * lambda_scale
  c_zeta_used <- if (is.null(c_zeta)) lambda_zeta / lambda_scale else c_zeta
  c_eta_used <- if (is.null(c_eta)) lambda_eta / lambda_scale else c_eta
  
  g_funcs <- lapply(seq_len(ncol(B0)), function(j){
    function(t){
      Bt <- splines::bs(t, degree = deg, knots = knots,
                        Boundary.knots = bkn, intercept = TRUE)
      as.numeric(Bt[, j])
    }
  })
  
  # ---------- G_mat over ALL samples ----------
  if (p == 1) {
    G_mat <- matrix(sapply(all.stime, function(ti) g_funcs[[1]](ti)), ncol = 1)
  } else {
    G_mat <- t(sapply(all.stime, function(ti) sapply(g_funcs, function(gj) gj(ti))))
  }
  
  # grids for integration (combined range)
  time_grid  <- seq(min(all_time), max(all_time), length.out = 500)
  g_time_mat <- if (p == 1) {
    matrix(g_funcs[[1]](time_grid), ncol = 1)
  } else {
    tmp <- sapply(g_funcs, function(gj) gj(time_grid))
    if (is.vector(tmp)) matrix(tmp, ncol = p) else tmp
  }
  
  # ---------- transfer step (unpenalized) ----------
  theta0 <- if (!is.null(transfer_start) && length(transfer_start) == (L + p)) {
    as.numeric(transfer_start)
  } else {
    rep(0, L + p)
  }  # [gamma(1:p), beta(1:L)]
  loss_wrap <- function(th)
    transfer_loss_fast(th, p, L, all.x, all.stime, all.type, G_mat, time_grid, g_time_mat)
  grad_wrap <- function(th)
    transfer_grad_fast(th, p, L, all.x, all.stime, all.type, G_mat, time_grid, g_time_mat)
  
  fit <- nlminb(start = theta0,
                objective = loss_wrap,
                gradient  = grad_wrap,
                control   = list(iter.max = 5000, eval.max = 5000))
  
  theta_hat_A <- fit$par
  gamma_hat_A <- theta_hat_A[1:p]
  beta_hat_A  <- theta_hat_A[(p + 1):(L + p)]
  
  # ---------- debias step (with smoothed L1) ----------
  # build target-specific design using the SAME spline params learned from all_time
  if (p == 1) {
    G_mat_target <- matrix(g_funcs[[1]](target$stime), ncol = 1)
  } else {
    G_mat_target <- t(sapply(target$stime, function(ti) sapply(g_funcs, function(gj) gj(ti))))
  }
  time_grid_target  <- seq(min(target$stime), max(target$stime), length.out = 500)
  g_time_mat_target <- if (p == 1) {
    matrix(g_funcs[[1]](time_grid_target), ncol = 1)
  } else {
    tmp <- sapply(g_funcs, function(gj) gj(time_grid_target))
    if (is.vector(tmp)) matrix(tmp, ncol = p) else tmp
  }
  
  delta0 <- if (!is.null(debias_start) && length(debias_start) == (L + p)) {
    as.numeric(debias_start)
  } else {
    rep(0, L + p)
  }  # [zeta(1:p), eta(1:L)]
  deb_loss <- function(d) debias_loss_fast(d, lambda_zeta, lambda_eta, L, p,
                                           target, G_mat_target, time_grid_target,
                                           g_time_mat_target, gamma_hat_A, beta_hat_A)
  deb_grad <- function(d) debias_grad_fast(d, lambda_zeta, lambda_eta, L, p,
                                           target, G_mat_target, time_grid_target,
                                           g_time_mat_target, gamma_hat_A, beta_hat_A)
  
  fit_debias <- optim(par = delta0, fn = deb_loss, gr = deb_grad,
                      method = "L-BFGS-B",
                      control = list(maxit = 2000, factr = 1e-10, pgtol = 1e-6))
  
  delta_hat_optim <- fit_debias$par
  zeta_hat_optim  <- delta_hat_optim[1:p]
  eta_hat_optim   <- delta_hat_optim[(p + 1):(L + p)]
  theta_hat_optim <- theta_hat_A + delta_hat_optim
  gamma_hat_optim <- theta_hat_optim[1:p]
  beta_hat_optim  <- theta_hat_optim[(p + 1):(L + p)]

  iter_fit <- run_one_sparse_iteration(
    delta_old = delta_hat_optim,
    theta_hat_A = theta_hat_A,
    lambda_zeta = lambda_zeta,
    lambda_eta = lambda_eta,
    L = L,
    p = p,
    target = target,
    g_funcs = g_funcs,
    G_mat_target = G_mat_target,
    time_grid_target = time_grid_target,
    g_time_mat_target = g_time_mat_target,
    gamma_hat_A = gamma_hat_A,
    beta_hat_A = beta_hat_A,
    iter_update_method = iter_update_method,
    ridge_iter = ridge_iter
  )
  theta_hat_iter <- iter_fit$theta_iter
  gamma_hat_iter <- iter_fit$gamma_hat_iter
  beta_hat_iter  <- iter_fit$beta_hat_iter
  delta_hat_iter <- theta_hat_iter - theta_hat_A
  zeta_hat_iter  <- delta_hat_iter[1:p]
  eta_hat_iter   <- delta_hat_iter[(p + 1):(L + p)]

  list(
    gamma_hat_A = gamma_hat_A,
    beta_hat_A  = beta_hat_A,
    theta_hat_A = theta_hat_A,
    zeta_hat_optim = zeta_hat_optim,
    eta_hat_optim = eta_hat_optim,
    delta_hat_optim = delta_hat_optim,
    theta_hat_optim = theta_hat_optim,
    gamma_hat_optim = gamma_hat_optim,
    beta_hat_optim = beta_hat_optim,
    zeta_hat = zeta_hat_iter,
    eta_hat = eta_hat_iter,
    delta_hat = delta_hat_iter,
    theta_hat = theta_hat_iter,
    gamma.hat = gamma_hat_iter,
    beta.hat = beta_hat_iter,
    gamma_hat = gamma_hat_iter,
    beta_hat = beta_hat_iter,
    lambda_zeta = lambda_zeta,
    lambda_eta = lambda_eta,
    c_zeta = c_zeta_used,
    c_eta = c_eta_used,
    d_penalty = d_penalty,
    n0_penalty = n0_penalty,
    p = p,
    g_funcs = g_funcs,
    estimation_engine = "old_optim_plus_one_sparse_iteration",
    iter_update_method = iter_update_method,
    ridge_iter = ridge_iter
  )
}



# ============== sievetl_ic_approx（样条基也只用 target1$stime 定义）================
sievetl_ic_approx <- function(
    target1, source1 = NULL,
    L = 6,
    c_grid = c(0.5, 1),
    c_zeta_grid = c(0, 0.05, 0.1, 0.2),
    c_eta_grid  = c(0, 0.05, 0.1, 0.2),
    lambda_fixed_stage1 = 0,
    tie_tol = 1e-8,               # AIC/BIC 并列容差
    tie_break = c("lz_small_le_large","both_large","both_small")  # 并列偏好
){
  tie_break  <- match.arg(tie_break)
  n <- nrow(target1$x)
  
  # ---------- combined stimes for spline knots (target + all sources) ----------
  combined_stimes <- {
    if (is.null(source1) || length(source1) == 0) {
      target1$stime
    } else {
      c(target1$stime, unlist(lapply(source1, function(s) s$stime)))
    }
  }
  combined_types <- {
    if (is.null(source1) || length(source1) == 0) {
      target1$type
    } else {
      c(target1$type, unlist(lapply(source1, function(s) s$type)))
    }
  }
  finite_mask <- is.finite(combined_stimes)
  combined_stimes <- combined_stimes[finite_mask]
  combined_types <- combined_types[finite_mask]
  
  # ---- helpers ----
  .get <- function(lst, ...) {
    keys <- c(...)
    for (k in keys) if (!is.null(lst[[k]])) return(lst[[k]])
    return(NULL)
  }
  
  build_g_funcs <- function(stimes_ref, types_ref, c = 1){
    # ---- spline knots / sieve dimension (pooled data) ----
    # N distinct pooled EVENT times, J = floor(c * N^(1/3))
    # Internal knots at pooled time quantiles, boundary knots at min/max pooled time
    # p_hat = J + 4 for cubic (degree=3) spline with intercept
    pooled_stime <- stimes_ref
    pooled_event_times <- pooled_stime[types_ref == 1]
    N <- length(unique(pooled_event_times))
    J <- floor(c * N^(1/3))
    if (J < 1) J <- 0
    knots <- if (J > 0) {
      stats::quantile(pooled_stime, probs = (1:J) / (J + 1),
                      type = 7, names = FALSE)
    } else {
      NULL
    }
    bkn <- range(pooled_stime)
    deg <- 3
    B0 <- splines::bs(pooled_stime, degree = deg, knots = knots,
                      Boundary.knots = bkn, intercept = TRUE)
    g_funcs <- lapply(seq_len(ncol(B0)), function(j){
      function(t){
        Bt <- splines::bs(t, degree = deg, knots = knots,
                          Boundary.knots = bkn, intercept = TRUE)
        as.numeric(Bt[, j])
      }
    })
    g_funcs
  }
  
  # ---- Stage 1：用“最终系数”（beta.hat/gamma.hat）计算未惩罚 nll 与 df_{beta,gamma} ----
  eval_ll_df_stage1 <- function(fit, g_funcs, target1){
    beta  <- .get(fit, "beta.hat",  "beta_hat_A",  "beta")
    gamma <- .get(fit, "gamma.hat", "gamma_hat_A", "gamma")
    stopifnot(!is.null(beta), !is.null(gamma))
    
    xb <- as.vector(target1$x %*% beta)
    p  <- length(gamma)
    
    # G_mat (n x p)
    if (p == 1){
      G_vec <- vapply(target1$stime, function(ti) g_funcs[[1]](ti), numeric(1))
      G_mat <- matrix(G_vec, ncol = 1)
    } else {
      G_mat <- t(vapply(
        target1$stime,
        function(ti) vapply(g_funcs, function(gj) gj(ti), numeric(1)),
        numeric(p)
      ))
    }
    g_gamma <- as.vector(G_mat %*% matrix(gamma, ncol = 1))
    
    term1 <- target1$type * (xb + g_gamma)
    
    # ∫_0^{Ti} exp(gamma^T g(t)) dt
    time_grid   <- seq(min(target1$stime), max(target1$stime), length.out = 500)
    g_time_mat  <- sapply(g_funcs, function(gj) gj(time_grid))
    if (is.vector(g_time_mat)) g_time_mat <- matrix(g_time_mat, ncol = p)
    exp_gamma_g <- exp(as.vector(g_time_mat %*% gamma))
    
    int_vec <- vapply(target1$stime, function(Ti){
      idx <- which(time_grid <= Ti)
      if (length(idx) < 2) return(0)
      pracma::trapz(time_grid[idx], exp_gamma_g[idx])
    }, numeric(1))
    
    term2  <- exp(xb) * int_vec
    loglik <- mean(term1 - term2, na.rm = TRUE)
    nll    <- -loglik
    
    # df_bg：最终（beta.hat, gamma.hat）的“非零计数”，只把 exact 0 当作零
    df_bg <- sum(beta != 0) + sum(gamma != 0)
    
    list(nll = nll, loglik = loglik, df_bg = df_bg)
  }
  
  # Stage 2 的 df：只数 (zeta, eta)
  eval_df_stage2 <- function(fit){
    zeta <- .get(fit, "zeta_hat", "zeta.hat", "zeta")
    eta  <- .get(fit, "eta_hat",  "eta.hat",  "eta")
    if (is.null(zeta) && is.null(eta)) return(NA_integer_)
    zc <- if (!is.null(zeta)) sum(zeta != 0) else 0L
    ec <- if (!is.null(eta))  sum(eta != 0) else 0L
    zc + ec
  }
  
  search_grid <- expand.grid(
    c = c_grid,
    c_zeta = c_zeta_grid,
    c_eta = c_eta_grid
  )
  search_grid$BIC <- NA_real_
  search_grid$nll <- NA_real_
  search_grid$df_bg <- NA_integer_
  search_grid$df_ze <- NA_integer_
  search_grid$p_hat <- NA_integer_
  search_grid$d_penalty <- NA_integer_
  search_grid$lambda_zeta <- NA_real_
  search_grid$lambda_eta <- NA_real_
  search_grid <- search_grid[order(search_grid$c, search_grid$c_zeta, search_grid$c_eta), , drop = FALSE]

  warm_starts <- setNames(vector("list", length(c_grid)), as.character(c_grid))

  for (k in seq_len(nrow(search_grid))) {
    c_val <- search_grid$c[k]
    c_zeta_val <- search_grid$c_zeta[k]
    c_eta_val <- search_grid$c_eta[k]
    c_key <- as.character(c_val)
    start_info <- warm_starts[[c_key]]

    fit <- tryCatch(
      sievetl_approx(
        target = target1, source = source1,
        L = L,
        c_zeta = c_zeta_val, c_eta = c_eta_val,
        n0_penalty = n,
        c = c_val,
        transfer_start = if (!is.null(start_info)) start_info$theta else NULL,
        debias_start = if (!is.null(start_info)) start_info$delta else NULL
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) next

    warm_starts[[c_key]] <- list(
      theta = fit$theta_hat_A,
      delta = fit$delta_hat_optim
    )

    ed1 <- eval_ll_df_stage1(fit, fit$g_funcs, target1)
    df_ze <- eval_df_stage2(fit)

    search_grid$nll[k] <- ed1$nll
    search_grid$df_bg[k] <- ed1$df_bg
    search_grid$df_ze[k] <- df_ze
    search_grid$BIC[k] <- ed1$nll + (log(n) / n) * df_ze
    search_grid$p_hat[k] <- fit$p
    search_grid$d_penalty[k] <- fit$p + L
    search_grid$lambda_zeta[k] <- fit$lambda_zeta
    search_grid$lambda_eta[k] <- fit$lambda_eta
  }

  minBIC <- min(search_grid$BIC, na.rm = TRUE)
  ties <- subset(search_grid, !is.na(BIC) & (BIC <= minBIC + tie_tol))
  
  if (nrow(ties) > 1) {
    if (tie_break == "lz_small_le_large") {
      ties <- ties[order(ties$BIC, ties$c_zeta, -ties$c_eta, ties$c), ]
    } else if (tie_break == "both_large") {
      ties <- ties[order(ties$BIC, -ties$c_zeta, -ties$c_eta, -ties$c), ]
    } else { # both_small
      ties <- ties[order(ties$BIC, ties$c_zeta, ties$c_eta, ties$c), ]
    }
  }

  valid_rows <- search_grid[is.finite(search_grid$BIC), , drop = FALSE]
  aic_tbl <- if (nrow(valid_rows) == 0L) {
    data.frame(c = numeric(0), p_hat = integer(0), d_penalty = integer(0),
               nll = numeric(0), df_bg = integer(0), BIC = numeric(0))
  } else {
    split_rows <- split(valid_rows, valid_rows$c, drop = TRUE)
    do.call(rbind, lapply(split_rows, function(df) {
      df <- df[order(df$BIC, df$c_zeta, df$c_eta), , drop = FALSE]
      df[1, c("c", "p_hat", "d_penalty", "nll", "df_bg", "BIC"), drop = FALSE]
    }))
  }

  list(
    best_p = ties$p_hat[1],
    best_d = ties$d_penalty[1],
    best_c = ties$c[1],
    best_c_zeta = ties$c_zeta[1],
    best_c_eta = ties$c_eta[1],
    best_lambda_zeta = ties$lambda_zeta[1],
    best_lambda_eta = ties$lambda_eta[1],
    aic_table = aic_tbl,
    bic_table = search_grid
  )
}



simulate_once_all <- function(n_target = 500, n_source = 2000,
                              beta2_source = 0.5, kappa_source = 1.5,
                              c1 = 1, covariate_shift = TRUE, L = 5) {
  if (!(L %in% c(5, 6))) stop("Only L = 5 or L = 50 is supported.")
  
  # --- Target Data ---
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
  ctime_t <- rexp(n_target, rate = 1 / c1)
  stime_t <- pmin(etime_t, ctime_t)
  type_t <- as.numeric(etime_t <= ctime_t)
  
  target <- list(x = x_t, stime = stime_t, type = type_t)
  
  # --- Source Data ---
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
  ctime_s <- rexp(n_source, rate = 1 / c1)
  stime_s <- pmin(etime_s, ctime_s)
  type_s <- as.numeric(etime_s <= ctime_s)
  
  source <- list(list(x = X_s, stime = stime_s, type = type_s))
  
  return(list(target = target, source = source))
}

make_baseline_from_gamma <- function(gamma.hat, g_funcs, t0) {
  # 1) 基本检查
  stopifnot(is.numeric(gamma.hat), is.list(g_funcs))
  if (length(gamma.hat) != length(g_funcs)) {
    stop("gamma.hat 的长度需要与 g_funcs 的长度一致。")
  }
  if (!is.numeric(t0)) stop("t0 必须是数值向量。")
  t0 <- sort(unique(as.numeric(t0)))
  if (length(t0) < 2L) stop("t0 至少需要包含两个不同的时间点以便数值积分。")
  
  # 2) 在 t0 上计算基函数矩阵 G：n_t0 × p
  G <- sapply(g_funcs, function(gf) {
    val <- gf(t0)
    if (!is.numeric(val)) stop("g_funcs 里的函数必须返回数值向量。")
    as.numeric(val)
  })
  # 若只有一个基函数，sapply 会退化成向量，这里强制成矩阵
  if (is.null(dim(G))) G <- matrix(G, ncol = 1)
  
  if (ncol(G) != length(gamma.hat)) {
    stop("在 t0 上评估后的基函数列数与 gamma.hat 长度不一致。")
  }
  
  # 3) 计算 log h0(t), h0(t)
  logh0 <- drop(G %*% gamma.hat)
  h0    <- exp(logh0)
  
  # 4) 梯形积分得到 H0(t)（累积基线风险）
  dt <- diff(t0)
  if (any(dt <= 0)) stop("t0 需要严格递增。")
  # 梯形：∫ h0 ≈ Σ (h0_{i}+h0_{i+1})/2 * Δt_i
  H0 <- c(0, cumsum(dt * (head(h0, -1) + tail(h0, -1)) / 2))
  
  # 返回
  data.frame(
    time  = t0,
    logh0 = logh0,
    h0    = h0,
    H0    = H0,
    row.names = NULL
  )
}

compute_cumhaz_if_curve <- function(gamma_hat, if_theta_mat, g_funcs, t0) {
  base_curve <- make_baseline_from_gamma(gamma_hat, g_funcs, t0)
  p <- length(gamma_hat)
  if_gamma_mat <- if_theta_mat[, seq_len(p), drop = FALSE]

  G <- sapply(g_funcs, function(gf) as.numeric(gf(base_curve$time)))
  if (is.null(dim(G))) G <- matrix(G, ncol = 1)
  h0 <- as.numeric(base_curve$h0)
  dt <- diff(base_curve$time)

  a_mat <- matrix(0, nrow = nrow(base_curve), ncol = p)
  for (j in seq_len(p)) {
    integrand <- G[, j] * h0
    a_mat[, j] <- c(0, cumsum(dt * (head(integrand, -1) + tail(integrand, -1)) / 2))
  }

  mean_if_gamma <- colMeans(if_gamma_mat)
  H0_if <- as.numeric(base_curve$H0 + a_mat %*% mean_if_gamma)

  data.frame(
    time = base_curve$time,
    H0 = H0_if,
    row.names = NULL
  )
}

compute_a_matrix <- function(gamma_eval, g_funcs, eval_times) {
  eval_times <- sort(unique(as.numeric(eval_times)))
  base_curve <- make_baseline_from_gamma(gamma_eval, g_funcs, eval_times)
  p <- length(gamma_eval)
  G <- sapply(g_funcs, function(gf) as.numeric(gf(base_curve$time)))
  if (is.null(dim(G))) G <- matrix(G, ncol = 1)
  h0 <- as.numeric(base_curve$h0)
  dt <- diff(base_curve$time)

  a_mat <- matrix(0, nrow = nrow(base_curve), ncol = p)
  for (j in seq_len(p)) {
    integrand <- G[, j] * h0
    a_mat[, j] <- c(0, cumsum(dt * (head(integrand, -1) + tail(integrand, -1)) / 2))
  }

  list(time = base_curve$time, a_mat = a_mat, base_curve = base_curve)
}

interpolate_a_matrix <- function(a_obj, xout) {
  xout <- as.numeric(xout)
  a_interp <- vapply(seq_len(ncol(a_obj$a_mat)), function(j) {
    stats::approx(a_obj$time, a_obj$a_mat[, j], xout = xout, rule = 2, ties = "ordered")$y
  }, numeric(length(xout)))
  if (is.null(dim(a_interp))) a_interp <- matrix(a_interp, ncol = 1)
  rownames(a_interp) <- NULL
  a_interp
}

compute_target_score_theta <- function(data, gamma_hat, beta_hat, g_funcs, time_grid) {
  n <- nrow(data$x)
  L <- length(beta_hat)
  p <- length(gamma_hat)
  d <- p + L  # order: theta = (gamma, beta)
  
  grad <- rep(0, d)
  
  for (i in 1:n) {
    x_i <- data$x[i, ]
    t_i <- data$stime[i]
    delta_i <- data$type[i]
    
    if (t_i <= 0) next
    
    # Event term
    g_t <- sapply(g_funcs, function(f) f(t_i))  # p-vector
    z_t <- c(g_t, x_i)  # order: gamma then beta
    grad <- grad - delta_i * z_t
    
    # Integral term
    idx <- which(time_grid <= t_i)
    if (length(idx) < 2) next
    
    t_grid_i <- time_grid[idx]
    G_mat <- t(sapply(g_funcs, function(f) f(t_grid_i)))  # p x m
    m <- length(t_grid_i)
    
    eta_g <- as.vector(gamma_hat %*% G_mat)   # length m
    eta_x <- sum(beta_hat * x_i)
    exp_eta <- exp(eta_g + eta_x)  # length m
    
    # Z_mat: (p+L) x m, gamma then beta
    Z_mat <- rbind(
      G_mat,
      matrix(rep(x_i, times = m), nrow = L)
    )
    
    integral <- apply(Z_mat * rep(exp_eta, each = d), 1, function(row) trapz(t_grid_i, row))
    grad <- grad + integral
  }
  
  grad / n
}

compute_target_score_theta_subjects <- function(data, gamma_hat, beta_hat, g_funcs, time_grid) {
  n <- nrow(data$x)
  L <- length(beta_hat)
  p <- length(gamma_hat)
  d <- p + L  # order: theta = (gamma, beta)

  score_mat <- matrix(0, nrow = n, ncol = d)

  for (i in 1:n) {
    x_i <- data$x[i, ]
    t_i <- data$stime[i]
    delta_i <- data$type[i]

    if (t_i <= 0) next

    g_t <- sapply(g_funcs, function(f) f(t_i))
    z_t <- c(g_t, x_i)  # order: gamma then beta
    score_i <- -delta_i * z_t

    idx <- which(time_grid <= t_i)
    if (length(idx) >= 2) {
      t_grid_i <- time_grid[idx]
      G_mat <- t(sapply(g_funcs, function(f) f(t_grid_i)))
      m <- length(t_grid_i)

      eta_g <- as.vector(gamma_hat %*% G_mat)
      eta_x <- sum(beta_hat * x_i)
      exp_eta <- exp(eta_g + eta_x)

      Z_mat <- rbind(
        G_mat,
        matrix(rep(x_i, times = m), nrow = L)
      )

      integral <- apply(Z_mat * rep(exp_eta, each = d), 1, function(row) trapz(t_grid_i, row))
      score_i <- score_i + integral
    }

    score_mat[i, ] <- score_i
  }

  score_mat
}

compute_target_hessian_theta <- function(data_k, gamma_hat, beta_hat, g_funcs, time_grid) {
  n_k <- nrow(data_k$x)
  p <- length(gamma_hat)
  L <- length(beta_hat)
  d <- p + L
  
  hessian <- matrix(0, d, d)
  
  for (i in 1:n_k) {
    x_i <- data_k$x[i, ]
    t_i <- data_k$stime[i]
    if (t_i <= 0) next
    
    t_grid_i <- time_grid[time_grid <= t_i]
    if (length(t_grid_i) < 2) next
    
    G_mat <- t(sapply(g_funcs, function(f) f(t_grid_i)))  # p x m
    m <- length(t_grid_i)
    
    eta_g_vec <- as.vector(gamma_hat %*% G_mat)
    exp_eta_g <- exp(eta_g_vec)
    
    eta_x <- sum(beta_hat * x_i)
    exp_eta_x <- exp(eta_x)
    
    integrand_array <- array(0, dim = c(d, d, m))
    for (j in 1:m) {
      g_j <- G_mat[, j]
      z_j <- c(g_j, x_i)  # order: gamma then beta
      integrand_array[, , j] <- tcrossprod(z_j) * exp_eta_g[j]
    }
    
    integral <- apply(integrand_array, c(1, 2), function(v) trapz(t_grid_i, v))
    hessian <- hessian + exp_eta_x * integral
  }
  
  hessian / n_k
}

combine_domains <- function(target, source_list = NULL) {
  if (is.null(source_list) || length(source_list) == 0L) return(target)
  list(
    x = do.call(rbind, c(list(target$x), lapply(source_list, `[[`, "x"))),
    stime = c(target$stime, unlist(lapply(source_list, `[[`, "stime"))),
    type = c(target$type, unlist(lapply(source_list, `[[`, "type")))
  )
}

compute_pooled_hessian_theta <- function(target, source_list, gamma_hat, beta_hat, g_funcs, time_grid) {
  pooled_data <- combine_domains(target, source_list)
  compute_target_hessian_theta(pooled_data, gamma_hat, beta_hat, g_funcs, time_grid)
}

fit_inverse_hessian_initializer_unpenalized <- function(H_mat, j, ridge = 1e-8) {
  d <- ncol(H_mat)
  e_j <- rep(0, d)
  e_j[j] <- 1
  tryCatch(
    as.vector(solve(H_mat, e_j)),
    error = function(e) as.vector(solve(H_mat + ridge * diag(d), e_j))
  )
}

fit_inverse_hessian_correction <- function(H_target, omega_A, j, lambda_l1, eps = 1e-4) {
  d <- ncol(H_target)
  e_j <- rep(0, d)
  e_j[j] <- 1
  if (lambda_l1 == 0) {
    omega_ls <- tryCatch(
      as.vector(solve(H_target, e_j)),
      error = function(e) as.vector(qr.solve(H_target, e_j))
    )
    return(omega_ls - omega_A)
  }
  obj <- function(delta) {
    omega <- omega_A + delta
    resid <- as.vector(H_target %*% omega - e_j)
    q <- 0.5 * sum(resid^2)
    q + lambda_l1 * sum(.smooth_abs(delta, eps))
  }
  grad <- function(delta) {
    omega <- omega_A + delta
    resid <- as.vector(H_target %*% omega - e_j)
    as.vector(crossprod(H_target, resid) + lambda_l1 * .smooth_abs_grad(delta, eps))
  }
  fit <- optim(par = rep(0, d), fn = obj, gr = grad, method = "BFGS",
               control = list(maxit = 2000, reltol = 1e-10))
  fit$par
}

estimate_transfer_inverse_hessian_by_column <- function(H_pooled, H_target, lambda_delta_vec, eps = 1e-4) {
  d <- ncol(H_target)
  stopifnot(length(lambda_delta_vec) == d)
  Omega_hat <- matrix(0, nrow = d, ncol = d)
  for (j in seq_len(d)) {
    omega_A <- fit_inverse_hessian_initializer_unpenalized(H_pooled, j)
    delta_j <- fit_inverse_hessian_correction(H_target, omega_A, j, lambda_delta_vec[j], eps = eps)
    Omega_hat[, j] <- omega_A + delta_j
  }
  Omega_hat
}

compute_target_inverse_hessian <- function(H_target, ridge = 1e-8) {
  # Non-transfer one-step benchmark: use only the target Hessian inverse.
  tryCatch(
    solve(H_target),
    error = function(e) solve(H_target + ridge * diag(ncol(H_target)))
  )
}
library(pracma)

# The old plugin IF one-step method has been retired.
# This driver now uses only the joint finite-dimensional sieve parameter
# theta = (gamma, beta) together with the target-domain score/Hessian.

read_env_string <- function(name, default = "") {
  val <- Sys.getenv(name, unset = "")
  if (!nzchar(val)) default else val
}

read_env_integer <- function(name, default) {
  raw <- read_env_string(name, "")
  parsed <- suppressWarnings(as.integer(raw))
  if (length(parsed) == 0L || is.na(parsed)) as.integer(default) else parsed
}

read_env_numeric <- function(name, default) {
  raw <- read_env_string(name, "")
  parsed <- suppressWarnings(as.numeric(raw))
  if (length(parsed) == 0L || is.na(parsed)) as.numeric(default) else parsed
}

format_kappa_tag <- function(kappa) {
  if (isTRUE(all.equal(kappa, round(kappa), tolerance = 1e-8))) {
    return(sprintf("%d", as.integer(round(kappa))))
  }
  formatC(kappa, format = "f", digits = 1)
}

read_by_column_lambda_file <- function(path, d_n) {
  if (!file.exists(path)) {
    stop("final_selected_c_by_column.csv not found: ", path)
  }
  tbl <- read.csv(path, stringsAsFactors = FALSE)
  required_cols <- c("j", "best_c", "best_lambda_value")
  if (!all(required_cols %in% names(tbl))) {
    stop("Missing required columns in by-column lambda file: ", path)
  }
  tbl <- tbl[order(tbl$j), ]
  if (!identical(as.integer(tbl$j), seq_len(d_n))) {
    stop("Column indices in by-column lambda file do not match 1:d_n for file: ", path)
  }
  list(
    c = as.numeric(tbl$best_c),
    lambda = as.numeric(tbl$best_lambda_value)
  )
}

read_sievetl_tuning_file <- function(path) {
  if (!file.exists(path)) {
    stop("SieveTL tuning file not found: ", path)
  }
  tbl <- read.csv(path, stringsAsFactors = FALSE)
  required_cols <- c("c_hat", "c_zeta_hat", "c_eta_hat",
                     "lambda_zeta_hat", "lambda_eta_hat")
  if (!all(required_cols %in% names(tbl))) {
    stop("Missing required columns in sievetl tuning file: ", path)
  }
  if (nrow(tbl) < 1L) {
    stop("Empty sievetl tuning file: ", path)
  }
  tbl[1, , drop = FALSE]
}

default_pilot_design_dir <- function(beta2, kappa, c1, n0, nA, L = 5,
                                     shift = FALSE, base_dir = getwd()) {
  file.path(
    path.expand(base_dir),
    "pilot_parameter_by_column",
    sprintf("n%dN%d", n0, nA),
    sprintf(
      "beta%.1f_kappa%s_L%d_c%.2f_shift%s",
      beta2, format_kappa_tag(kappa), L, c1,
      ifelse(isTRUE(shift), "TRUE", "FALSE")
    )
  )
}

default_sievetl_tuning_path <- function(beta2, kappa, c1, n0, nA, L = 5,
                                        shift = FALSE, base_dir = getwd()) {
  design_dir <- default_pilot_design_dir(beta2, kappa, c1, n0, nA, L, shift, base_dir)
  file.path(design_dir, "final_selected_tuning.csv")
}

# -------- Hoffman2 job-array controls --------
nsim <- read_env_integer("NSIM", 100L)
task_id <- read_env_integer("SGE_TASK_ID", 1L)
pack_size <- read_env_integer("PACK_SIZE", 1L)

start_id <- task_id
ids <- start_id:(start_id + pack_size - 1)
ids <- ids[ids >= 1 & ids <= nsim]
if (length(ids) == 0) stop("No valid ids to run for this task.")

# -------- Simulation settings (single combination from env) --------
beta2 <- read_env_numeric("BETA2", 0.5)
kappa <- read_env_numeric("KAPPA", 2)
c1 <- read_env_numeric("C1", 1)
cspline <- read_env_numeric("C", 1)
covariate_shift <- read_env_integer("SHIFT", 1L) == 1
L <- read_env_integer("LVAL", 5L)

n_target <- read_env_integer("N_TARGET", 200L)
n_source <- read_env_integer("N_SOURCE", 1000L)

# -------- Output paths --------
if (n_target == 500 && n_source == 2000) {
  base_out <- path.expand(Sys.getenv("OUT_BASE", unset = file.path(getwd(), "results_largen")))
} else {
  base_out <- path.expand(Sys.getenv("OUT_BASE", unset = file.path(getwd(), "results")))
}
out_dir <- file.path(
  base_out,
  sprintf("n%dN%d", n_target, n_source),
  sprintf("beta%.1f_kappa%s_L%d_c%.2f_C%.2f_shift%s",
          beta2, format_kappa_tag(kappa), L, c1, cspline, ifelse(covariate_shift, "TRUE", "FALSE"))
)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

simulate_reference_target_event_times <- function(L = 5, n_ref = 200000L,
                                                  seed = 20260421L) {
  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_seed_exists) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (old_seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(as.integer(seed + L))
  if (L == 5) {
    x_t <- cbind(
      runif(n_ref),
      rbinom(n_ref, 1, 0.5),
      runif(n_ref),
      runif(n_ref),
      runif(n_ref)
    )
    beta_t <- c(-0.5, 0.5, 0.2, 0.1, 0.1)
  } else {
    Sigma <- outer(seq_len(L), seq_len(L), function(i, j) 0.5^abs(i - j))
    x_t <- MASS::mvrnorm(n_ref, mu = rep(0, L), Sigma = Sigma)
    beta_t <- c(rep(0.5, 5), rep(0, L - 5))
  }
  linpred_t <- x_t %*% beta_t
  scale_t <- sqrt(2 / 2) * exp(-linpred_t / 2)
  as.numeric(rweibull(n_ref, shape = 2, scale = scale_t))
}

reference_target_quant_times <- local({
  cache <- new.env(parent = emptyenv())

  function(L = 5, probs = seq(0.1, 0.9, by = 0.1),
           n_ref = as.integer(Sys.getenv("REF_TARGET_N", "200000"))) {
    key <- sprintf("L%d_n%d", as.integer(L), as.integer(n_ref))
    if (!exists(key, envir = cache, inherits = FALSE)) {
      assign(
        key,
        simulate_reference_target_event_times(L = L, n_ref = n_ref),
        envir = cache
      )
    }
    ref_times <- get(key, envir = cache, inherits = FALSE)
    as.numeric(stats::quantile(ref_times, probs = probs, type = 7, names = FALSE))
  }
})

run_one_setting <- function(beta2, kappa, c1, L, covariate_shift, seed_id) {
  D <- simulate_once_all(
    n_target = n_target, n_source = n_source,
    beta2_source = beta2, kappa_source = kappa,
    c1 = c1, covariate_shift = covariate_shift, L = L
  )

  tuning_file <- default_sievetl_tuning_path(
    beta2 = beta2, kappa = kappa, c1 = c1,
    n0 = n_target, nA = n_source, L = L,
    shift = covariate_shift, base_dir = getwd()
  )
  tuning_row <- read_sievetl_tuning_file(tuning_file)
  cspline_use <- as.numeric(tuning_row$c_hat[1])
  cz_fixed <- as.numeric(tuning_row$c_zeta_hat[1])
  ce_fixed <- as.numeric(tuning_row$c_eta_hat[1])

  result <- sievetl_approx(
    target = D$target, source = D$source,
    L = L,
    c_zeta = cz_fixed, c_eta = ce_fixed,
    n0_penalty = n_target,
    c = cspline_use
  )

  pooled_stime <- c(D$target$stime, D$source[[1]]$stime)
  pooled_type  <- c(D$target$type,  D$source[[1]]$type)
  pooled_event_times <- pooled_stime[pooled_type == 1]
  N <- length(unique(pooled_event_times))
  J <- floor(cspline_use * N^(1/3))
  if (J < 1) J <- 0
  knots <- if (J > 0) {
    stats::quantile(pooled_stime, probs = (1:J) / (J + 1),
                    type = 7, names = FALSE)
  } else {
    NULL
  }
  bkn <- range(pooled_stime)
  deg <- 3
  B0 <- splines::bs(pooled_stime, degree = deg, knots = knots,
                    Boundary.knots = bkn, intercept = TRUE)
  p_hat <- ncol(B0)
  d_n <- L + p_hat
  by_column_lambda_file <- file.path(
    default_pilot_design_dir(
      beta2 = beta2, kappa = kappa, c1 = c1,
      n0 = n_target, nA = n_source, L = L,
      shift = covariate_shift, base_dir = getwd()
    ),
    "final_selected_c_by_column.csv"
  )
  g_funcs <- lapply(seq_len(ncol(B0)), function(j) {
    function(t) {
      Bt <- splines::bs(t, degree = deg, knots = knots,
                        Boundary.knots = bkn, intercept = TRUE)
      as.numeric(Bt[, j])
    }
  })

  time_grid <- seq(min(D$target$stime, na.rm = TRUE),
                   max(D$target$stime, na.rm = TRUE), length.out = 500)

  # New joint sieve transfer debiasing method.
  # Parameter ordering is fixed as theta = (gamma, beta) everywhere below.
  theta_hat <- c(result$gamma.hat, result$beta.hat)

  H_target <- compute_target_hessian_theta(
    D$target, result$gamma.hat, result$beta.hat,
    g_funcs = g_funcs, time_grid = time_grid
  )
  pooled_time_grid <- seq(min(pooled_stime, na.rm = TRUE),
                          max(pooled_stime, na.rm = TRUE), length.out = 500)
  H_pooled <- compute_pooled_hessian_theta(
    D$target, D$source, result$gamma.hat, result$beta.hat,
    g_funcs = g_funcs, time_grid = pooled_time_grid
  )
  score_vec <- compute_target_score_theta(
    D$target, result$gamma.hat, result$beta.hat,
    g_funcs = g_funcs, time_grid = time_grid
  )
  score_subject_mat <- compute_target_score_theta_subjects(
    D$target, result$gamma.hat, result$beta.hat,
    g_funcs = g_funcs, time_grid = time_grid
  )

  by_column_lambda <- read_by_column_lambda_file(by_column_lambda_file, d_n)
  Omega_hat_by_column <- estimate_transfer_inverse_hessian_by_column(
    H_pooled = H_pooled,
    H_target = H_target,
    lambda_delta_vec = by_column_lambda$lambda
  )
  hessian_inverse_product_by_column <- H_target %*% Omega_hat_by_column
  cat("Using by-column lambda file for seed", seed_id, ":\n")
  cat(by_column_lambda_file, "\n")
  cat("H_target %*% Omega_hat_by_column for seed", seed_id, ":\n")
  print(round(hessian_inverse_product_by_column, 4))
  theta_db_by_column <- as.vector(theta_hat - Omega_hat_by_column %*% score_vec)
  gamma_db_by_column <- theta_db_by_column[1:p_hat]
  beta_db_by_column  <- theta_db_by_column[(p_hat + 1):(p_hat + L)]
  if_theta_db_by_column_mat <- -score_subject_mat %*% t(Omega_hat_by_column)
  mean_if_theta_db_by_column <- colMeans(if_theta_db_by_column_mat)
  theta_db_col_if <- as.vector(theta_hat + mean_if_theta_db_by_column)
  gamma_db_col_if <- theta_db_col_if[1:p_hat]
  beta_db_col_if <- theta_db_col_if[(p_hat + 1):(p_hat + L)]
  theta_db_col_plugin <- theta_db_by_column
  gamma_db_col_plugin <- gamma_db_by_column
  beta_db_col_plugin <- beta_db_by_column

  Sigma_theta_hat_omega <- (Omega_hat_by_column + t(Omega_hat_by_column)) / 2
  Sigma_gamma_hat_omega <- Sigma_theta_hat_omega[1:p_hat, 1:p_hat, drop = FALSE]
  Sigma_beta_hat_omega <- Sigma_theta_hat_omega[(p_hat + 1):(p_hat + L), (p_hat + 1):(p_hat + L), drop = FALSE]
  se_beta_omega <- sqrt(pmax(diag(Sigma_beta_hat_omega) / n0, 0))

  Omega_target <- compute_target_inverse_hessian(H_target)
  hessian_inverse_target_product <- H_target %*% Omega_target
  cat("H_target %*% Omega_target for seed", seed_id, ":\n")
  print(round(hessian_inverse_target_product, 4))
  theta_os <- as.vector(theta_hat - Omega_target %*% score_vec)
  gamma_os <- theta_os[1:p_hat]
  beta_os  <- theta_os[(p_hat + 1):(p_hat + L)]
  if_theta_os_mat <- -score_subject_mat %*% t(Omega_target)
  mean_if_theta_os <- colMeans(if_theta_os_mat)
  theta_os_if <- as.vector(theta_hat + mean_if_theta_os)
  gamma_os_if <- theta_os_if[1:p_hat]
  beta_os_if <- theta_os_if[(p_hat + 1):(p_hat + L)]
  theta_os_plugin <- theta_os
  gamma_os_plugin <- gamma_os
  beta_os_plugin <- beta_os

  # TODO: add standard-error / covariance estimation for the transfer inverse-Hessian
  # debiased estimator only after the theory is finalized. The retired plugin-IF
  # variance formulas are intentionally not reused here.

  pooled_observed_times <- sort(unique(pooled_stime[is.finite(pooled_stime)]))
  if (length(pooled_observed_times) < 2) {
    pooled_observed_times <- sort(unique(c(bkn[1], bkn[2])))
  }

  t_q <- reference_target_quant_times(L = L, probs = seq(0, 1, by = 0.1))
  curve_eval_times <- sort(unique(c(pooled_observed_times, t_q[is.finite(t_q)])))

  bh_combined_curve <- make_baseline_from_gamma(result$gamma_hat_A, g_funcs, curve_eval_times)
  bh_sieveTL_curve  <- make_baseline_from_gamma(result$gamma.hat,   g_funcs, curve_eval_times)
  bh_db_col_full_curve <- make_baseline_from_gamma(gamma_db_by_column, g_funcs, curve_eval_times)
  bh_os_full_curve  <- make_baseline_from_gamma(gamma_os, g_funcs, curve_eval_times)
  # Plug-in one-step cumulative-hazard estimators rebuilt from the updated gamma.
  cumhaz_db_col_plugin_curve <- bh_db_col_full_curve[, c("time", "H0")]
  cumhaz_os_plugin_curve <- bh_os_full_curve[, c("time", "H0")]
  # IF-based one-step cumulative-hazard estimators on the same pooled-observed grid.
  cumhaz_db_col_curve <- compute_cumhaz_if_curve(result$gamma.hat, if_theta_db_by_column_mat, g_funcs, curve_eval_times)
  cumhaz_os_curve <- compute_cumhaz_if_curve(result$gamma.hat, if_theta_os_mat, g_funcs, curve_eval_times)
  baseline_db_col_curve <- bh_db_col_full_curve[, c("time", "logh0", "h0")]
  baseline_os_curve <- bh_os_full_curve[, c("time", "logh0", "h0")]

  logh0_combined_at <- approx(bh_combined_curve$time, bh_combined_curve$logh0,
                              xout = t_q, rule = 2, ties = "ordered")$y
  h0_combined_at <- approx(bh_combined_curve$time, bh_combined_curve$h0,
                           xout = t_q, rule = 2, ties = "ordered")$y
  H0_combined_at <- approx(bh_combined_curve$time, bh_combined_curve$H0,
                           xout = t_q, rule = 2, ties = "ordered")$y

  logh0_sieveTL_at <- approx(bh_sieveTL_curve$time, bh_sieveTL_curve$logh0,
                             xout = t_q, rule = 2, ties = "ordered")$y
  h0_sieveTL_at <- approx(bh_sieveTL_curve$time, bh_sieveTL_curve$h0,
                          xout = t_q, rule = 2, ties = "ordered")$y
  H0_sieveTL_at <- approx(bh_sieveTL_curve$time, bh_sieveTL_curve$H0,
                          xout = t_q, rule = 2, ties = "ordered")$y

  logh0_db_col_at <- approx(bh_db_col_full_curve$time, bh_db_col_full_curve$logh0,
                            xout = t_q, rule = 2, ties = "ordered")$y
  h0_db_col_at <- approx(bh_db_col_full_curve$time, bh_db_col_full_curve$h0,
                         xout = t_q, rule = 2, ties = "ordered")$y
  H0_db_col_plugin_at <- approx(cumhaz_db_col_plugin_curve$time, cumhaz_db_col_plugin_curve$H0,
                                xout = t_q, rule = 2, ties = "ordered")$y
  H0_db_col_at <- approx(cumhaz_db_col_curve$time, cumhaz_db_col_curve$H0,
                         xout = t_q, rule = 2, ties = "ordered")$y

  logh0_os_at <- approx(bh_os_full_curve$time, bh_os_full_curve$logh0,
                        xout = t_q, rule = 2, ties = "ordered")$y
  h0_os_at <- approx(bh_os_full_curve$time, bh_os_full_curve$h0,
                     xout = t_q, rule = 2, ties = "ordered")$y
  H0_os_plugin_at <- approx(cumhaz_os_plugin_curve$time, cumhaz_os_plugin_curve$H0,
                            xout = t_q, rule = 2, ties = "ordered")$y
  H0_os_at <- approx(cumhaz_os_curve$time, cumhaz_os_curve$H0,
                     xout = t_q, rule = 2, ties = "ordered")$y

  a_grid_obj <- compute_a_matrix(result$gamma.hat, g_funcs, curve_eval_times)
  a_quant_mat_omega <- interpolate_a_matrix(a_grid_obj, t_q)
  a_grid_mat_omega <- a_grid_obj$a_mat
  var_H0_db_col_plugin_omega_at <- rowSums((a_quant_mat_omega %*% Sigma_gamma_hat_omega) * a_quant_mat_omega) / n0
  var_H0_db_col_plugin_omega_at <- pmax(var_H0_db_col_plugin_omega_at, 0)
  se_H0_db_col_plugin_omega_at <- sqrt(var_H0_db_col_plugin_omega_at)
  var_H0_db_col_plugin_curve_omega <- rowSums((a_grid_mat_omega %*% Sigma_gamma_hat_omega) * a_grid_mat_omega) / n0
  var_H0_db_col_plugin_curve_omega <- pmax(var_H0_db_col_plugin_curve_omega, 0)

  risk_sieveTL <- exp(sum(result$beta.hat * x0))
  S_sieveTL_at <- exp(-H0_sieveTL_at * risk_sieveTL)
  grad_S_quant <- lapply(seq_along(H0_db_col_plugin_at), function(j) {
    c(
      -S_sieveTL_at[j] * risk_sieveTL * a_quant_mat_omega[j, ],
      -S_sieveTL_at[j] * risk_sieveTL * H0_sieveTL_at[j] * x0
    )
  })
  var_S_db_col_plugin_omega_at <- vapply(grad_S_quant, function(g) {
    as.numeric(t(g) %*% Sigma_theta_hat_omega %*% g)
  }, numeric(1)) / n0
  var_S_db_col_plugin_omega_at <- pmax(var_S_db_col_plugin_omega_at, 0)
  se_S_db_col_plugin_omega_at <- sqrt(var_S_db_col_plugin_omega_at)

  list(
    sievetl = list(
      result = result,
      time_quantiles = t_q,
      logh0_combined_at = logh0_combined_at,
      h0_combined_at = h0_combined_at,
      H0_combined_at = H0_combined_at,
      logh0_sieveTL_at = logh0_sieveTL_at,
      h0_sieveTL_at = h0_sieveTL_at,
      H0_sieveTL_at = H0_sieveTL_at,
      logh0_db_col_at = logh0_db_col_at,
      h0_db_col_at = h0_db_col_at,
      H0_db_col_plugin_at = H0_db_col_plugin_at,
      H0_db_col_at = H0_db_col_at,
      logh0_os_at = logh0_os_at,
      h0_os_at = h0_os_at,
      H0_os_plugin_at = H0_os_plugin_at,
      H0_os_at = H0_os_at,
      logh0_combined_curve = bh_combined_curve,
      h0_combined_curve = bh_combined_curve,
      H0_combined_curve = bh_combined_curve,
      logh0_sieveTL_curve = bh_sieveTL_curve,
      h0_sieveTL_curve = bh_sieveTL_curve,
      H0_sieveTL_curve = bh_sieveTL_curve,
      baseline_db_col_curve = baseline_db_col_curve,
      cumhaz_db_col_plugin_curve = cumhaz_db_col_plugin_curve,
      cumhaz_db_col_curve = cumhaz_db_col_curve,
      baseline_os_curve = baseline_os_curve,
      cumhaz_os_plugin_curve = cumhaz_os_plugin_curve,
      cumhaz_os_curve = cumhaz_os_curve,
      by_column_lambda_file = by_column_lambda_file,
      c_db_by_column = by_column_lambda$c,
      lambda_db_by_column = by_column_lambda$lambda,
      d_n = d_n,
      hessian_inverse_product_by_column = hessian_inverse_product_by_column,
      hessian_inverse_target_product = hessian_inverse_target_product,
      theta_db_col_if = as.numeric(theta_db_col_if),
      gamma_db_col_if = as.numeric(gamma_db_col_if),
      beta_db_col_if = as.numeric(beta_db_col_if),
      theta_db_col_plugin = as.numeric(theta_db_col_plugin),
      gamma_db_col_plugin = as.numeric(gamma_db_col_plugin),
      beta_db_col_plugin = as.numeric(beta_db_col_plugin),
      Sigma_theta_hat_omega = Sigma_theta_hat_omega,
      Sigma_gamma_hat_omega = Sigma_gamma_hat_omega,
      Sigma_beta_hat_omega = Sigma_beta_hat_omega,
      se_beta_omega = as.numeric(se_beta_omega),
      var_H0_db_col_plugin_omega_at = as.numeric(var_H0_db_col_plugin_omega_at),
      se_H0_db_col_plugin_omega_at = as.numeric(se_H0_db_col_plugin_omega_at),
      var_H0_db_col_plugin_curve_omega = as.numeric(var_H0_db_col_plugin_curve_omega),
      var_S_db_col_plugin_omega_at = as.numeric(var_S_db_col_plugin_omega_at),
      se_S_db_col_plugin_omega_at = as.numeric(se_S_db_col_plugin_omega_at),
      gamma_db_by_column = as.numeric(gamma_db_by_column),
      beta_db_by_column = as.numeric(beta_db_by_column),
      theta_os_if = as.numeric(theta_os_if),
      gamma_os_if = as.numeric(gamma_os_if),
      beta_os_if = as.numeric(beta_os_if),
      theta_os_plugin = as.numeric(theta_os_plugin),
      gamma_os_plugin = as.numeric(gamma_os_plugin),
      beta_os_plugin = as.numeric(beta_os_plugin),
      gamma_os = as.numeric(gamma_os),
      beta_os = as.numeric(beta_os)
    ),
    selection_fixed = list(
      c = cspline_use,
      c_zeta = cz_fixed,
      c_eta = ce_fixed,
      lambda_zeta = result$lambda_zeta,
      lambda_eta = result$lambda_eta,
      d_penalty = result$d_penalty,
      n0_penalty = result$n0_penalty,
      tuning_file = tuning_file
    )
  )
}

run_seed_batch <- function(ids, beta2, kappa, c1, L, covariate_shift, out_dir) {
  for (id in ids) {
    set.seed(100000 + id)
    results <- list()

    key <- sprintf("beta%.1f_kappa%s_L%d_c%.2f_C%.2f_shift%s",
                   beta2, format_kappa_tag(kappa), L, c1, cspline, ifelse(covariate_shift, "TRUE", "FALSE"))
    res <- tryCatch(
      run_one_setting(beta2, kappa, c1, L, covariate_shift, id),
      error = function(e) list(error = conditionMessage(e))
    )
    results[[key]] <- res

    out_file <- file.path(out_dir, sprintf("seed_%04d.rds", id))
    saveRDS(list(
      seed = id,
      settings = results,
      params = list(beta2 = beta2, kappa = kappa, c1 = c1, L = L, shift = covariate_shift)
    ), out_file)
  }
}

if (sys.nframe() == 0) {
  run_seed_batch(ids, beta2, kappa, c1, L, covariate_shift, out_dir)
}
