# Pointwise Monte Carlo coverage for the SieveTL transfer estimator.
#
# This file now constructs influence functions from the joint finite-dimensional
# sieve parameter theta = (gamma, beta) using the transfer inverse Hessian.
# The cumulative baseline hazard is treated only as a smooth plug-in functional
# of gamma_os. There is no separate one-step correction solved directly for the
# cumulative hazard.

suppressPackageStartupMessages({
  library(MASS)
  library(pracma)
})

# ------------------------ simulation + fit helpers ------------------------

transfer_loss_fast <- function(theta, p, L,
                               all.x, all.stime, all.type,
                               G_mat, time_grid, g_time_mat) {
  gamma <- theta[1:p]
  beta <- theta[(p + 1):(L + p)]

  xb <- as.vector(all.x %*% beta)
  g_gamma <- if (p == 1) as.vector(G_mat * gamma) else as.vector(G_mat %*% gamma)
  term1 <- all.type * (xb + g_gamma)

  exp_gamma_g <- exp(g_time_mat %*% gamma)
  int_vec <- vapply(all.stime, function(Ti) {
    idx <- which(time_grid <= Ti)
    if (length(idx) < 2L) return(0)
    trapz(time_grid[idx], exp_gamma_g[idx])
  }, numeric(1))

  term2 <- exp(xb) * int_vec
  -mean(term1 - term2, na.rm = TRUE)
}

transfer_grad_fast <- function(theta, p, L,
                               all.x, all.stime, all.type,
                               G_mat, time_grid, g_time_mat) {
  gamma <- theta[1:p]
  beta <- theta[(p + 1):(L + p)]

  xb <- as.vector(all.x %*% beta)
  exp_xb <- exp(xb)
  eg_time <- as.vector(g_time_mat %*% gamma)
  exp_eg <- exp(eg_time)

  int_vec <- vapply(all.stime, function(Ti) {
    idx <- which(time_grid <= Ti)
    if (length(idx) < 2L) return(0)
    pracma::trapz(time_grid[idx], exp_eg[idx])
  }, numeric(1))

  term_mat <- matrix(0, nrow = length(all.stime), ncol = p)
  for (i in seq_along(all.stime)) {
    idx <- which(time_grid <= all.stime[i])
    if (length(idx) < 2L) next
    for (j in seq_len(p)) {
      term_mat[i, j] <- pracma::trapz(time_grid[idx], g_time_mat[idx, j] * exp_eg[idx])
    }
  }

  grad_gamma <- -colMeans(G_mat * all.type, na.rm = TRUE) +
    colMeans(exp_xb * term_mat, na.rm = TRUE)
  grad_beta <- -colMeans(all.x * all.type, na.rm = TRUE) +
    colMeans(all.x * (exp_xb * int_vec), na.rm = TRUE)

  c(grad_gamma, grad_beta)
}

.smooth_abs <- function(u, eps) sqrt(u * u + eps * eps) - eps
.smooth_abs_grad <- function(u, eps) u / sqrt(u * u + eps * eps)

debias_loss_fast <- function(delta, lambda_zeta, lambda_eta, L, p,
                             target, G_mat_target, time_grid_target,
                             g_time_mat_target, gamma_hat_A, beta_hat_A,
                             eps = 1e-4) {
  zeta <- delta[1:p]
  eta <- delta[(p + 1):(L + p)]

  gamma <- gamma_hat_A + zeta
  beta <- beta_hat_A + eta

  xb <- as.vector(target$x %*% beta)
  g_gamma <- if (p == 1) as.vector(G_mat_target * gamma) else as.vector(G_mat_target %*% gamma)
  term1 <- target$type * (xb + g_gamma)

  eg_time <- as.vector(g_time_mat_target %*% gamma)
  exp_eg <- exp(eg_time)

  int_vec <- vapply(target$stime, function(Ti) {
    idx <- which(time_grid_target <= Ti)
    if (length(idx) < 2L) return(0)
    pracma::trapz(time_grid_target[idx], exp_eg[idx])
  }, numeric(1))

  term2 <- exp(xb) * int_vec
  nll <- -mean(term1 - term2, na.rm = TRUE)

  pen <- lambda_zeta * sum(.smooth_abs(zeta, eps)) +
    lambda_eta * sum(.smooth_abs(eta, eps))

  nll + pen
}

debias_grad_fast <- function(delta, lambda_zeta, lambda_eta, L, p,
                             target, G_mat_target, time_grid_target,
                             g_time_mat_target, gamma_hat_A, beta_hat_A,
                             eps = 1e-4) {
  zeta <- delta[1:p]
  eta <- delta[(p + 1):(L + p)]

  gamma <- gamma_hat_A + zeta
  beta <- beta_hat_A + eta

  xb <- as.vector(target$x %*% beta)
  exp_xb <- exp(xb)
  eg_time <- as.vector(g_time_mat_target %*% gamma)
  exp_eg <- exp(eg_time)

  int_vec <- vapply(target$stime, function(Ti) {
    idx <- which(time_grid_target <= Ti)
    if (length(idx) < 2L) return(0)
    pracma::trapz(time_grid_target[idx], exp_eg[idx])
  }, numeric(1))

  term_mat <- matrix(0, nrow = length(target$stime), ncol = p)
  for (i in seq_along(target$stime)) {
    idx <- which(time_grid_target <= target$stime[i])
    if (length(idx) < 2L) next
    for (j in seq_len(p)) {
      term_mat[i, j] <- pracma::trapz(time_grid_target[idx], g_time_mat_target[idx, j] * exp_eg[idx])
    }
  }

  grad_zeta <- -colMeans(G_mat_target * target$type, na.rm = TRUE) +
    colMeans(exp_xb * term_mat, na.rm = TRUE)
  grad_eta <- -colMeans(target$x * target$type, na.rm = TRUE) +
    colMeans(target$x * (exp_xb * int_vec), na.rm = TRUE)

  grad_zeta <- grad_zeta + lambda_zeta * .smooth_abs_grad(zeta, eps)
  grad_eta <- grad_eta + lambda_eta * .smooth_abs_grad(eta, eps)

  c(grad_zeta, grad_eta)
}

simulate_once_all <- function(n_target = 500, n_source = 2000,
                              beta2_source = 0.5, kappa_source = 1.5,
                              c1 = 1, covariate_shift = TRUE, L = 5) {
  if (!(L %in% c(5, 6))) stop("Only L = 5 or L = 50 is supported.")

  if (L == 5) {
    X1_t <- runif(n_target)
    X2_t <- rbinom(n_target, 1, 0.5)
    X3_t <- runif(n_target)
    X4_t <- runif(n_target)
    X5_t <- runif(n_target)
    x_t <- cbind(X1_t, X2_t, X3_t, X4_t, X5_t)
    beta_t <- c(-0.5, 0.5, 0.2, 0.1, 0.1)
  } else {
    Sigma <- outer(seq_len(L), seq_len(L), function(i, j) 0.5^abs(i - j))
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

  if (L == 5) {
    X1_s <- runif(n_source)
    X2_s <- rbinom(n_source, 1, 0.5)
    X3_s <- if (covariate_shift) rbeta(n_source, 1, 2) else runif(n_source)
    X4_s <- runif(n_source)
    X5_s <- runif(n_source)
    x_s <- cbind(X1_s, X2_s, X3_s, X4_s, X5_s)
    beta_s <- c(-0.5, beta2_source, 0.2, 0.1, 0.1)
  } else {
    Sigma <- outer(seq_len(L), seq_len(L), function(i, j) 0.5^abs(i - j))
    if (covariate_shift) {
      eps <- matrix(rnorm(L, sd = 0.3))
      Sigma <- Sigma + eps %*% t(eps)
    }
    x_s <- MASS::mvrnorm(n_source, mu = rep(0, L), Sigma = Sigma)
    beta_s <- c(0.5, beta2_source, 0.5, 0.5, 0.5, rep(0, L - 5))
  }
  linpred_s <- x_s %*% beta_s
  scale_s <- sqrt(2 / kappa_source) * exp(-linpred_s / 2)
  etime_s <- rweibull(n_source, shape = 2, scale = scale_s)
  ctime_s <- rexp(n_source, rate = 1 / c1)
  stime_s <- pmin(etime_s, ctime_s)
  type_s <- as.numeric(etime_s <= ctime_s)
  source <- list(list(x = x_s, stime = stime_s, type = type_s))

  list(target = target, source = source, beta_true = beta_t)
}

build_spline_helpers <- function(target, source, c_sieve = 1) {
  pooled_stime <- c(target$stime, unlist(lapply(source, function(s) s$stime)))
  pooled_type <- c(target$type, unlist(lapply(source, function(s) s$type)))
  pooled_event_times <- pooled_stime[pooled_type == 1]
  N <- length(unique(pooled_event_times))
  J <- floor(c_sieve * N^(1 / 3))
  if (J < 1) J <- 0
  knots <- if (J > 0) {
    stats::quantile(pooled_stime, probs = (1:J) / (J + 1), type = 7, names = FALSE)
  } else {
    NULL
  }
  bkn <- range(pooled_stime)
  deg <- 3
  B0 <- splines::bs(pooled_stime, degree = deg, knots = knots,
                    Boundary.knots = bkn, intercept = TRUE)
  p <- ncol(B0)
  g_funcs <- lapply(seq_len(p), function(j) {
    function(t) {
      Bt <- splines::bs(t, degree = deg, knots = knots,
                        Boundary.knots = bkn, intercept = TRUE)
      as.numeric(Bt[, j])
    }
  })
  list(g_funcs = g_funcs, p = p, pooled_stime = pooled_stime, pooled_type = pooled_type)
}

make_baseline_from_gamma <- function(gamma_hat, g_funcs, t_grid) {
  t_grid <- sort(unique(as.numeric(t_grid)))
  if (length(t_grid) < 2L) stop("Need at least two grid points to form baseline increments.")
  G <- sapply(g_funcs, function(gf) as.numeric(gf(t_grid)))
  if (is.null(dim(G))) G <- matrix(G, ncol = 1)
  logh0 <- drop(G %*% gamma_hat)
  h0 <- exp(logh0)
  dt <- diff(t_grid)
  H0 <- c(0, cumsum(dt * (head(h0, -1) + tail(h0, -1)) / 2))
  data.frame(time = t_grid, logh0 = logh0, h0 = h0, H0 = H0)
}

format_kappa_tag <- function(kappa) {
  if (isTRUE(all.equal(kappa, round(kappa), tolerance = 1e-8))) {
    return(sprintf("%d", as.integer(round(kappa))))
  }
  formatC(kappa, format = "f", digits = 1)
}

fit_sievetl_transfer <- function(target, source, L, c_sieve, lambda_zeta, lambda_eta) {
  helpers <- build_spline_helpers(target, source, c_sieve = c_sieve)
  g_funcs <- helpers$g_funcs
  p <- helpers$p

  all.x <- as.matrix(do.call(rbind, c(list(target$x), lapply(source, function(s) s$x))))
  all.stime <- c(target$stime, unlist(lapply(source, function(s) s$stime)))
  all.type <- c(target$type, unlist(lapply(source, function(s) s$type)))

  if (p == 1) {
    G_mat <- matrix(sapply(all.stime, function(ti) g_funcs[[1]](ti)), ncol = 1)
  } else {
    G_mat <- t(sapply(all.stime, function(ti) sapply(g_funcs, function(gj) gj(ti))))
  }
  time_grid <- seq(min(all.stime), max(all.stime), length.out = 500)
  g_time_mat <- if (p == 1) {
    matrix(g_funcs[[1]](time_grid), ncol = 1)
  } else {
    tmp <- sapply(g_funcs, function(gj) gj(time_grid))
    if (is.vector(tmp)) matrix(tmp, ncol = p) else tmp
  }

  theta0 <- rep(0, L + p)
  fit_stage1 <- nlminb(
    start = theta0,
    objective = function(th) transfer_loss_fast(th, p, L, all.x, all.stime, all.type, G_mat, time_grid, g_time_mat),
    gradient = function(th) transfer_grad_fast(th, p, L, all.x, all.stime, all.type, G_mat, time_grid, g_time_mat),
    control = list(iter.max = 5000, eval.max = 5000)
  )

  theta_hat_A <- fit_stage1$par
  gamma_hat_A <- theta_hat_A[1:p]
  beta_hat_A <- theta_hat_A[(p + 1):(L + p)]

  if (p == 1) {
    G_mat_target <- matrix(g_funcs[[1]](target$stime), ncol = 1)
  } else {
    G_mat_target <- t(sapply(target$stime, function(ti) sapply(g_funcs, function(gj) gj(ti))))
  }
  time_grid_target <- seq(min(target$stime), max(target$stime), length.out = 500)
  g_time_mat_target <- if (p == 1) {
    matrix(g_funcs[[1]](time_grid_target), ncol = 1)
  } else {
    tmp <- sapply(g_funcs, function(gj) gj(time_grid_target))
    if (is.vector(tmp)) matrix(tmp, ncol = p) else tmp
  }

  delta0 <- rep(0, L + p)
  fit_stage2 <- optim(
    par = delta0,
    fn = function(d) debias_loss_fast(d, lambda_zeta, lambda_eta, L, p,
                                      target, G_mat_target, time_grid_target,
                                      g_time_mat_target, gamma_hat_A, beta_hat_A),
    gr = function(d) debias_grad_fast(d, lambda_zeta, lambda_eta, L, p,
                                      target, G_mat_target, time_grid_target,
                                      g_time_mat_target, gamma_hat_A, beta_hat_A),
    method = "L-BFGS-B",
    control = list(maxit = 2000, factr = 1e-10, pgtol = 1e-6)
  )

  delta_hat <- fit_stage2$par
  zeta_hat <- delta_hat[1:p]
  eta_hat <- delta_hat[(p + 1):(L + p)]

  fit <- list(
    gamma_hat_A = gamma_hat_A,
    beta_hat_A = beta_hat_A,
    zeta_hat = zeta_hat,
    eta_hat = eta_hat,
    gamma.hat = gamma_hat_A + zeta_hat,
    beta.hat = beta_hat_A + eta_hat
  )

  list(fit = fit, g_funcs = g_funcs)
}

# ------------------------ plug-in IF coverage helpers ------------------------

compute_target_score_theta <- function(data, gamma_hat, beta_hat, g_funcs, time_grid) {
  n <- nrow(data$x)
  L <- length(beta_hat)
  p <- length(gamma_hat)
  d <- p + L  # order: theta = (gamma, beta)

  grad <- rep(0, d)

  for (i in seq_len(n)) {
    x_i <- data$x[i, ]
    t_i <- data$stime[i]
    delta_i <- data$type[i]

    if (t_i <= 0) next

    g_t <- sapply(g_funcs, function(f) f(t_i))
    z_t <- c(g_t, x_i)
    grad <- grad - delta_i * z_t

    idx <- which(time_grid <= t_i)
    if (length(idx) < 2L) next

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

  for (i in seq_len(n)) {
    x_i <- data$x[i, ]
    t_i <- data$stime[i]
    delta_i <- data$type[i]

    if (t_i <= 0) next

    g_t <- sapply(g_funcs, function(f) f(t_i))
    z_t <- c(g_t, x_i)
    score_i <- -delta_i * z_t

    idx <- which(time_grid <= t_i)
    if (length(idx) >= 2L) {
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

compute_pooled_score_theta <- function(target, source_list, gamma_hat, beta_hat, g_funcs, time_grid) {
  pooled_data <- combine_domains(target, source_list)
  compute_target_score_theta(pooled_data, gamma_hat, beta_hat, g_funcs, time_grid)
}

compute_pooled_score_theta_subjects <- function(target, source_list, gamma_hat, beta_hat, g_funcs, time_grid) {
  pooled_data <- combine_domains(target, source_list)
  compute_target_score_theta_subjects(pooled_data, gamma_hat, beta_hat, g_funcs, time_grid)
}

apply_one_step_update <- function(theta_hat, Omega_hat_col, score_vec) {
  as.vector(theta_hat - Omega_hat_col %*% score_vec)
}

compute_target_hessian_theta <- function(data_k, gamma_hat, beta_hat, g_funcs, time_grid) {
  n_k <- nrow(data_k$x)
  p <- length(gamma_hat)
  L <- length(beta_hat)
  d <- p + L

  hessian <- matrix(0, d, d)

  for (i in seq_len(n_k)) {
    x_i <- data_k$x[i, ]
    t_i <- data_k$stime[i]
    if (t_i <= 0) next

    t_grid_i <- time_grid[time_grid <= t_i]
    if (length(t_grid_i) < 2L) next

    G_mat <- t(sapply(g_funcs, function(f) f(t_grid_i)))
    m <- length(t_grid_i)

    eta_g_vec <- as.vector(gamma_hat %*% G_mat)
    exp_eta_g <- exp(eta_g_vec)
    eta_x <- sum(beta_hat * x_i)
    exp_eta_x <- exp(eta_x)

    integrand_array <- array(0, dim = c(d, d, m))
    for (j in seq_len(m)) {
      g_j <- G_mat[, j]
      z_j <- c(g_j, x_i)
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

default_by_column_lambda_path <- function(design_row, c_sieve, base_dir = getwd(), L = 5) {
  old_path <- file.path(
    default_pilot_design_dir(
      beta2 = design_row$beta2,
      kappa = design_row$kappa,
      c1 = design_row$c,
      n0 = design_row$n0,
      nA = design_row$nA,
      L = L,
      shift = isTRUE(design_row$shift),
      base_dir = base_dir
    ),
    "final_selected_c_by_column.csv"
  )
  if (file.exists(old_path)) {
    return(old_path)
  }

  new_path <- file.path(
    path.expand(base_dir),
    "pilot_parameter_by_column",
    sprintf("n%dN%d", design_row$n0, design_row$nA),
    sprintf(
      "beta%.1f_kappa%s_L%d_c%.2f_C%.2f_shift%s",
      design_row$beta2,
      format_kappa_tag(design_row$kappa),
      L,
      design_row$c,
      c_sieve,
      ifelse(isTRUE(design_row$shift), "TRUE", "FALSE")
    ),
    "final_selected_c_by_column.csv"
  )
  new_path
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

current_script_dir <- function() {
  for (i in rev(seq_len(sys.nframe()))) {
    ofile <- sys.frame(i)$ofile
    if (!is.null(ofile)) return(dirname(normalizePath(ofile)))
  }
  getwd()
}

load_exact_run_script_env <- function() {
  script_dir <- current_script_dir()
  run_script <- file.path(script_dir, "run_PLC_CR_spline_hoffman2.R")
  if (!file.exists(run_script)) {
    stop("Cannot find run_PLC_CR_spline_hoffman2.R next to coverage_plugin_wald_sievetl.R")
  }

  env_names <- c(
    "NSIM", "SGE_TASK_ID", "PACK_SIZE", "BETA2", "KAPPA", "C1", "C",
    "SHIFT", "LVAL", "N_TARGET", "N_SOURCE", "OUT_BASE"
  )
  old_env <- Sys.getenv(env_names, unset = NA_character_)
  on.exit({
    for (k in seq_along(env_names)) {
      if (is.na(old_env[k])) {
        Sys.unsetenv(env_names[k])
      } else {
        do.call(Sys.setenv, setNames(list(old_env[k]), env_names[k]))
      }
    }
  }, add = TRUE)

  Sys.setenv(
    NSIM = "1",
    SGE_TASK_ID = "1",
    PACK_SIZE = "1",
    BETA2 = "0.5",
    KAPPA = "2",
    C1 = "1",
    C = "1",
    SHIFT = "0",
    LVAL = "5",
    N_TARGET = "100",
    N_SOURCE = "1000",
    OUT_BASE = tempdir()
  )

  exact_env <- new.env(parent = globalenv())
  sys.source(run_script, envir = exact_env)
  exact_env
}

.exact_run_env <- load_exact_run_script_env()

for (fn_name in c(
  "transfer_loss_fast",
  "transfer_grad_fast",
  ".smooth_abs",
  ".smooth_abs_grad",
  "debias_loss_fast",
  "debias_grad_fast",
  "sievetl_approx",
  "simulate_once_all",
  "make_baseline_from_gamma",
  "compute_cumhaz_if_curve",
  "compute_target_score_theta",
  "compute_target_score_theta_subjects",
  "compute_target_hessian_theta",
  "combine_domains",
  "compute_pooled_hessian_theta",
  "fit_inverse_hessian_initializer_unpenalized",
  "fit_inverse_hessian_correction",
  "estimate_transfer_inverse_hessian_by_column",
  "compute_target_inverse_hessian",
  "read_by_column_lambda_file",
  "read_sievetl_tuning_file",
  "default_pilot_design_dir",
  "default_sievetl_tuning_path",
  "lambda_to_penalty_c",
  "reference_target_quant_times"
)) {
  assign(fn_name, get(fn_name, envir = .exact_run_env), envir = environment())
}

resolve_sievetl_tuning_row <- function(beta2, kappa, c1, n0, nA, L = 5,
                                       shift = FALSE, base_dir = getwd()) {
  tuning_file <- default_sievetl_tuning_path(beta2, kappa, c1, n0, nA, L, shift, base_dir)
  read_sievetl_tuning_file(tuning_file)
}

build_g_funcs_exact_from_run <- function(target, source, cspline) {
  pooled_stime <- c(target$stime, source[[1]]$stime)
  pooled_type  <- c(target$type,  source[[1]]$type)
  pooled_event_times <- pooled_stime[pooled_type == 1]
  N <- length(unique(pooled_event_times))
  J <- floor(cspline * N^(1/3))
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
  g_funcs <- lapply(seq_len(ncol(B0)), function(j) {
    function(t) {
      Bt <- splines::bs(t, degree = deg, knots = knots,
                        Boundary.knots = bkn, intercept = TRUE)
      as.numeric(Bt[, j])
    }
  })
  list(g_funcs = g_funcs, p_hat = p_hat)
}

true_beta_target <- function(L) {
  if (L == 5) return(c(-0.5, 0.5, 0.2, 0.1, 0.1))
  c(rep(0.5, 5), rep(0, L - 5))
}

compute_plugin_inference <- function(target, source, beta_hat, gamma_hat, g_funcs,
                                     lambda_delta_vec,
                                     delta_hat,
                                     beta_hat_se_base = beta_hat,
                                     gamma_hat_se_base = gamma_hat,
                                     x0 = NULL,
                                     quant_times = NULL,
                                     quant_probs = seq(0.1, 0.9, by = 0.1),
                                     integration_grid_size = 500,
                                     debug_q10 = FALSE) {
  x <- as.matrix(target$x)
  stime <- as.numeric(target$stime)
  delta <- as.numeric(target$type)
  n0 <- nrow(x)
  L <- ncol(x)
  p <- length(gamma_hat)
  d <- p + L
  if (is.null(x0)) {
    x0 <- rep(0.5, L)
  } else {
    x0 <- as.numeric(x0)
  }
  if (length(x0) != L) {
    stop("x0 must have length L.")
  }

  pooled_stime <- c(target$stime, unlist(lapply(source, `[[`, "stime")))
  pooled_type <- c(target$type, unlist(lapply(source, `[[`, "type")))
  pooled_observed_times <- sort(unique(as.numeric(pooled_stime[is.finite(pooled_stime)])))
  if (length(pooled_observed_times) == 0L) {
    pooled_observed_times <- sort(unique(stime[is.finite(stime)]))
  }
  if (length(pooled_observed_times) == 0L) {
    stop("No finite observed times available to build the inference grid.")
  }

  test_grid_bounds_10_90 <- reference_target_quant_times(L = L, probs = c(0.1, 0.9))
  test_grid_bounds_20_90 <- reference_target_quant_times(L = L, probs = c(0.2, 0.9))
  grid_start <- min(c(pooled_observed_times, test_grid_bounds_10_90), na.rm = TRUE)
  grid_end <- max(c(pooled_observed_times, test_grid_bounds_10_90), na.rm = TRUE)
  if (grid_end < grid_start) {
    tmp <- grid_start
    grid_start <- grid_end
    grid_end <- tmp
  }
  n_grid <- 500L
  scb_grid <- seq(grid_start, grid_end, length.out = n_grid)
  test_grid_10_90_start <- test_grid_bounds_10_90[1]
  test_grid_10_90_end <- test_grid_bounds_10_90[2]
  if (test_grid_10_90_end < test_grid_10_90_start) {
    test_grid_10_90_end <- test_grid_10_90_start
  }
  test_grid_10_90 <- seq(test_grid_10_90_start, test_grid_10_90_end, length.out = 50L)
  test_grid_20_90_start <- test_grid_bounds_20_90[1]
  test_grid_20_90_end <- test_grid_bounds_20_90[2]
  if (test_grid_20_90_end < test_grid_20_90_start) {
    test_grid_20_90_end <- test_grid_20_90_start
  }
  test_grid_20_90 <- seq(test_grid_20_90_start, test_grid_20_90_end, length.out = 50L)
  time_grid_target <- seq(min(stime, na.rm = TRUE), max(stime, na.rm = TRUE), length.out = max(as.integer(integration_grid_size), 2L))
  if (length(pooled_observed_times) < 2L) {
    pooled_observed_times <- sort(unique(c(grid_start, grid_end)))
  }
  pooled_time_grid <- seq(min(pooled_stime, na.rm = TRUE), max(pooled_stime, na.rm = TRUE), length.out = max(as.integer(integration_grid_size), 2L))

  H_target <- compute_target_hessian_theta(target, gamma_hat, beta_hat, g_funcs, time_grid_target)
  H_pooled <- compute_pooled_hessian_theta(target, source, gamma_hat, beta_hat, g_funcs, pooled_time_grid)
  Omega_hat_by_column <- estimate_transfer_inverse_hessian_by_column(H_pooled, H_target, lambda_delta_vec)
  Omega_hat_by_column_zero <- estimate_transfer_inverse_hessian_by_column(
    H_pooled, H_target, rep(0, d)
  )

  theta_hat <- c(gamma_hat, beta_hat)
  if (!is.numeric(delta_hat) || length(delta_hat) != d) {
    stop("delta_hat must be a numeric vector of length p + L.")
  }
  score_vec <- compute_target_score_theta(target, gamma_hat, beta_hat, g_funcs, time_grid_target)
  score_subject_mat <- compute_target_score_theta_subjects(target, gamma_hat, beta_hat, g_funcs, time_grid_target)
  # Match run_PLC_CR_spline_hoffman2.R exactly for the active transfer-by-column path.
  if_theta_db_by_column_mat <- -score_subject_mat %*% t(Omega_hat_by_column)
  IF_theta_hat_mat <- if_theta_db_by_column_mat
  IF_gamma_hat_mat <- IF_theta_hat_mat[, seq_len(p), drop = FALSE]
  IF_beta_hat_mat <- IF_theta_hat_mat[, (p + 1):(p + L), drop = FALSE]

  # Second variance track based directly on the transfer inverse Hessian:
  # Varhat(theta_db) = Omega_hat / n0.
  Sigma_theta_hat_omega <- (Omega_hat_by_column + t(Omega_hat_by_column)) / 2
  Sigma_gamma_hat_omega <- Sigma_theta_hat_omega[seq_len(p), seq_len(p), drop = FALSE]
  Sigma_beta_hat_omega <- Sigma_theta_hat_omega[(p + 1):(p + L), (p + 1):(p + L), drop = FALSE]
  Sigma_theta_hat_omega_zero <- (Omega_hat_by_column_zero + t(Omega_hat_by_column_zero)) / 2
  Sigma_gamma_hat_omega_zero <- Sigma_theta_hat_omega_zero[seq_len(p), seq_len(p), drop = FALSE]
  Sigma_beta_hat_omega_zero <- Sigma_theta_hat_omega_zero[(p + 1):(p + L), (p + 1):(p + L), drop = FALSE]

  theta_db_by_column <- apply_one_step_update(theta_hat, Omega_hat_by_column, score_vec)
  gamma_db_by_column <- theta_db_by_column[seq_len(p)]
  beta_db_by_column <- theta_db_by_column[(p + 1):(p + L)]

  theta_os <- theta_db_by_column
  gamma_os <- gamma_db_by_column
  beta_os <- beta_db_by_column

  beta_hat_os <- beta_os
  avar_beta <- crossprod(IF_beta_hat_mat) / n0
  avar_beta <- (avar_beta + t(avar_beta)) / 2
  se_beta <- sqrt(pmax(diag(avar_beta) / n0, 0))
  se_beta_omega <- sqrt(pmax(diag(Sigma_beta_hat_omega) / n0, 0))
  se_beta_omega_zero <- sqrt(pmax(diag(Sigma_beta_hat_omega_zero) / n0, 0))

  if (is.null(quant_times)) {
    quant_times <- reference_target_quant_times(L = L, probs = quant_probs)
  } else {
    quant_times <- as.numeric(quant_times)
  }
  quant_labels <- paste0("Lambda_q", sprintf("%02d", round(100 * quant_probs)))
  curve_eval_times <- sort(unique(c(pooled_observed_times, quant_times[is.finite(quant_times)])))

  cumhaz_tl_curve <- make_baseline_from_gamma(gamma_hat, g_funcs, curve_eval_times)
  Lambda0_tl_quant <- stats::approx(cumhaz_tl_curve$time, cumhaz_tl_curve$H0,
                                    xout = quant_times, rule = 2, ties = "ordered")$y
  Lambda0_tl_grid <- stats::approx(cumhaz_tl_curve$time, cumhaz_tl_curve$H0,
                                   xout = scb_grid, rule = 2, ties = "ordered")$y
  Lambda0_tl_test_grid_10_90 <- stats::approx(cumhaz_tl_curve$time, cumhaz_tl_curve$H0,
                                              xout = test_grid_10_90, rule = 2, ties = "ordered")$y
  Lambda0_tl_test_grid_20_90 <- stats::approx(cumhaz_tl_curve$time, cumhaz_tl_curve$H0,
                                              xout = test_grid_20_90, rule = 2, ties = "ordered")$y

  cumhaz_os_linearized_curve <- compute_cumhaz_if_curve(gamma_hat_se_base, if_theta_db_by_column_mat, g_funcs, curve_eval_times)
  cumhaz_os_curve <- cumhaz_os_linearized_curve
  Lambda0_os_quant <- stats::approx(cumhaz_os_curve$time, cumhaz_os_curve$H0,
                                    xout = quant_times, rule = 2, ties = "ordered")$y
  Lambda0_os_grid <- stats::approx(cumhaz_os_curve$time, cumhaz_os_curve$H0,
                                   xout = scb_grid, rule = 2, ties = "ordered")$y
  Lambda0_os_test_grid_10_90 <- stats::approx(cumhaz_os_curve$time, cumhaz_os_curve$H0,
                                              xout = test_grid_10_90, rule = 2, ties = "ordered")$y
  Lambda0_os_test_grid_20_90 <- stats::approx(cumhaz_os_curve$time, cumhaz_os_curve$H0,
                                              xout = test_grid_20_90, rule = 2, ties = "ordered")$y

  a_pooled_obj <- compute_a_matrix(gamma_hat, g_funcs, curve_eval_times)
  IF_Lambda_hat_pooled <- IF_gamma_hat_mat %*% t(a_pooled_obj$a_mat)
  IF_Lambda_hat_mat <- t(vapply(seq_len(nrow(IF_Lambda_hat_pooled)), function(i) {
    stats::approx(a_pooled_obj$time, IF_Lambda_hat_pooled[i, ], xout = quant_times, rule = 2, ties = "ordered")$y
  }, numeric(length(quant_times))))
  a_grid_obj <- compute_a_matrix(gamma_hat, g_funcs, scb_grid)
  IF_Lambda_hat_grid <- IF_gamma_hat_mat %*% t(a_grid_obj$a_mat)
  a_pooled_obj_os <- compute_a_matrix(gamma_os, g_funcs, curve_eval_times)
  a_grid_obj_os <- compute_a_matrix(gamma_os, g_funcs, scb_grid)
  IF_Lambda_os_pooled <- IF_gamma_hat_mat %*% t(a_pooled_obj_os$a_mat)
  IF_Lambda_os_mat <- t(vapply(seq_len(nrow(IF_Lambda_os_pooled)), function(i) {
    stats::approx(a_pooled_obj_os$time, IF_Lambda_os_pooled[i, ], xout = quant_times, rule = 2, ties = "ordered")$y
  }, numeric(length(quant_times))))
  IF_Lambda_os_grid <- IF_gamma_hat_mat %*% t(a_grid_obj_os$a_mat)

  a_quant_mat_omega <- interpolate_a_matrix(a_pooled_obj, quant_times)
  a_grid_mat_omega <- interpolate_a_matrix(a_pooled_obj, scb_grid)
  a_quant_mat_omega_os <- interpolate_a_matrix(a_pooled_obj_os, quant_times)
  a_grid_mat_omega_os <- interpolate_a_matrix(a_pooled_obj_os, scb_grid)
  var_hat_Lambda_omega <- rowSums((a_quant_mat_omega %*% Sigma_gamma_hat_omega) * a_quant_mat_omega) / n0
  var_hat_Lambda_omega <- pmax(var_hat_Lambda_omega, 0)
  se_Lambda_omega <- sqrt(var_hat_Lambda_omega)
  var_hat_Lambda_grid_omega <- rowSums((a_grid_mat_omega %*% Sigma_gamma_hat_omega) * a_grid_mat_omega) / n0
  var_hat_Lambda_grid_omega <- pmax(var_hat_Lambda_grid_omega, 0)
  var_hat_Lambda_omega_zero <- rowSums((a_quant_mat_omega %*% Sigma_gamma_hat_omega_zero) * a_quant_mat_omega) / n0
  var_hat_Lambda_omega_zero <- pmax(var_hat_Lambda_omega_zero, 0)
  se_Lambda_omega_zero <- sqrt(var_hat_Lambda_omega_zero)
  var_hat_Lambda_grid_omega_zero <- rowSums((a_grid_mat_omega %*% Sigma_gamma_hat_omega_zero) * a_grid_mat_omega) / n0
  var_hat_Lambda_grid_omega_zero <- pmax(var_hat_Lambda_grid_omega_zero, 0)
  var_hat_Lambda_omega_os <- rowSums((a_quant_mat_omega_os %*% Sigma_gamma_hat_omega) * a_quant_mat_omega_os) / n0
  var_hat_Lambda_omega_os <- pmax(var_hat_Lambda_omega_os, 0)
  se_Lambda_omega_os <- sqrt(var_hat_Lambda_omega_os)
  var_hat_Lambda_grid_omega_os <- rowSums((a_grid_mat_omega_os %*% Sigma_gamma_hat_omega) * a_grid_mat_omega_os) / n0
  var_hat_Lambda_grid_omega_os <- pmax(var_hat_Lambda_grid_omega_os, 0)
  var_hat_Lambda_omega_zero_os <- rowSums((a_quant_mat_omega_os %*% Sigma_gamma_hat_omega_zero) * a_quant_mat_omega_os) / n0
  var_hat_Lambda_omega_zero_os <- pmax(var_hat_Lambda_omega_zero_os, 0)
  se_Lambda_omega_zero_os <- sqrt(var_hat_Lambda_omega_zero_os)
  var_hat_Lambda_grid_omega_zero_os <- rowSums((a_grid_mat_omega_os %*% Sigma_gamma_hat_omega_zero) * a_grid_mat_omega_os) / n0
  var_hat_Lambda_grid_omega_zero_os <- pmax(var_hat_Lambda_grid_omega_zero_os, 0)

  risk_os <- exp(sum(beta_os * x0))
  risk_tl <- exp(sum(beta_hat * x0))
  R_hat_os_quant <- Lambda0_os_quant * risk_os
  R_hat_os_grid <- Lambda0_os_grid * risk_os
  R_hat_os_test_grid_10_90 <- Lambda0_os_test_grid_10_90 * risk_os
  R_hat_os_test_grid_20_90 <- Lambda0_os_test_grid_20_90 * risk_os
  R_hat_tl_quant <- Lambda0_tl_quant * risk_tl
  R_hat_tl_grid <- Lambda0_tl_grid * risk_tl
  R_hat_tl_test_grid_10_90 <- Lambda0_tl_test_grid_10_90 * risk_tl
  R_hat_tl_test_grid_20_90 <- Lambda0_tl_test_grid_20_90 * risk_tl
  S_hat_os_quant <- exp(-Lambda0_os_quant * risk_os)
  S_hat_os_grid <- exp(-Lambda0_os_grid * risk_os)
  S_hat_os_test_grid_10_90 <- exp(-Lambda0_os_test_grid_10_90 * risk_os)
  S_hat_os_test_grid_20_90 <- exp(-Lambda0_os_test_grid_20_90 * risk_os)
  S_hat_tl_quant <- exp(-Lambda0_tl_quant * risk_tl)
  S_hat_tl_grid <- exp(-Lambda0_tl_grid * risk_tl)
  S_hat_tl_test_grid_10_90 <- exp(-Lambda0_tl_test_grid_10_90 * risk_tl)
  S_hat_tl_test_grid_20_90 <- exp(-Lambda0_tl_test_grid_20_90 * risk_tl)
  IF_beta_x0 <- as.vector(IF_beta_hat_mat %*% x0)
  IF_R_hat_mat <- vapply(seq_along(Lambda0_tl_quant), function(j) {
    risk_tl * (IF_Lambda_hat_mat[, j] + Lambda0_tl_quant[j] * IF_beta_x0)
  }, numeric(n0))
  IF_R_hat_grid <- vapply(seq_along(Lambda0_tl_grid), function(j) {
    risk_tl * (IF_Lambda_hat_grid[, j] + Lambda0_tl_grid[j] * IF_beta_x0)
  }, numeric(n0))
  IF_R_os_mat <- vapply(seq_along(Lambda0_os_quant), function(j) {
    risk_os * (IF_Lambda_os_mat[, j] + Lambda0_os_quant[j] * IF_beta_x0)
  }, numeric(n0))
  IF_R_os_grid <- vapply(seq_along(Lambda0_os_grid), function(j) {
    risk_os * (IF_Lambda_os_grid[, j] + Lambda0_os_grid[j] * IF_beta_x0)
  }, numeric(n0))
  if (is.null(dim(IF_R_hat_mat))) IF_R_hat_mat <- matrix(IF_R_hat_mat, ncol = length(Lambda0_tl_quant))
  if (is.null(dim(IF_R_hat_grid))) IF_R_hat_grid <- matrix(IF_R_hat_grid, ncol = length(Lambda0_tl_grid))
  if (is.null(dim(IF_Lambda_os_mat))) IF_Lambda_os_mat <- matrix(IF_Lambda_os_mat, ncol = length(Lambda0_os_quant))
  if (is.null(dim(IF_Lambda_os_grid))) IF_Lambda_os_grid <- matrix(IF_Lambda_os_grid, ncol = length(Lambda0_os_grid))
  if (is.null(dim(IF_R_os_mat))) IF_R_os_mat <- matrix(IF_R_os_mat, ncol = length(Lambda0_os_quant))
  if (is.null(dim(IF_R_os_grid))) IF_R_os_grid <- matrix(IF_R_os_grid, ncol = length(Lambda0_os_grid))
  grad_S_quant <- lapply(seq_along(Lambda0_tl_quant), function(j) {
    c(
      -S_hat_tl_quant[j] * risk_tl * a_quant_mat_omega[j, ],
      -S_hat_tl_quant[j] * risk_tl * Lambda0_tl_quant[j] * x0
    )
  })
  grad_S_grid <- lapply(seq_along(Lambda0_tl_grid), function(j) {
    c(
      -S_hat_tl_grid[j] * risk_tl * a_grid_mat_omega[j, ],
      -S_hat_tl_grid[j] * risk_tl * Lambda0_tl_grid[j] * x0
    )
  })
  grad_S_quant_os <- lapply(seq_along(Lambda0_os_quant), function(j) {
    c(
      -S_hat_os_quant[j] * risk_os * a_quant_mat_omega_os[j, ],
      -S_hat_os_quant[j] * risk_os * Lambda0_os_quant[j] * x0
    )
  })
  grad_S_grid_os <- lapply(seq_along(Lambda0_os_grid), function(j) {
    c(
      -S_hat_os_grid[j] * risk_os * a_grid_mat_omega_os[j, ],
      -S_hat_os_grid[j] * risk_os * Lambda0_os_grid[j] * x0
    )
  })
  var_hat_S_omega <- vapply(grad_S_quant, function(g) {
    as.numeric(t(g) %*% Sigma_theta_hat_omega %*% g)
  }, numeric(1)) / n0
  var_hat_S_grid_omega <- vapply(grad_S_grid, function(g) {
    as.numeric(t(g) %*% Sigma_theta_hat_omega %*% g)
  }, numeric(1)) / n0
  var_hat_S_omega_zero <- vapply(grad_S_quant, function(g) {
    as.numeric(t(g) %*% Sigma_theta_hat_omega_zero %*% g)
  }, numeric(1)) / n0
  var_hat_S_grid_omega_zero <- vapply(grad_S_grid, function(g) {
    as.numeric(t(g) %*% Sigma_theta_hat_omega_zero %*% g)
  }, numeric(1)) / n0
  var_hat_S_omega_os <- vapply(grad_S_quant_os, function(g) {
    as.numeric(t(g) %*% Sigma_theta_hat_omega %*% g)
  }, numeric(1)) / n0
  var_hat_S_grid_omega_os <- vapply(grad_S_grid_os, function(g) {
    as.numeric(t(g) %*% Sigma_theta_hat_omega %*% g)
  }, numeric(1)) / n0
  var_hat_S_omega_zero_os <- vapply(grad_S_quant_os, function(g) {
    as.numeric(t(g) %*% Sigma_theta_hat_omega_zero %*% g)
  }, numeric(1)) / n0
  var_hat_S_grid_omega_zero_os <- vapply(grad_S_grid_os, function(g) {
    as.numeric(t(g) %*% Sigma_theta_hat_omega_zero %*% g)
  }, numeric(1)) / n0
  grad_R_quant <- lapply(seq_along(Lambda0_tl_quant), function(j) {
    risk_tl * c(a_quant_mat_omega[j, ], Lambda0_tl_quant[j] * x0)
  })
  grad_R_grid <- lapply(seq_along(Lambda0_tl_grid), function(j) {
    risk_tl * c(a_grid_mat_omega[j, ], Lambda0_tl_grid[j] * x0)
  })
  grad_R_quant_os <- lapply(seq_along(Lambda0_os_quant), function(j) {
    risk_os * c(a_quant_mat_omega_os[j, ], Lambda0_os_quant[j] * x0)
  })
  grad_R_grid_os <- lapply(seq_along(Lambda0_os_grid), function(j) {
    risk_os * c(a_grid_mat_omega_os[j, ], Lambda0_os_grid[j] * x0)
  })
  var_hat_R_omega <- vapply(grad_R_quant, function(g) {
    as.numeric(t(g) %*% Sigma_theta_hat_omega %*% g)
  }, numeric(1)) / n0
  var_hat_R_grid_omega <- vapply(grad_R_grid, function(g) {
    as.numeric(t(g) %*% Sigma_theta_hat_omega %*% g)
  }, numeric(1)) / n0
  var_hat_R_omega_zero <- vapply(grad_R_quant, function(g) {
    as.numeric(t(g) %*% Sigma_theta_hat_omega_zero %*% g)
  }, numeric(1)) / n0
  var_hat_R_grid_omega_zero <- vapply(grad_R_grid, function(g) {
    as.numeric(t(g) %*% Sigma_theta_hat_omega_zero %*% g)
  }, numeric(1)) / n0
  var_hat_R_omega_os <- vapply(grad_R_quant_os, function(g) {
    as.numeric(t(g) %*% Sigma_theta_hat_omega %*% g)
  }, numeric(1)) / n0
  var_hat_R_grid_omega_os <- vapply(grad_R_grid_os, function(g) {
    as.numeric(t(g) %*% Sigma_theta_hat_omega %*% g)
  }, numeric(1)) / n0
  var_hat_R_omega_zero_os <- vapply(grad_R_quant_os, function(g) {
    as.numeric(t(g) %*% Sigma_theta_hat_omega_zero %*% g)
  }, numeric(1)) / n0
  var_hat_R_grid_omega_zero_os <- vapply(grad_R_grid_os, function(g) {
    as.numeric(t(g) %*% Sigma_theta_hat_omega_zero %*% g)
  }, numeric(1)) / n0
  var_hat_S_omega <- pmax(var_hat_S_omega, 0)
  var_hat_S_grid_omega <- pmax(var_hat_S_grid_omega, 0)
  var_hat_S_omega_zero <- pmax(var_hat_S_omega_zero, 0)
  var_hat_S_grid_omega_zero <- pmax(var_hat_S_grid_omega_zero, 0)
  var_hat_S_omega_os <- pmax(var_hat_S_omega_os, 0)
  var_hat_S_grid_omega_os <- pmax(var_hat_S_grid_omega_os, 0)
  var_hat_S_omega_zero_os <- pmax(var_hat_S_omega_zero_os, 0)
  var_hat_S_grid_omega_zero_os <- pmax(var_hat_S_grid_omega_zero_os, 0)
  var_hat_R_omega <- pmax(var_hat_R_omega, 0)
  var_hat_R_grid_omega <- pmax(var_hat_R_grid_omega, 0)
  var_hat_R_omega_zero <- pmax(var_hat_R_omega_zero, 0)
  var_hat_R_grid_omega_zero <- pmax(var_hat_R_grid_omega_zero, 0)
  var_hat_R_omega_os <- pmax(var_hat_R_omega_os, 0)
  var_hat_R_grid_omega_os <- pmax(var_hat_R_grid_omega_os, 0)
  var_hat_R_omega_zero_os <- pmax(var_hat_R_omega_zero_os, 0)
  var_hat_R_grid_omega_zero_os <- pmax(var_hat_R_grid_omega_zero_os, 0)
  se_S_omega <- sqrt(var_hat_S_omega)
  se_S_omega_zero <- sqrt(var_hat_S_omega_zero)
  se_R_omega <- sqrt(var_hat_R_omega)
  se_R_omega_zero <- sqrt(var_hat_R_omega_zero)
  se_S_omega_os <- sqrt(var_hat_S_omega_os)
  se_S_omega_zero_os <- sqrt(var_hat_S_omega_zero_os)
  se_R_omega_os <- sqrt(var_hat_R_omega_os)
  se_R_omega_zero_os <- sqrt(var_hat_R_omega_zero_os)
  survival_labels <- paste0("S_q", sprintf("%02d", round(100 * quant_probs)))

  internal_checks <- list(
    x0_length = length(x0),
    dim_IF_theta_hat_mat = dim(IF_theta_hat_mat),
    dim_IF_gamma_hat_mat = dim(IF_gamma_hat_mat),
    dim_IF_Lambda_hat_mat = dim(IF_Lambda_hat_mat),
    theta_score_if_max_abs_diff = max(abs((theta_hat - Omega_hat_by_column %*% score_vec) - (theta_hat + colMeans(IF_theta_hat_mat)))),
    theta_update_max_abs_diff = max(abs(colMeans(IF_theta_hat_mat) - (theta_os - theta_hat))),
    gamma_update_max_abs_diff = max(abs(colMeans(IF_gamma_hat_mat) - (gamma_os - gamma_hat))),
    omega_vs_if_beta_var_max_abs_diff = max(abs(Sigma_beta_hat_omega / n0 - avar_beta / n0)),
    lambda_os_increasing = all(diff(Lambda0_os_grid) >= -1e-10),
    survival_in_unit_interval = all(S_hat_os_quant >= -1e-12 & S_hat_os_quant <= 1 + 1e-12) &&
      all(S_hat_os_grid >= -1e-12 & S_hat_os_grid <= 1 + 1e-12),
    survival_nonincreasing = all(diff(S_hat_os_grid) <= 1e-10)
  )

  if (isTRUE(debug_q10) && length(quant_times) > 0L) {
    q10_truth <- quant_times[1]^2
    cat("DEBUG q10 diagnostics:\n")
    cat("  quant_time[1] =", format(quant_times[1], digits = 16), "\n")
    cat("  true_Lambda_q10 =", format(q10_truth, digits = 16), "\n")
    cat("  Lambda0_os_quant[1] =", format(Lambda0_os_quant[1], digits = 16), "\n")
    cat("  pooled_curve_first_time = ", format(cumhaz_os_curve$time[1], digits = 16), "\n", sep = "")
    cat("  mean_IF_Lambda_q10 = ", format(mean(IF_Lambda_hat_mat[, 1]), digits = 16), "\n", sep = "")
    cat("  min/max_IF_Lambda_q10 = [",
        format(min(IF_Lambda_hat_mat[, 1]), digits = 16), ", ",
        format(max(IF_Lambda_hat_mat[, 1]), digits = 16), "]\n", sep = "")
  }

  list(
    n0 = n0,
    theta_hat = theta_hat,
    theta_db_by_column = theta_db_by_column,
    theta_os = theta_os,
    gamma_hat = gamma_hat,
    gamma_db_by_column = gamma_db_by_column,
    gamma_os = gamma_os,
    beta_hat = beta_hat,
    beta_db_by_column = beta_db_by_column,
    beta_os = beta_os,
    beta_hat_os = beta_hat_os,
    beta_hat_tl = beta_hat,
    se_beta = se_beta,
    se_beta_os = se_beta,
    Omega_hat_by_column = Omega_hat_by_column,
    Omega_hat_col = Omega_hat_by_column,
    delta_hat = delta_hat,
    U_target = score_vec,
    if_theta_db_by_column_mat = if_theta_db_by_column_mat,
    IF_theta_hat_mat = IF_theta_hat_mat,
    IF_beta_hat_mat = IF_beta_hat_mat,
    IF_gamma_hat_mat = IF_gamma_hat_mat,
    Sigma_theta_hat_omega = Sigma_theta_hat_omega,
    Sigma_gamma_hat_omega = Sigma_gamma_hat_omega,
    Sigma_theta_hat_omega_zero = Sigma_theta_hat_omega_zero,
    Sigma_gamma_hat_omega_zero = Sigma_gamma_hat_omega_zero,
    Sigma_hat_beta = avar_beta,
    Sigma_hat_beta_omega = Sigma_beta_hat_omega,
    Sigma_hat_beta_omega_zero = Sigma_beta_hat_omega_zero,
    x0 = x0,
    integration_grid = scb_grid,
    scb_test_grid = test_grid_10_90,
    scb_test_grid_10_90 = test_grid_10_90,
    scb_test_grid_20_90 = test_grid_20_90,
    a_grid_mat = a_grid_obj$a_mat,
    quant_times = quant_times,
    quant_labels = quant_labels,
    survival_labels = survival_labels,
    cumhaz_os_linearized_curve = cumhaz_os_linearized_curve,
    Lambda0_tl_quant = Lambda0_tl_quant,
    Lambda0_tl_grid = Lambda0_tl_grid,
    Lambda0_tl_test_grid = Lambda0_tl_test_grid_10_90,
    Lambda0_tl_test_grid_10_90 = Lambda0_tl_test_grid_10_90,
    Lambda0_tl_test_grid_20_90 = Lambda0_tl_test_grid_20_90,
    Lambda0_os_quant = Lambda0_os_quant,
    Lambda0_os_grid = Lambda0_os_grid,
    Lambda0_os_test_grid = Lambda0_os_test_grid_10_90,
    Lambda0_os_test_grid_10_90 = Lambda0_os_test_grid_10_90,
    Lambda0_os_test_grid_20_90 = Lambda0_os_test_grid_20_90,
    Lambda_hat_tl_quant = Lambda0_tl_quant,
    Lambda_hat_tl_grid = Lambda0_tl_grid,
    Lambda_hat_tl_test_grid = Lambda0_tl_test_grid_10_90,
    Lambda_hat_tl_test_grid_10_90 = Lambda0_tl_test_grid_10_90,
    Lambda_hat_tl_test_grid_20_90 = Lambda0_tl_test_grid_20_90,
    Lambda_hat_os_quant = Lambda0_os_quant,
    Lambda_hat_os_grid = Lambda0_os_grid,
    Lambda_hat_os_test_grid = Lambda0_os_test_grid_10_90,
    Lambda_hat_os_test_grid_10_90 = Lambda0_os_test_grid_10_90,
    Lambda_hat_os_test_grid_20_90 = Lambda0_os_test_grid_20_90,
    IF_Lambda_hat_mat = IF_Lambda_hat_mat,
    IF_Lambda_hat_grid = IF_Lambda_hat_grid,
    IF_Lambda_os_mat = IF_Lambda_os_mat,
    IF_Lambda_os_grid = IF_Lambda_os_grid,
    IF_Lambda_os_plugin_mat = IF_Lambda_os_mat,
    IF_Lambda_os_plugin_grid = IF_Lambda_os_grid,
    se_beta_omega = se_beta_omega,
    se_beta_omega_zero = se_beta_omega_zero,
    var_hat_Lambda_omega = var_hat_Lambda_omega,
    se_Lambda_omega = se_Lambda_omega,
    var_hat_Lambda_grid_omega = var_hat_Lambda_grid_omega,
    var_hat_Lambda_omega_zero = var_hat_Lambda_omega_zero,
    se_Lambda_omega_zero = se_Lambda_omega_zero,
    var_hat_Lambda_grid_omega_zero = var_hat_Lambda_grid_omega_zero,
    risk_os = risk_os,
    risk_tl = risk_tl,
    R_hat_os_quant = R_hat_os_quant,
    R_hat_os_grid = R_hat_os_grid,
    R_hat_os_test_grid = R_hat_os_test_grid_10_90,
    R_hat_os_test_grid_10_90 = R_hat_os_test_grid_10_90,
    R_hat_os_test_grid_20_90 = R_hat_os_test_grid_20_90,
    R_hat_tl_quant = R_hat_tl_quant,
    R_hat_tl_grid = R_hat_tl_grid,
    R_hat_tl_test_grid = R_hat_tl_test_grid_10_90,
    R_hat_tl_test_grid_10_90 = R_hat_tl_test_grid_10_90,
    R_hat_tl_test_grid_20_90 = R_hat_tl_test_grid_20_90,
    S_hat_os_quant = S_hat_os_quant,
    S_hat_os_grid = S_hat_os_grid,
    S_hat_os_test_grid = S_hat_os_test_grid_10_90,
    S_hat_os_test_grid_10_90 = S_hat_os_test_grid_10_90,
    S_hat_os_test_grid_20_90 = S_hat_os_test_grid_20_90,
    S_hat_tl_quant = S_hat_tl_quant,
    S_hat_tl_grid = S_hat_tl_grid,
    S_hat_tl_test_grid = S_hat_tl_test_grid_10_90,
    S_hat_tl_test_grid_10_90 = S_hat_tl_test_grid_10_90,
    S_hat_tl_test_grid_20_90 = S_hat_tl_test_grid_20_90,
    IF_beta_x0 = IF_beta_x0,
    IF_R_hat_mat = IF_R_hat_mat,
    IF_R_hat_grid = IF_R_hat_grid,
    IF_R_os_plugin_mat = IF_R_os_mat,
    IF_R_os_plugin_grid = IF_R_os_grid,
    var_hat_S_omega = var_hat_S_omega,
    se_S_omega = se_S_omega,
    var_hat_S_grid_omega = var_hat_S_grid_omega,
    var_hat_S_omega_zero = var_hat_S_omega_zero,
    se_S_omega_zero = se_S_omega_zero,
    var_hat_S_grid_omega_zero = var_hat_S_grid_omega_zero,
    var_hat_S_omega_os = var_hat_S_omega_os,
    se_S_omega_os = se_S_omega_os,
    var_hat_S_grid_omega_os = var_hat_S_grid_omega_os,
    var_hat_S_omega_zero_os = var_hat_S_omega_zero_os,
    se_S_omega_zero_os = se_S_omega_zero_os,
    var_hat_S_grid_omega_zero_os = var_hat_S_grid_omega_zero_os,
    var_hat_R_omega = var_hat_R_omega,
    se_R_omega = se_R_omega,
    var_hat_R_grid_omega = var_hat_R_grid_omega,
    var_hat_R_omega_zero = var_hat_R_omega_zero,
    se_R_omega_zero = se_R_omega_zero,
    var_hat_R_grid_omega_zero = var_hat_R_grid_omega_zero,
    var_hat_R_omega_os = var_hat_R_omega_os,
    se_R_omega_os = se_R_omega_os,
    var_hat_R_grid_omega_os = var_hat_R_grid_omega_os,
    var_hat_R_omega_zero_os = var_hat_R_omega_zero_os,
    se_R_omega_zero_os = se_R_omega_zero_os,
    var_hat_R_grid_omega_zero_os = var_hat_R_grid_omega_zero_os,
    var_hat_Lambda_omega_os = var_hat_Lambda_omega_os,
    se_Lambda_omega_os = se_Lambda_omega_os,
    var_hat_Lambda_grid_omega_os = var_hat_Lambda_grid_omega_os,
    var_hat_Lambda_omega_zero_os = var_hat_Lambda_omega_zero_os,
    se_Lambda_omega_zero_os = se_Lambda_omega_zero_os,
    var_hat_Lambda_grid_omega_zero_os = var_hat_Lambda_grid_omega_zero_os,
    checks = internal_checks
  )
}

compute_multiplier_process_intervals_from_if <- function(point_quant,
                                                         point_grid,
                                                         point_test_grid,
                                                         if_quant_mat,
                                                         if_grid_mat,
                                                         scb_grid,
                                                         scb_test_grid,
                                                         alpha = 0.05,
                                                         multiplier_boot = 500,
                                                         multiplier = c("normal", "rademacher"),
                                                         seed = NULL,
                                                         sd_floor = 1e-8,
                                                         xi = NULL) {
  multiplier <- match.arg(multiplier)
  n0 <- nrow(if_quant_mat)
  if (!is.finite(n0) || n0 < 1L) stop("Invalid target sample size in inference object.")
  if (!is.finite(multiplier_boot) || multiplier_boot < 50L) {
    stop("multiplier_boot should be at least 50.")
  }
  if (is.null(xi)) {
    if (!is.null(seed)) set.seed(seed)
    if (multiplier == "normal") {
      xi <- matrix(rnorm(multiplier_boot * n0), nrow = multiplier_boot, ncol = n0)
    } else {
      xi <- matrix(sample(c(-1, 1), size = multiplier_boot * n0, replace = TRUE), nrow = multiplier_boot, ncol = n0)
    }
  } else {
    if (!all(dim(xi) == c(multiplier_boot, n0))) {
      stop("Provided xi must have dimension multiplier_boot x n0.")
    }
  }

  z_quant <- (xi %*% if_quant_mat) / sqrt(n0)
  s_quant_boot <- sqrt(pmax(colMeans(z_quant^2), 0))
  se_quant_boot <- s_quant_boot / sqrt(n0)
  z_alpha <- stats::qnorm(1 - alpha / 2)
  ci_lower_quant <- point_quant - z_alpha * se_quant_boot
  ci_upper_quant <- point_quant + z_alpha * se_quant_boot

  z_grid <- (xi %*% if_grid_mat) / sqrt(n0)
  s_grid_boot <- sqrt(pmax(colMeans(z_grid^2), 0))
  s_grid_safe <- pmax(s_grid_boot, sd_floor)
  studentized <- abs(sweep(z_grid, 2, s_grid_safe, "/"))
  sup_stats <- apply(studentized, 1, max, na.rm = TRUE)
  scb_crit <- as.numeric(stats::quantile(sup_stats, probs = 1 - alpha, names = FALSE, na.rm = TRUE))
  scb_halfwidth <- scb_crit * s_grid_safe / sqrt(n0)
  scb_lower_grid <- point_grid - scb_halfwidth
  scb_upper_grid <- point_grid + scb_halfwidth
  scb_lower_test <- stats::approx(scb_grid, scb_lower_grid,
                                  xout = scb_test_grid, rule = 2, ties = "ordered")$y
  scb_upper_test <- stats::approx(scb_grid, scb_upper_grid,
                                  xout = scb_test_grid, rule = 2, ties = "ordered")$y
  pointwise_width <- ci_upper_quant - ci_lower_quant
  scb_width_grid <- scb_upper_grid - scb_lower_grid
  scb_width_mean <- mean(scb_upper_test - scb_lower_test, na.rm = TRUE)

  list(
    pointwise = list(
      lower = ci_lower_quant,
      upper = ci_upper_quant,
      se = se_quant_boot,
      width = pointwise_width
    ),
    scb = list(
      s_grid = s_grid_boot,
      crit = scb_crit,
      lower = scb_lower_grid,
      upper = scb_upper_grid,
      lower_test = scb_lower_test,
      upper_test = scb_upper_test,
      width_grid = scb_width_grid,
      width_mean = scb_width_mean
    ),
    checks = list(
      multiplier_quant_dim = dim(z_quant),
      multiplier_grid_dim = dim(z_grid),
      xi_dim = dim(xi)
    )
  )
}

compute_multiplier_cumh_intervals_from_if <- function(...) {
  compute_multiplier_process_intervals_from_if(...)
}

compute_multiplier_survival_intervals_from_risk_if <- function(point_quant_risk,
                                                               point_grid_risk,
                                                               point_test_grid_risk,
                                                               if_quant_mat_risk,
                                                               if_grid_mat_risk,
                                                               scb_grid,
                                                               scb_test_grid,
                                                               alpha = 0.05,
                                                               multiplier_boot = 500,
                                                               multiplier = c("normal", "rademacher"),
                                                               seed = NULL,
                                                               sd_floor = 1e-8,
                                                               xi = NULL) {
  risk_intervals <- compute_multiplier_process_intervals_from_if(
    point_quant = point_quant_risk,
    point_grid = point_grid_risk,
    point_test_grid = point_test_grid_risk,
    if_quant_mat = if_quant_mat_risk,
    if_grid_mat = if_grid_mat_risk,
    scb_grid = scb_grid,
    scb_test_grid = scb_test_grid,
    alpha = alpha,
    multiplier_boot = multiplier_boot,
    multiplier = match.arg(multiplier),
    seed = seed,
    sd_floor = sd_floor,
    xi = xi
  )
  z_alpha <- stats::qnorm(1 - alpha / 2)

  survival_lower_quant <- exp(-risk_intervals$pointwise$upper)
  survival_upper_quant <- exp(-risk_intervals$pointwise$lower)
  survival_lower_quant <- pmax(pmin(survival_lower_quant, 1), 0)
  survival_upper_quant <- pmax(pmin(survival_upper_quant, 1), 0)
  survival_width_quant <- survival_upper_quant - survival_lower_quant
  survival_se_equiv <- survival_width_quant / (2 * z_alpha)

  survival_lower_grid <- exp(-risk_intervals$scb$upper)
  survival_upper_grid <- exp(-risk_intervals$scb$lower)
  survival_lower_grid <- pmax(pmin(survival_lower_grid, 1), 0)
  survival_upper_grid <- pmax(pmin(survival_upper_grid, 1), 0)
  survival_lower_test <- exp(-risk_intervals$scb$upper_test)
  survival_upper_test <- exp(-risk_intervals$scb$lower_test)
  survival_lower_test <- pmax(pmin(survival_lower_test, 1), 0)
  survival_upper_test <- pmax(pmin(survival_upper_test, 1), 0)
  survival_width_grid <- survival_upper_grid - survival_lower_grid
  survival_width_mean <- mean(survival_upper_test - survival_lower_test, na.rm = TRUE)

  list(
    pointwise = list(
      lower = survival_lower_quant,
      upper = survival_upper_quant,
      se = survival_se_equiv,
      se_risk = risk_intervals$pointwise$se,
      width = survival_width_quant,
      lower_risk = risk_intervals$pointwise$lower,
      upper_risk = risk_intervals$pointwise$upper
    ),
    scb = list(
      s_grid_risk = risk_intervals$scb$s_grid,
      crit_R = risk_intervals$scb$crit,
      lower = survival_lower_grid,
      upper = survival_upper_grid,
      lower_test = survival_lower_test,
      upper_test = survival_upper_test,
      width_grid = survival_width_grid,
      width_mean = survival_width_mean,
      lower_risk = risk_intervals$scb$lower,
      upper_risk = risk_intervals$scb$upper,
      lower_test_risk = risk_intervals$scb$lower_test,
      upper_test_risk = risk_intervals$scb$upper_test
    ),
    checks = risk_intervals$checks
  )
}

compute_multiplier_process_intervals <- function(inference,
                                                 alpha = 0.05,
                                                 multiplier_boot = 500,
                                                 multiplier = c("normal", "rademacher"),
                                                 seed = NULL,
                                                 sd_floor = 1e-8) {
  multiplier <- match.arg(multiplier)
  n0 <- nrow(inference$IF_Lambda_os_mat)
  if (!is.null(seed)) set.seed(seed)
  xi <- if (multiplier == "normal") {
    matrix(rnorm(multiplier_boot * n0), nrow = multiplier_boot, ncol = n0)
  } else {
    matrix(sample(c(-1, 1), size = multiplier_boot * n0, replace = TRUE), nrow = multiplier_boot, ncol = n0)
  }
  list(
    xi = xi,
    cumh = list(
      os = compute_multiplier_process_intervals_from_if(
        point_quant = inference$Lambda_hat_os_quant,
        point_grid = inference$Lambda_hat_os_grid,
        point_test_grid = inference$Lambda_hat_os_test_grid,
        if_quant_mat = inference$IF_Lambda_os_mat,
        if_grid_mat = inference$IF_Lambda_os_grid,
        scb_grid = inference$integration_grid,
        scb_test_grid = inference$scb_test_grid,
        alpha = alpha,
        multiplier_boot = multiplier_boot,
        multiplier = multiplier,
        sd_floor = sd_floor,
        xi = xi
      ),
      os_plugin = compute_multiplier_process_intervals_from_if(
        point_quant = inference$Lambda_hat_os_quant,
        point_grid = inference$Lambda_hat_os_grid,
        point_test_grid = inference$Lambda_hat_os_test_grid,
        if_quant_mat = inference$IF_Lambda_os_plugin_mat,
        if_grid_mat = inference$IF_Lambda_os_plugin_grid,
        scb_grid = inference$integration_grid,
        scb_test_grid = inference$scb_test_grid,
        alpha = alpha,
        multiplier_boot = multiplier_boot,
        multiplier = multiplier,
        sd_floor = sd_floor,
        xi = xi
      )
    ),
    survival = list(
      os = compute_multiplier_survival_intervals_from_risk_if(
        point_quant_risk = inference$R_hat_os_quant,
        point_grid_risk = inference$R_hat_os_grid,
        point_test_grid_risk = inference$R_hat_os_test_grid,
        if_quant_mat_risk = inference$IF_R_hat_mat,
        if_grid_mat_risk = inference$IF_R_hat_grid,
        scb_grid = inference$integration_grid,
        scb_test_grid = inference$scb_test_grid,
        alpha = alpha,
        multiplier_boot = multiplier_boot,
        multiplier = multiplier,
        sd_floor = sd_floor,
        xi = xi
      ),
      os_plugin = compute_multiplier_survival_intervals_from_risk_if(
        point_quant_risk = inference$R_hat_os_quant,
        point_grid_risk = inference$R_hat_os_grid,
        point_test_grid_risk = inference$R_hat_os_test_grid,
        if_quant_mat_risk = inference$IF_R_os_plugin_mat,
        if_grid_mat_risk = inference$IF_R_os_plugin_grid,
        scb_grid = inference$integration_grid,
        scb_test_grid = inference$scb_test_grid,
        alpha = alpha,
        multiplier_boot = multiplier_boot,
        multiplier = multiplier,
        sd_floor = sd_floor,
        xi = xi
      )
    )
  )
}

compute_multiplier_cumh_intervals <- function(...) {
  compute_multiplier_process_intervals(...)
}

build_replication_records <- function(inference, beta_true, design_row, rep_id, c_sieve,
                                      z_alpha = 1.96, cumh_bootstrap = NULL,
                                      survival_bootstrap = NULL) {
  beta_truth <- beta_true[seq_along(inference$beta_hat)]
  beta_records_base <- data.frame(
    replication = rep_id,
    parameter_type = "beta",
    parameter_name = paste0("beta_", seq_along(inference$beta_hat)),
    n0 = design_row$n0,
    nA = design_row$nA,
    beta2 = design_row$beta2,
    kappa = design_row$kappa,
    c = design_row$c,
    shift = design_row$shift,
    c_sieve = c_sieve,
    estimate = inference$beta_hat_os,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    covered = NA_real_,
    se = NA_real_,
    bias = inference$beta_hat_os - beta_truth,
    estimate_sievetl = inference$beta_hat_tl,
    bias_sievetl = inference$beta_hat_tl - beta_truth,
    ci_width = NA_real_,
    variance_estimator = NA_character_
  )
  beta_lower_omega <- inference$beta_hat_os - z_alpha * inference$se_beta_omega
  beta_upper_omega <- inference$beta_hat_os + z_alpha * inference$se_beta_omega
  beta_records_omega <- beta_records_base
  beta_records_omega$ci_lower <- beta_lower_omega
  beta_records_omega$ci_upper <- beta_upper_omega
  beta_records_omega$covered <- as.numeric(beta_truth >= beta_lower_omega & beta_truth <= beta_upper_omega)
  beta_records_omega$se <- inference$se_beta_omega
  beta_records_omega$ci_width <- beta_upper_omega - beta_lower_omega
  beta_records_omega$variance_estimator <- "omega"
  beta_lower_omega_zero <- inference$beta_hat_os - z_alpha * inference$se_beta_omega_zero
  beta_upper_omega_zero <- inference$beta_hat_os + z_alpha * inference$se_beta_omega_zero
  beta_records_omega_zero <- beta_records_base
  beta_records_omega_zero$ci_lower <- beta_lower_omega_zero
  beta_records_omega_zero$ci_upper <- beta_upper_omega_zero
  beta_records_omega_zero$covered <- as.numeric(beta_truth >= beta_lower_omega_zero & beta_truth <= beta_upper_omega_zero)
  beta_records_omega_zero$se <- inference$se_beta_omega_zero
  beta_records_omega_zero$ci_width <- beta_upper_omega_zero - beta_lower_omega_zero
  beta_records_omega_zero$variance_estimator <- "omega_zero"

  lambda_truth <- inference$quant_times^2
  if (is.null(cumh_bootstrap)) {
    stop("cumh_bootstrap must be provided for cumulative-hazard coverage records.")
  }
  lambda_lower_plugin <- cumh_bootstrap$os$pointwise$lower
  lambda_upper_plugin <- cumh_bootstrap$os$pointwise$upper
  lambda_width_os <- cumh_bootstrap$os$pointwise$width
  lambda_records_plugin_base <- data.frame(
    replication = rep_id,
    parameter_type = "cumh_plugin",
    parameter_name = inference$quant_labels,
    n0 = design_row$n0,
    nA = design_row$nA,
    beta2 = design_row$beta2,
    kappa = design_row$kappa,
    c = design_row$c,
    shift = design_row$shift,
    c_sieve = c_sieve,
    estimate = inference$Lambda_hat_os_quant,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    covered = NA_real_,
    se = NA_real_,
    bias = inference$Lambda_hat_os_quant - lambda_truth,
    estimate_sievetl = inference$Lambda_hat_tl_quant,
    bias_sievetl = inference$Lambda_hat_tl_quant - lambda_truth,
    ci_width = NA_real_,
    variance_estimator = NA_character_
  )
  lambda_records_plugin_omega <- lambda_records_plugin_base
  lambda_records_plugin_omega$ci_lower <- inference$Lambda_hat_os_quant - z_alpha * inference$se_Lambda_omega
  lambda_records_plugin_omega$ci_upper <- inference$Lambda_hat_os_quant + z_alpha * inference$se_Lambda_omega
  lambda_lower_omega <- lambda_records_plugin_omega$ci_lower
  lambda_upper_omega <- lambda_records_plugin_omega$ci_upper
  lambda_records_plugin_omega$ci_lower <- lambda_lower_omega
  lambda_records_plugin_omega$ci_upper <- lambda_upper_omega
  lambda_records_plugin_omega$covered <- as.numeric(lambda_truth >= lambda_lower_omega & lambda_truth <= lambda_upper_omega)
  lambda_records_plugin_omega$se <- inference$se_Lambda_omega
  lambda_records_plugin_omega$ci_width <- lambda_upper_omega - lambda_lower_omega
  lambda_records_plugin_omega$variance_estimator <- "omega"
  lambda_records_plugin_omega_zero <- lambda_records_plugin_base
  lambda_lower_omega_zero <- inference$Lambda_hat_os_quant - z_alpha * inference$se_Lambda_omega_zero
  lambda_upper_omega_zero <- inference$Lambda_hat_os_quant + z_alpha * inference$se_Lambda_omega_zero
  lambda_records_plugin_omega_zero$ci_lower <- lambda_lower_omega_zero
  lambda_records_plugin_omega_zero$ci_upper <- lambda_upper_omega_zero
  lambda_records_plugin_omega_zero$covered <- as.numeric(lambda_truth >= lambda_lower_omega_zero & lambda_truth <= lambda_upper_omega_zero)
  lambda_records_plugin_omega_zero$se <- inference$se_Lambda_omega_zero
  lambda_records_plugin_omega_zero$ci_width <- lambda_upper_omega_zero - lambda_lower_omega_zero
  lambda_records_plugin_omega_zero$variance_estimator <- "omega_zero"
  lambda_lower_omega_os <- inference$Lambda_hat_os_quant - z_alpha * inference$se_Lambda_omega_os
  lambda_upper_omega_os <- inference$Lambda_hat_os_quant + z_alpha * inference$se_Lambda_omega_os
  lambda_records_plugin_omega_os <- lambda_records_plugin_base
  lambda_records_plugin_omega_os$ci_lower <- lambda_lower_omega_os
  lambda_records_plugin_omega_os$ci_upper <- lambda_upper_omega_os
  lambda_records_plugin_omega_os$covered <- as.numeric(lambda_truth >= lambda_lower_omega_os & lambda_truth <= lambda_upper_omega_os)
  lambda_records_plugin_omega_os$se <- inference$se_Lambda_omega_os
  lambda_records_plugin_omega_os$ci_width <- lambda_upper_omega_os - lambda_lower_omega_os
  lambda_records_plugin_omega_os$variance_estimator <- "omega_os"
  lambda_lower_omega_zero_os <- inference$Lambda_hat_os_quant - z_alpha * inference$se_Lambda_omega_zero_os
  lambda_upper_omega_zero_os <- inference$Lambda_hat_os_quant + z_alpha * inference$se_Lambda_omega_zero_os
  lambda_records_plugin_omega_zero_os <- lambda_records_plugin_base
  lambda_records_plugin_omega_zero_os$ci_lower <- lambda_lower_omega_zero_os
  lambda_records_plugin_omega_zero_os$ci_upper <- lambda_upper_omega_zero_os
  lambda_records_plugin_omega_zero_os$covered <- as.numeric(lambda_truth >= lambda_lower_omega_zero_os & lambda_truth <= lambda_upper_omega_zero_os)
  lambda_records_plugin_omega_zero_os$se <- inference$se_Lambda_omega_zero_os
  lambda_records_plugin_omega_zero_os$ci_width <- lambda_upper_omega_zero_os - lambda_lower_omega_zero_os
  lambda_records_plugin_omega_zero_os$variance_estimator <- "omega_zero_os"

  lambda_lower_boot <- cumh_bootstrap$os$pointwise$lower
  lambda_upper_boot <- cumh_bootstrap$os$pointwise$upper
  lambda_records_boot <- data.frame(
    replication = rep_id,
    parameter_type = "cumh_boot",
    parameter_name = inference$quant_labels,
    n0 = design_row$n0,
    nA = design_row$nA,
    beta2 = design_row$beta2,
    kappa = design_row$kappa,
    c = design_row$c,
    shift = design_row$shift,
    c_sieve = c_sieve,
    estimate = NA_real_,
    ci_lower = lambda_lower_boot,
    ci_upper = lambda_upper_boot,
    covered = as.numeric(lambda_truth >= lambda_lower_boot & lambda_truth <= lambda_upper_boot),
    se = NA_real_,
    bias = NA_real_,
    estimate_sievetl = NA_real_,
    ci_width = NA_real_,
    bias_sievetl = NA_real_,
    variance_estimator = "if"
  )

  make_scb_record <- function(parameter_type, parameter_name, truth_grid, test_grid, scb_obj,
                              variance_estimator = "if") {
    lower_test <- stats::approx(inference$integration_grid, scb_obj$lower,
                                xout = test_grid, rule = 2, ties = "ordered")$y
    upper_test <- stats::approx(inference$integration_grid, scb_obj$upper,
                                xout = test_grid, rule = 2, ties = "ordered")$y
    data.frame(
      replication = rep_id,
      parameter_type = parameter_type,
      parameter_name = parameter_name,
      n0 = design_row$n0,
      nA = design_row$nA,
      beta2 = design_row$beta2,
      kappa = design_row$kappa,
      c = design_row$c,
      shift = design_row$shift,
      c_sieve = c_sieve,
      estimate = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      covered = as.numeric(all(truth_grid >= lower_test & truth_grid <= upper_test)),
      se = NA_real_,
      bias = NA_real_,
      estimate_sievetl = NA_real_,
      ci_width = mean(upper_test - lower_test, na.rm = TRUE),
      bias_sievetl = NA_real_,
      variance_estimator = variance_estimator
    )
  }
  scb_record_10_90 <- make_scb_record(
    parameter_type = "cumh_scb",
    parameter_name = "scb_t10_t90",
    truth_grid = inference$scb_test_grid_10_90^2,
    test_grid = inference$scb_test_grid_10_90,
    scb_obj = cumh_bootstrap$os$scb
  )
  scb_record_20_90 <- make_scb_record(
    parameter_type = "cumh_scb",
    parameter_name = "scb_t20_t90",
    truth_grid = inference$scb_test_grid_20_90^2,
    test_grid = inference$scb_test_grid_20_90,
    scb_obj = cumh_bootstrap$os$scb
  )
  scb_record_10_90_os_plugin <- make_scb_record(
    parameter_type = "cumh_scb",
    parameter_name = "scb_t10_t90",
    truth_grid = inference$scb_test_grid_10_90^2,
    test_grid = inference$scb_test_grid_10_90,
    scb_obj = cumh_bootstrap$os_plugin$scb,
    variance_estimator = "if_os"
  )
  scb_record_20_90_os_plugin <- make_scb_record(
    parameter_type = "cumh_scb",
    parameter_name = "scb_t20_t90",
    truth_grid = inference$scb_test_grid_20_90^2,
    test_grid = inference$scb_test_grid_20_90,
    scb_obj = cumh_bootstrap$os_plugin$scb,
    variance_estimator = "if_os"
  )

  if (is.null(survival_bootstrap)) {
    stop("survival_bootstrap must be provided for survival coverage records.")
  }
  risk_true <- exp(sum(beta_true * inference$x0))
  S_truth_quant <- exp(-(inference$quant_times^2) * risk_true)
  survival_pointwise_records_base <- data.frame(
    replication = rep_id,
    parameter_type = "survival_pointwise",
    parameter_name = inference$survival_labels,
    n0 = design_row$n0,
    nA = design_row$nA,
    beta2 = design_row$beta2,
    kappa = design_row$kappa,
    c = design_row$c,
    shift = design_row$shift,
    c_sieve = c_sieve,
    estimate = inference$S_hat_os_quant,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    covered = NA_real_,
    se = NA_real_,
    bias = inference$S_hat_os_quant - S_truth_quant,
    estimate_sievetl = inference$S_hat_tl_quant,
    bias_sievetl = inference$S_hat_tl_quant - S_truth_quant,
    ci_width = NA_real_,
    variance_estimator = NA_character_
  )
  R_lower_omega <- inference$R_hat_os_quant - z_alpha * inference$se_R_omega
  R_upper_omega <- inference$R_hat_os_quant + z_alpha * inference$se_R_omega
  survival_lower_omega <- exp(-R_upper_omega)
  survival_upper_omega <- exp(-R_lower_omega)
  survival_lower_omega <- pmax(pmin(survival_lower_omega, 1), 0)
  survival_upper_omega <- pmax(pmin(survival_upper_omega, 1), 0)
  survival_se_omega_equiv <- (survival_upper_omega - survival_lower_omega) / (2 * z_alpha)
  survival_pointwise_records_omega <- survival_pointwise_records_base
  survival_pointwise_records_omega$ci_lower <- survival_lower_omega
  survival_pointwise_records_omega$ci_upper <- survival_upper_omega
  survival_pointwise_records_omega$covered <- as.numeric(S_truth_quant >= survival_lower_omega & S_truth_quant <= survival_upper_omega)
  survival_pointwise_records_omega$se <- survival_se_omega_equiv
  survival_pointwise_records_omega$ci_width <- survival_upper_omega - survival_lower_omega
  survival_pointwise_records_omega$variance_estimator <- "omega"
  R_lower_omega_zero <- inference$R_hat_os_quant - z_alpha * inference$se_R_omega_zero
  R_upper_omega_zero <- inference$R_hat_os_quant + z_alpha * inference$se_R_omega_zero
  survival_lower_omega_zero <- exp(-R_upper_omega_zero)
  survival_upper_omega_zero <- exp(-R_lower_omega_zero)
  survival_lower_omega_zero <- pmax(pmin(survival_lower_omega_zero, 1), 0)
  survival_upper_omega_zero <- pmax(pmin(survival_upper_omega_zero, 1), 0)
  survival_se_omega_zero_equiv <- (survival_upper_omega_zero - survival_lower_omega_zero) / (2 * z_alpha)
  survival_pointwise_records_omega_zero <- survival_pointwise_records_base
  survival_pointwise_records_omega_zero$ci_lower <- survival_lower_omega_zero
  survival_pointwise_records_omega_zero$ci_upper <- survival_upper_omega_zero
  survival_pointwise_records_omega_zero$covered <- as.numeric(S_truth_quant >= survival_lower_omega_zero & S_truth_quant <= survival_upper_omega_zero)
  survival_pointwise_records_omega_zero$se <- survival_se_omega_zero_equiv
  survival_pointwise_records_omega_zero$ci_width <- survival_upper_omega_zero - survival_lower_omega_zero
  survival_pointwise_records_omega_zero$variance_estimator <- "omega_zero"
  R_lower_omega_os <- inference$R_hat_os_quant - z_alpha * inference$se_R_omega_os
  R_upper_omega_os <- inference$R_hat_os_quant + z_alpha * inference$se_R_omega_os
  survival_lower_omega_os <- exp(-R_upper_omega_os)
  survival_upper_omega_os <- exp(-R_lower_omega_os)
  survival_lower_omega_os <- pmax(pmin(survival_lower_omega_os, 1), 0)
  survival_upper_omega_os <- pmax(pmin(survival_upper_omega_os, 1), 0)
  survival_se_omega_os_equiv <- (survival_upper_omega_os - survival_lower_omega_os) / (2 * z_alpha)
  survival_pointwise_records_omega_os <- survival_pointwise_records_base
  survival_pointwise_records_omega_os$ci_lower <- survival_lower_omega_os
  survival_pointwise_records_omega_os$ci_upper <- survival_upper_omega_os
  survival_pointwise_records_omega_os$covered <- as.numeric(S_truth_quant >= survival_lower_omega_os & S_truth_quant <= survival_upper_omega_os)
  survival_pointwise_records_omega_os$se <- survival_se_omega_os_equiv
  survival_pointwise_records_omega_os$ci_width <- survival_upper_omega_os - survival_lower_omega_os
  survival_pointwise_records_omega_os$variance_estimator <- "omega_os"
  R_lower_omega_zero_os <- inference$R_hat_os_quant - z_alpha * inference$se_R_omega_zero_os
  R_upper_omega_zero_os <- inference$R_hat_os_quant + z_alpha * inference$se_R_omega_zero_os
  survival_lower_omega_zero_os <- exp(-R_upper_omega_zero_os)
  survival_upper_omega_zero_os <- exp(-R_lower_omega_zero_os)
  survival_lower_omega_zero_os <- pmax(pmin(survival_lower_omega_zero_os, 1), 0)
  survival_upper_omega_zero_os <- pmax(pmin(survival_upper_omega_zero_os, 1), 0)
  survival_se_omega_zero_os_equiv <- (survival_upper_omega_zero_os - survival_lower_omega_zero_os) / (2 * z_alpha)
  survival_pointwise_records_omega_zero_os <- survival_pointwise_records_base
  survival_pointwise_records_omega_zero_os$ci_lower <- survival_lower_omega_zero_os
  survival_pointwise_records_omega_zero_os$ci_upper <- survival_upper_omega_zero_os
  survival_pointwise_records_omega_zero_os$covered <- as.numeric(S_truth_quant >= survival_lower_omega_zero_os & S_truth_quant <= survival_upper_omega_zero_os)
  survival_pointwise_records_omega_zero_os$se <- survival_se_omega_zero_os_equiv
  survival_pointwise_records_omega_zero_os$ci_width <- survival_upper_omega_zero_os - survival_lower_omega_zero_os
  survival_pointwise_records_omega_zero_os$variance_estimator <- "omega_zero_os"
  survival_scb_record_10_90 <- make_scb_record(
    parameter_type = "survival_scb",
    parameter_name = "scb_t10_t90",
    truth_grid = exp(-(inference$scb_test_grid_10_90^2) * risk_true),
    test_grid = inference$scb_test_grid_10_90,
    scb_obj = survival_bootstrap$os$scb
  )
  survival_scb_record_20_90 <- make_scb_record(
    parameter_type = "survival_scb",
    parameter_name = "scb_t20_t90",
    truth_grid = exp(-(inference$scb_test_grid_20_90^2) * risk_true),
    test_grid = inference$scb_test_grid_20_90,
    scb_obj = survival_bootstrap$os$scb
  )
  survival_scb_record_10_90_os_plugin <- make_scb_record(
    parameter_type = "survival_scb",
    parameter_name = "scb_t10_t90",
    truth_grid = exp(-(inference$scb_test_grid_10_90^2) * risk_true),
    test_grid = inference$scb_test_grid_10_90,
    scb_obj = survival_bootstrap$os_plugin$scb,
    variance_estimator = "if_os"
  )
  survival_scb_record_20_90_os_plugin <- make_scb_record(
    parameter_type = "survival_scb",
    parameter_name = "scb_t20_t90",
    truth_grid = exp(-(inference$scb_test_grid_20_90^2) * risk_true),
    test_grid = inference$scb_test_grid_20_90,
    scb_obj = survival_bootstrap$os_plugin$scb,
    variance_estimator = "if_os"
  )

  dplyr::bind_rows(
    beta_records_omega,
    beta_records_omega_zero,
    lambda_records_plugin_omega,
    lambda_records_plugin_omega_zero,
    lambda_records_plugin_omega_os,
    lambda_records_plugin_omega_zero_os,
    lambda_records_boot,
    scb_record_10_90,
    scb_record_20_90,
    scb_record_10_90_os_plugin,
    scb_record_20_90_os_plugin,
    survival_pointwise_records_omega,
    survival_pointwise_records_omega_zero,
    survival_pointwise_records_omega_os,
    survival_pointwise_records_omega_zero_os,
    survival_scb_record_10_90,
    survival_scb_record_20_90,
    survival_scb_record_10_90_os_plugin,
    survival_scb_record_20_90_os_plugin
  )
}

default_sievetl_coverage_design <- function(beta2 = 0.5,
                                            kappa = 2,
                                            c = 1.0,
                                            shift = FALSE,
                                            n0 = 100,
                                            nA = 600) {
  data.frame(
    beta2 = beta2,
    kappa = kappa,
    shift = shift,
    n0 = n0,
    nA = nA,
    c = c
  )
}

summarize_coverage_records <- function(records) {
  if (!("variance_estimator" %in% names(records))) {
    records$variance_estimator <- "if"
  }
  group_keys <- c("variance_estimator", "parameter_type", "parameter_name", "n0", "nA", "beta2", "kappa", "c", "shift", "c_sieve")
  split_keys <- interaction(records[group_keys], drop = TRUE, lex.order = TRUE)
  grouped_records <- split(records, split_keys)

  summary_rows <- lapply(grouped_records, function(df) {
    has_estimate <- any(is.finite(df$estimate))
    has_se <- any(is.finite(df$se))
    has_bias <- any(is.finite(df$bias))
    has_estimate_sievetl <- "estimate_sievetl" %in% names(df) && any(is.finite(df$estimate_sievetl))
    has_bias_sievetl <- "bias_sievetl" %in% names(df) && any(is.finite(df$bias_sievetl))
    has_width <- "ci_width" %in% names(df) && any(is.finite(df$ci_width))

    data.frame(
      variance_estimator = df$variance_estimator[1],
      parameter_type = df$parameter_type[1],
      parameter_name = df$parameter_name[1],
      n0 = df$n0[1],
      nA = df$nA[1],
      beta2 = df$beta2[1],
      kappa = df$kappa[1],
      c = df$c[1],
      shift = df$shift[1],
      c_sieve = df$c_sieve[1],
      bias = if (has_bias) mean(df$bias, na.rm = TRUE) else NA_real_,
      empirical_se = if (has_estimate) sd(df$estimate, na.rm = TRUE) else NA_real_,
      bias_sievetl = if (has_bias_sievetl) mean(df$bias_sievetl, na.rm = TRUE) else NA_real_,
      empirical_se_sievetl = if (has_estimate_sievetl) sd(df$estimate_sievetl, na.rm = TRUE) else NA_real_,
      avg_estimated_se = if (has_se) stats::median(df$se, na.rm = TRUE) else NA_real_,
      coverage = mean(df$covered, na.rm = TRUE),
      avg_ci_width = if (has_width) mean(df$ci_width, na.rm = TRUE) else NA_real_,
      n_used = length(unique(df$replication))
    )
  })

  do.call(rbind, summary_rows)
}

# ------------------------ public main entry point ------------------------

run_sievetl_plugin_coverage_study <- function(nsim = 100,
                                              design = default_sievetl_coverage_design(),
                                              L = 5,
                                              c_sieve = 1,
                                              lambda_base_dir = getwd(),
                                              x0 = NULL,
                                              quant_probs = seq(0.1, 0.9, by = 0.1),
                                              integration_grid_size = 500,
                                              conf_level = 0.95,
                                              multiplier_boot = 500,
                                              multiplier = "normal",
                                              seed_base = 100000,
                                              output_csv = NULL,
                                              verbose = TRUE,
                                              debug_q10 = FALSE) {
  z_alpha <- qnorm(1 - (1 - conf_level) / 2)
  if (is.null(x0)) {
    x0 <- rep(0.5, L)
  } else {
    x0 <- as.numeric(x0)
  }
  if (length(x0) != L) {
    stop("x0 must have length L.")
  }
  record_list <- list()
  replication_list <- list()
  rec_idx <- 1L

  for (d in seq_len(nrow(design))) {
    row <- design[d, ]
    tuning_row <- tryCatch(
      resolve_sievetl_tuning_row(
        row$beta2, row$kappa, row$c, row$n0, row$nA,
        L = L, shift = isTRUE(row$shift), base_dir = lambda_base_dir
      ),
      error = function(e) NULL
    )
    if (is.null(tuning_row)) next
    penalties <- list(
      c_zeta = as.numeric(tuning_row$c_zeta_hat[1]),
      c_eta = as.numeric(tuning_row$c_eta_hat[1]),
      lambda_zeta = as.numeric(tuning_row$lambda_zeta_hat[1]),
      lambda_eta = as.numeric(tuning_row$lambda_eta_hat[1])
    )
    c_sieve_use <- as.numeric(tuning_row$c_hat[1])

    if (verbose) {
      message(sprintf(
        "Running beta2=%.2f, kappa=%s, shift=%s, n0=%d, nA=%d, c=%.2f, C_hat=%.2f, nsim=%d",
        row$beta2, format_kappa_tag(row$kappa), ifelse(isTRUE(row$shift), "TRUE", "FALSE"),
        row$n0, row$nA, row$c, c_sieve_use, nsim
      ))
    }

    for (rep_id in seq_len(nsim)) {
      set.seed(seed_base + 100000L * d + rep_id)

      sim_data <- tryCatch(
        simulate_once_all(
          n_target = row$n0,
          n_source = row$nA,
          beta2_source = row$beta2,
          kappa_source = row$kappa,
          c1 = row$c,
          covariate_shift = isTRUE(row$shift),
          L = L
        ),
        error = function(e) {
          cat(sprintf("Replication %d: simulate_once_all failed: %s\n", rep_id, conditionMessage(e)))
          NULL
        }
      )
      if (is.null(sim_data)) next
      sim_data$beta_true <- true_beta_target(L)

      fit_obj <- tryCatch(
        {
          fit_exact <- sievetl_approx(
            target = sim_data$target,
            source = sim_data$source,
            L = L,
            c_zeta = penalties$c_zeta,
            c_eta = penalties$c_eta,
            n0_penalty = row$n0,
            c = c_sieve_use
          )
          spline_obj <- build_g_funcs_exact_from_run(
            target = sim_data$target,
            source = sim_data$source,
            cspline = c_sieve_use
          )
          list(
            fit = fit_exact,
            g_funcs = spline_obj$g_funcs,
            p_hat = spline_obj$p_hat
          )
        },
        error = function(e) {
          cat(sprintf("Replication %d: sievetl_approx/build_g_funcs_exact_from_run failed: %s\n", rep_id, conditionMessage(e)))
          NULL
        }
      )
      if (is.null(fit_obj)) next

      d_n <- fit_obj$p_hat + L
      by_column_lambda_file <- default_by_column_lambda_path(row[1, ], c_sieve = c_sieve_use, base_dir = lambda_base_dir, L = L)
      by_column_lambda <- tryCatch(
        read_by_column_lambda_file(by_column_lambda_file, d_n),
        error = function(e) {
          cat(sprintf("Replication %d: read_by_column_lambda_file failed: %s\n", rep_id, conditionMessage(e)))
          NULL
        }
      )
      if (is.null(by_column_lambda)) next
      delta_hat <- c(fit_obj$fit$zeta_hat, fit_obj$fit$eta_hat)

      quant_times <- reference_target_quant_times(L = L, probs = quant_probs)

      inference <- tryCatch(
        compute_plugin_inference(
          target = sim_data$target,
          source = sim_data$source,
          beta_hat = fit_obj$fit$beta.hat,
          gamma_hat = fit_obj$fit$gamma.hat,
          beta_hat_se_base = fit_obj$fit$beta_hat_optim,
          gamma_hat_se_base = fit_obj$fit$gamma_hat_optim,
          g_funcs = fit_obj$g_funcs,
          lambda_delta_vec = by_column_lambda$lambda,
          delta_hat = delta_hat,
          x0 = x0,
          quant_times = quant_times,
          quant_probs = quant_probs,
          integration_grid_size = integration_grid_size,
          debug_q10 = debug_q10
        ),
        error = function(e) {
          cat(sprintf("Replication %d: compute_plugin_inference failed: %s\n", rep_id, conditionMessage(e)))
          NULL
        }
      )
      if (is.null(inference)) next

      process_bootstrap <- tryCatch(
        compute_multiplier_process_intervals(
          inference = inference,
          alpha = 1 - conf_level,
          multiplier_boot = multiplier_boot,
          multiplier = multiplier,
          seed = seed_base + 500000L * d + rep_id
        ),
        error = function(e) {
          cat(sprintf("Replication %d: compute_multiplier_process_intervals failed: %s\n", rep_id, conditionMessage(e)))
          NULL
        }
      )
      if (is.null(process_bootstrap)) next
      cumh_bootstrap <- process_bootstrap$cumh
      survival_bootstrap <- process_bootstrap$survival

      if (isTRUE(debug_q10) && length(inference$quant_times) > 0L) {
        cat("  multiplier_CI_q10 = [",
            format(cumh_bootstrap$os$pointwise$lower[1], digits = 16), ", ",
            format(cumh_bootstrap$os$pointwise$upper[1], digits = 16), "]\n", sep = "")
      }

      beta_ci_lower_os <- inference$beta_hat_os - z_alpha * inference$se_beta_os
      beta_ci_upper_os <- inference$beta_hat_os + z_alpha * inference$se_beta_os
      beta_truth <- sim_data$beta_true[seq_along(inference$beta_hat)]
      inference$beta_ci_lower_os <- beta_ci_lower_os
      inference$beta_ci_upper_os <- beta_ci_upper_os
      inference$beta_cover_os <- as.numeric(beta_truth >= beta_ci_lower_os & beta_truth <= beta_ci_upper_os)
      inference$beta_ci_width_os <- beta_ci_upper_os - beta_ci_lower_os

      lambda_truth <- inference$quant_times^2
      inference$lambda_lower_os <- cumh_bootstrap$os$pointwise$lower
      inference$lambda_upper_os <- cumh_bootstrap$os$pointwise$upper
      inference$lambda_cover_os <- as.numeric(lambda_truth >= inference$lambda_lower_os & lambda_truth <= inference$lambda_upper_os)
      inference$lambda_ci_width_os <- cumh_bootstrap$os$pointwise$width
      inference$scb_lower_os <- cumh_bootstrap$os$scb$lower_test
      inference$scb_upper_os <- cumh_bootstrap$os$scb$upper_test
      lambda_truth_grid_10_90 <- inference$scb_test_grid_10_90^2
      inference$scb_lower_os_10_90 <- stats::approx(
        inference$integration_grid, cumh_bootstrap$os$scb$lower,
        xout = inference$scb_test_grid_10_90, rule = 2, ties = "ordered"
      )$y
      inference$scb_upper_os_10_90 <- stats::approx(
        inference$integration_grid, cumh_bootstrap$os$scb$upper,
        xout = inference$scb_test_grid_10_90, rule = 2, ties = "ordered"
      )$y
      inference$scb_cover_os_10_90 <- as.numeric(all(
        lambda_truth_grid_10_90 >= inference$scb_lower_os_10_90 &
          lambda_truth_grid_10_90 <= inference$scb_upper_os_10_90
      ))
      inference$scb_width_os_mean_10_90 <- mean(
        inference$scb_upper_os_10_90 - inference$scb_lower_os_10_90,
        na.rm = TRUE
      )
      lambda_truth_grid_20_90 <- inference$scb_test_grid_20_90^2
      inference$scb_lower_os_20_90 <- stats::approx(
        inference$integration_grid, cumh_bootstrap$os$scb$lower,
        xout = inference$scb_test_grid_20_90, rule = 2, ties = "ordered"
      )$y
      inference$scb_upper_os_20_90 <- stats::approx(
        inference$integration_grid, cumh_bootstrap$os$scb$upper,
        xout = inference$scb_test_grid_20_90, rule = 2, ties = "ordered"
      )$y
      inference$scb_cover_os_20_90 <- as.numeric(all(
        lambda_truth_grid_20_90 >= inference$scb_lower_os_20_90 &
          lambda_truth_grid_20_90 <= inference$scb_upper_os_20_90
      ))
      inference$scb_width_os_mean_20_90 <- mean(
        inference$scb_upper_os_20_90 - inference$scb_lower_os_20_90,
        na.rm = TRUE
      )
      inference$scb_cover_os <- inference$scb_cover_os_10_90
      inference$scb_width_os_mean <- inference$scb_width_os_mean_10_90
      survival_truth <- exp(-(inference$quant_times^2) * exp(sum(sim_data$beta_true * x0)))
      inference$survival_lower_os <- survival_bootstrap$os$pointwise$lower
      inference$survival_upper_os <- survival_bootstrap$os$pointwise$upper
      inference$survival_cover_os <- as.numeric(survival_truth >= inference$survival_lower_os & survival_truth <= inference$survival_upper_os)
      inference$survival_ci_width_os <- survival_bootstrap$os$pointwise$width
      inference$survival_scb_lower_os <- survival_bootstrap$os$scb$lower_test
      inference$survival_scb_upper_os <- survival_bootstrap$os$scb$upper_test
      survival_truth_grid_10_90 <- exp(-(inference$scb_test_grid_10_90^2) * exp(sum(sim_data$beta_true * x0)))
      inference$survival_scb_lower_os_10_90 <- stats::approx(
        inference$integration_grid, survival_bootstrap$os$scb$lower,
        xout = inference$scb_test_grid_10_90, rule = 2, ties = "ordered"
      )$y
      inference$survival_scb_upper_os_10_90 <- stats::approx(
        inference$integration_grid, survival_bootstrap$os$scb$upper,
        xout = inference$scb_test_grid_10_90, rule = 2, ties = "ordered"
      )$y
      inference$survival_scb_cover_os_10_90 <- as.numeric(all(
        survival_truth_grid_10_90 >= inference$survival_scb_lower_os_10_90 &
          survival_truth_grid_10_90 <= inference$survival_scb_upper_os_10_90
      ))
      inference$survival_scb_width_os_mean_10_90 <- mean(
        inference$survival_scb_upper_os_10_90 - inference$survival_scb_lower_os_10_90,
        na.rm = TRUE
      )
      survival_truth_grid_20_90 <- exp(-(inference$scb_test_grid_20_90^2) * exp(sum(sim_data$beta_true * x0)))
      inference$survival_scb_lower_os_20_90 <- stats::approx(
        inference$integration_grid, survival_bootstrap$os$scb$lower,
        xout = inference$scb_test_grid_20_90, rule = 2, ties = "ordered"
      )$y
      inference$survival_scb_upper_os_20_90 <- stats::approx(
        inference$integration_grid, survival_bootstrap$os$scb$upper,
        xout = inference$scb_test_grid_20_90, rule = 2, ties = "ordered"
      )$y
      inference$survival_scb_cover_os_20_90 <- as.numeric(all(
        survival_truth_grid_20_90 >= inference$survival_scb_lower_os_20_90 &
          survival_truth_grid_20_90 <= inference$survival_scb_upper_os_20_90
      ))
      inference$survival_scb_width_os_mean_20_90 <- mean(
        inference$survival_scb_upper_os_20_90 - inference$survival_scb_lower_os_20_90,
        na.rm = TRUE
      )
      inference$survival_scb_cover_os <- inference$survival_scb_cover_os_10_90
      inference$survival_scb_width_os_mean <- inference$survival_scb_width_os_mean_10_90

      record_list[[rec_idx]] <- build_replication_records(
        inference = inference,
        beta_true = sim_data$beta_true,
        design_row = row,
        rep_id = rep_id,
        c_sieve = c_sieve_use,
        z_alpha = z_alpha,
        cumh_bootstrap = cumh_bootstrap,
        survival_bootstrap = survival_bootstrap
      )
      replication_list[[rec_idx]] <- list(
        design_index = d,
        replication = rep_id,
        inference = inference,
        cumh_bootstrap = cumh_bootstrap,
        survival_bootstrap = survival_bootstrap,
        lambda_delta_vec = by_column_lambda$lambda
      )
      rec_idx <- rec_idx + 1L
    }
  }

  if (length(record_list) == 0L) {
    stop("No successful replications were completed.")
  }

  records <- do.call(rbind, record_list)
  summary_df <- summarize_coverage_records(records)

  if (!is.null(output_csv)) {
    write.csv(summary_df, output_csv, row.names = FALSE)
  }

  list(summary = summary_df, records = records, replications = replication_list)
}
