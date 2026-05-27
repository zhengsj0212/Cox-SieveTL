#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(pracma)
  library(MASS)
})

# =========================================================
# Paths / setup
# =========================================================
setwd('/Users/yuxisong/Library/CloudStorage/Box-Box/Cox-SieveTL/simulation/results/new_icb_response_result')

source("./codes/cv_source_selection.R")
source("/Users/yuxisong/Library/CloudStorage/Box-Box/Cox-SieveTL/simulation/results/new_icb_response_result/codes/cv_eval_selected_source_bic.R")

results_dir   <- "./results_allTMB_new_hessian_lasso"
selected_file <- file.path(results_dir, "selected_sources_summary.csv")
data_file     <- "./codes/extracted_cancer_data_by_type.xlsx"

pilot_dir <- file.path(results_dir, "pilot_parameter_by_column_realdata")
out_dir   <- file.path(results_dir, "realdata_transfer_hessian_inference")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# Controls
# =========================================================
lambda_zeta_default <- 0.005
lambda_eta_default  <- 0
quant_probs <- seq(0.1, 0.9, by = 0.1)
integration_grid_size <- 500
multiplier_boot <- 500
conf_level <- 0.95
z_alpha <- qnorm(1 - (1 - conf_level) / 2)


# =========================================================
# Target list / lambda map
# =========================================================
cv_tbl <- read.csv("./results_allTMB_new_hessian_lasso/cv_metrics_summary_table_with_c_lambda_bic.csv",
                   stringsAsFactors = FALSE)

cv_tbl2 <- cv_tbl %>%
  mutate(
    best_lambda_zeta = ifelse(is.na(best_lambda_zeta) | best_lambda_zeta == 0,
                              lambda_zeta_default, best_lambda_zeta),
    best_lambda_eta  = ifelse(is.na(best_lambda_eta)  | best_lambda_eta  == 0,
                              lambda_eta_default,  best_lambda_eta)
  )
c_map <- cv_tbl$c_mult
all_sheets <- paste0("cancer", c(1:10, 12:16))
names(c_map) <- as.character(all_sheets)
c_default <- 1
lambda_map <- setNames(
  lapply(seq_len(nrow(cv_tbl2)), function(i) {
    list(
      lambda_zeta = cv_tbl2$best_lambda_zeta[i],
      lambda_eta  = cv_tbl2$best_lambda_eta[i]
    )
  }),
  cv_tbl2$target_cancer
)

target_sheets <- cv_tbl2$target_cancer

pretty_name_map <- c(
  cancer1  = "Bladder",
  cancer2  = "Breast",
  cancer3  = "Colorectal",
  cancer4  = "Endometrial",
  cancer5  = "Esophageal",
  cancer6  = "Gastric",
  cancer7  = "Head & Neck",
  cancer8  = "Hepatobiliary",
  cancer9  = "Melanoma",
  cancer10 = "Mesothelioma",
  cancer11 = "NSCLC",
  cancer12 = "Ovarian",
  cancer13 = "Pancreatic",
  cancer14 = "Renal",
  cancer15 = "Sarcoma",
  cancer16 = "SCLC",
  cancer18 = "CNS"
)

# =========================================================
# Helpers for transfer inverse-Hessian inference
# =========================================================
.smooth_abs <- function(u, eps) sqrt(u * u + eps * eps) - eps
.smooth_abs_grad <- function(u, eps) u / sqrt(u * u + eps * eps)

read_by_column_lambda_file <- function(path, d_n) {
  if (!file.exists(path)) {
    stop("pilot_selected_c_by_column.csv not found: ", path)
  }
  tbl <- read.csv(path, stringsAsFactors = FALSE)
  required_cols <- c("j", "best_c", "best_lambda_value")
  if (!all(required_cols %in% names(tbl))) {
    stop("Missing required columns in by-column lambda file: ", path)
  }
  tbl <- tbl[order(tbl$j), ]
  if (!identical(as.integer(tbl$j), seq_len(d_n))) {
    stop("Column indices do not match 1:d_n for file: ", path)
  }
  list(
    c = as.numeric(tbl$best_c),
    lambda = as.numeric(tbl$best_lambda_value)
  )
}

combine_domains <- function(target, source_list = NULL) {
  if (is.null(source_list) || length(source_list) == 0L) return(target)
  list(
    x = do.call(rbind, c(list(target$x), lapply(source_list, `[[`, "x"))),
    stime = c(target$stime, unlist(lapply(source_list, `[[`, "stime"))),
    type = c(target$type, unlist(lapply(source_list, `[[`, "type")))
  )
}

compute_target_score_theta <- function(data, gamma_hat, beta_hat, g_funcs, time_grid) {
  n <- nrow(data$x)
  L <- length(beta_hat)
  p <- length(gamma_hat)
  d <- p + L
  
  grad <- rep(0, d)
  
  for (i in seq_len(n)) {
    x_i <- data$x[i, ]
    t_i <- data$stime[i]
    delta_i <- data$type[i]
    
    if (!is.finite(t_i) || t_i <= 0) next
    
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
    
    integral <- apply(
      Z_mat * rep(exp_eta, each = d),
      1,
      function(row) trapz(t_grid_i, row)
    )
    
    grad <- grad + integral
  }
  
  grad / n
}

compute_target_score_theta_subjects <- function(data, gamma_hat, beta_hat, g_funcs, time_grid) {
  n <- nrow(data$x)
  L <- length(beta_hat)
  p <- length(gamma_hat)
  d <- p + L
  
  score_mat <- matrix(0, nrow = n, ncol = d)
  
  for (i in seq_len(n)) {
    x_i <- data$x[i, ]
    t_i <- data$stime[i]
    delta_i <- data$type[i]
    
    if (!is.finite(t_i) || t_i <= 0) next
    
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
      
      integral <- apply(
        Z_mat * rep(exp_eta, each = d),
        1,
        function(row) trapz(t_grid_i, row)
      )
      
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
  
  for (i in seq_len(n_k)) {
    x_i <- data_k$x[i, ]
    t_i <- data_k$stime[i]
    if (!is.finite(t_i) || t_i <= 0) next
    
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

fit_inverse_hessian_correction <- function(H_target, omega_A, j, lambda_l1,
                                           eps = 1e-4, ridge = 1e-6) {
  d <- ncol(H_target)
  e_j <- rep(0, d)
  e_j[j] <- 1
  
  if (lambda_l1 == 0) {
    return(
      tryCatch(
        as.vector(solve(H_target, e_j) - omega_A),
        error = function(e) {
          as.vector(solve(H_target + ridge * diag(d), e_j) - omega_A)
        }
      )
    )
  }
  
  obj <- function(delta) {
    omega <- omega_A + delta
    r <- as.vector(H_target %*% omega - e_j)
    
    0.5 * sum(r^2) +
      lambda_l1 * sum(.smooth_abs(delta, eps))
  }
  
  grad <- function(delta) {
    omega <- omega_A + delta
    r <- as.vector(H_target %*% omega - e_j)
    
    as.vector(
      t(H_target) %*% r +
        lambda_l1 * .smooth_abs_grad(delta, eps)
    )
  }
  
  fit <- optim(
    par = rep(0, d),
    fn = obj,
    gr = grad,
    method = "BFGS",
    control = list(maxit = 2000, reltol = 1e-10)
  )
  
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

compute_plugin_inference <- function(target, source, beta_hat, gamma_hat, g_funcs,
                                     lambda_delta_vec,
                                     quant_times = NULL,
                                     grid_range = NULL,
                                     quant_probs = seq(0.1, 0.9, by = 0.1),
                                     integration_grid_size = 500) {
  x <- as.matrix(target$x)
  stime <- as.numeric(target$stime)
  n0 <- nrow(x)
  L <- ncol(x)
  p <- length(gamma_hat)
  
  if (is.null(grid_range)) {
    grid_start <- min(stime, na.rm = TRUE)
    grid_end <- max(stime, na.rm = TRUE)
  } else {
    grid_start <- as.numeric(grid_range[1])
    grid_end <- as.numeric(grid_range[2])
  }
  
  scb_grid <- seq(grid_start, grid_end, length.out = max(as.integer(integration_grid_size), 2L))
  time_grid_target <- seq(
    min(stime, na.rm = TRUE),
    max(stime, na.rm = TRUE),
    length.out = max(as.integer(integration_grid_size), 2L)
  )
  
  pooled_stime <- c(target$stime, unlist(lapply(source, `[[`, "stime")))
  pooled_observed_times <- sort(unique(pooled_stime[is.finite(pooled_stime)]))
  if (length(pooled_observed_times) < 2L) {
    pooled_observed_times <- sort(unique(c(grid_start, grid_end)))
  }
  
  pooled_time_grid <- seq(
    min(pooled_stime, na.rm = TRUE),
    max(pooled_stime, na.rm = TRUE),
    length.out = max(as.integer(integration_grid_size), 2L)
  )
  
  # Hessian / score evaluated at final lasso-corrected estimator
  H_target <- compute_target_hessian_theta(
    data_k = target,
    gamma_hat = gamma_hat,
    beta_hat = beta_hat,
    g_funcs = g_funcs,
    time_grid = time_grid_target
  )
  
  H_pooled <- compute_pooled_hessian_theta(
    target = target,
    source_list = source,
    gamma_hat = gamma_hat,
    beta_hat = beta_hat,
    g_funcs = g_funcs,
    time_grid = pooled_time_grid
  )
  H_target <- regularize_hessian(H_target, ridge = 1e-6)
  H_pooled <- regularize_hessian(H_pooled, ridge = 1e-6)
  Omega_hat_by_column <- estimate_transfer_inverse_hessian_by_column(
    H_pooled = H_pooled,
    H_target = H_target,
    lambda_delta_vec = lambda_delta_vec
  )
  
  score_subject_mat <- compute_target_score_theta_subjects(
    data = target,
    gamma_hat = gamma_hat,
    beta_hat = beta_hat,
    g_funcs = g_funcs,
    time_grid = time_grid_target
  )
  
  # Influence function approximation
  IF_theta_hat_mat <- -score_subject_mat %*% t(Omega_hat_by_column)
  IF_gamma_hat_mat <- IF_theta_hat_mat[, seq_len(p), drop = FALSE]
  IF_beta_hat_mat  <- IF_theta_hat_mat[, (p + 1):(p + L), drop = FALSE]
  
  # IMPORTANT:
  # Point estimator is the lasso-corrected estimator from fit_algorithm1().
  # Do NOT apply another theta - Omega * score update here.
  gamma_os <- gamma_hat
  beta_os  <- beta_hat
  
  avar_beta <- crossprod(IF_beta_hat_mat) / n0
  avar_beta <- (avar_beta + t(avar_beta)) / 2
  se_beta <- sqrt(pmax(diag(avar_beta) / n0, 0))
  
  if (is.null(quant_times)) {
    quant_times <- as.numeric(stats::quantile(
      stime[is.finite(stime)],
      probs = quant_probs,
      type = 7,
      names = FALSE
    ))
  } else {
    quant_times <- as.numeric(quant_times)
  }
  
  quant_labels <- paste0("Lambda_q", sprintf("%02d", round(100 * quant_probs)))
  
  # Plug-in cumulative baseline hazard from final lasso-corrected gamma
  base_pooled_curve <- make_baseline_from_gamma(
    gamma_hat = gamma_hat,
    g_funcs = g_funcs,
    t_grid = pooled_observed_times
  )
  
  Lambda0_os_quant <- stats::approx(
    base_pooled_curve$time,
    base_pooled_curve$H0,
    xout = quant_times,
    rule = 2,
    ties = "ordered"
  )$y
  
  base_grid_curve <- make_baseline_from_gamma(
    gamma_hat = gamma_hat,
    g_funcs = g_funcs,
    t_grid = scb_grid
  )
  
  Lambda0_os_grid <- base_grid_curve$H0
  
  # IF for cumulative hazard
  a_pooled_obj <- compute_a_matrix(gamma_hat, g_funcs, pooled_observed_times)
  IF_Lambda_hat_pooled <- IF_gamma_hat_mat %*% t(a_pooled_obj$a_mat)
  
  IF_Lambda_hat_mat <- t(vapply(seq_len(nrow(IF_Lambda_hat_pooled)), function(i) {
    stats::approx(
      a_pooled_obj$time,
      IF_Lambda_hat_pooled[i, ],
      xout = quant_times,
      rule = 2,
      ties = "ordered"
    )$y
  }, numeric(length(quant_times))))
  
  a_grid_obj <- compute_a_matrix(gamma_hat, g_funcs, scb_grid)
  IF_Lambda_hat_grid <- IF_gamma_hat_mat %*% t(a_grid_obj$a_mat)
  
  var_hat_Lambda <- colSums(IF_Lambda_hat_mat^2) / n0
  se_Lambda <- sqrt(pmax(var_hat_Lambda / n0, 0))
  
  list(
    n0 = n0,
    theta_hat = c(gamma_hat, beta_hat),
    gamma_hat = gamma_hat,
    beta_hat = beta_hat,
    gamma_os = gamma_os,
    beta_os = beta_os,
    beta_hat_os = beta_os,
    se_beta = se_beta,
    Omega_hat_by_column = Omega_hat_by_column,
    IF_theta_hat_mat = IF_theta_hat_mat,
    IF_beta_hat_mat = IF_beta_hat_mat,
    IF_gamma_hat_mat = IF_gamma_hat_mat,
    Sigma_hat_beta = avar_beta,
    integration_grid = scb_grid,
    quant_times = quant_times,
    quant_labels = quant_labels,
    Lambda_hat_os_quant = Lambda0_os_quant,
    Lambda_hat_os_grid = Lambda0_os_grid,
    se_Lambda = se_Lambda,
    var_hat_Lambda = var_hat_Lambda,
    IF_Lambda_hat_mat = IF_Lambda_hat_mat,
    IF_Lambda_hat_grid = IF_Lambda_hat_grid
  )
}
compute_multiplier_cumh_intervals <- function(inference,
                                              alpha = 0.05,
                                              multiplier_boot = 500,
                                              multiplier = c("normal", "rademacher"),
                                              seed = NULL,
                                              sd_floor = 1e-8) {
  multiplier <- match.arg(multiplier)
  n0 <- inference$n0
  if (!is.null(seed)) set.seed(seed)
  
  if (multiplier == "normal") {
    xi <- matrix(rnorm(multiplier_boot * n0), nrow = multiplier_boot, ncol = n0)
  } else {
    xi <- matrix(sample(c(-1, 1), size = multiplier_boot * n0, replace = TRUE),
                 nrow = multiplier_boot, ncol = n0)
  }
  
  z_quant <- (xi %*% inference$IF_Lambda_hat_mat) / sqrt(n0)
  s_quant_boot <- sqrt(pmax(colMeans(z_quant^2), 0))
  se_quant_boot <- s_quant_boot / sqrt(n0)
  z_alpha <- stats::qnorm(1 - alpha / 2)
  ci_lower_quant <- inference$Lambda_hat_os_quant - z_alpha * se_quant_boot
  ci_upper_quant <- inference$Lambda_hat_os_quant + z_alpha * se_quant_boot
  
  z_grid <- (xi %*% inference$IF_Lambda_hat_grid) / sqrt(n0)
  s_grid_boot <- sqrt(pmax(colMeans(z_grid^2), 0))
  s_grid_safe <- pmax(s_grid_boot, sd_floor)
  studentized <- abs(sweep(z_grid, 2, s_grid_safe, "/"))
  sup_stats <- apply(studentized, 1, max, na.rm = TRUE)
  scb_crit <- as.numeric(stats::quantile(sup_stats, probs = 1 - alpha, names = FALSE, na.rm = TRUE))
  scb_halfwidth <- scb_crit * s_grid_safe / sqrt(n0)
  scb_lower_grid <- inference$Lambda_hat_os_grid - scb_halfwidth
  scb_upper_grid <- inference$Lambda_hat_os_grid + scb_halfwidth
  
  list(
    pointwise = list(
      lower = ci_lower_quant,
      upper = ci_upper_quant,
      se = se_quant_boot
    ),
    scb = list(
      s_grid = s_grid_boot,
      crit = scb_crit,
      lower = scb_lower_grid,
      upper = scb_upper_grid
    )
  )
}

# =========================================================
# Main function for one target
# =========================================================
run_realdata_transfer_hessian_inference_one <- function(tname) {
  
  cancers <- lapply(all_sheets, function(sh) {
    df <- read_excel(data_file, sheet = sh)
    df <- cap_cancer_df(df)
    na.omit(df)
  })
  names(cancers) <- all_sheets
  proc <- lapply(cancers, prep_df)
  
  selected_df <- read.csv(selected_file, stringsAsFactors = FALSE)
  
  target <- proc[[tname]]
  target_df_raw <- target$df
  
  lam <- lambda_map[[tname]]
  if (is.null(lam)) {
    lam <- list(lambda_zeta = lambda_zeta_default, lambda_eta = lambda_eta_default)
  }
  
  lambda_zeta_t <- lam$lambda_zeta
  lambda_eta_t  <- lam$lambda_eta
  c_mult <- if (tname %in% names(c_map)) c_map[[tname]] else c_default
  
  sel_row <- selected_df[selected_df$target_cancer == tname, , drop = FALSE]
  sel_sources <- character(0)
  if (nrow(sel_row) > 0 &&
      !is.na(sel_row$n_selected[1]) && sel_row$n_selected[1] > 0) {
    sel_sources <- parse_sources(sel_row$selected_sources[1])
    sel_sources <- sel_sources[sel_sources != "" & !is.na(sel_sources)]
    sel_sources <- sel_sources[sel_sources %in% names(proc) & sel_sources != tname]
  }
  
  if (length(sel_sources) == 0) {
    stop("No selected sources found for ", tname)
  }
  
  source_df_pool_raw <- bind_rows(lapply(sel_sources, function(sn) proc[[sn]]$df))
  
  apply_list <- list(src = source_df_pool_raw)
  pre <- preprocess_features(
    target_df_raw,
    apply_list,
    target$feature_cols,
    target$time_col,
    target$event_col
  )
  
  target_df <- pre$train_df
  source_df_pool <- pre$apply_df_list$src
  
  # fit your real-data SieveTL
  fit <- fit_algorithm1(
    target_train_df = target_df,
    source_df = source_df_pool,
    feature_cols = target$feature_cols,
    time_col = target$time_col,
    event_col = target$event_col,
    lambda_zeta = lambda_zeta_t,
    lambda_eta = lambda_eta_t,
    c_mult = c_mult,
    maxit = 5000,
    use_hessian_lasso_update = TRUE,
    ridge_iter = 1e-4
  )
  
  conv_ok <- is.list(fit) &&
    !is.null(fit$conv_transfer) && fit$conv_transfer == 1L &&
    !is.null(fit$conv_debias) && fit$conv_debias == 1L &&
    !is.null(fit$gamma_hat) && !is.null(fit$beta_hat) &&
    !is.null(fit$basis)
  
  if (!conv_ok) stop("SieveTL fit failed for ", tname)
  
  target_list <- list(
    x = as.matrix(target_df[, target$feature_cols, drop = FALSE]),
    stime = as.numeric(target_df[[target$time_col]]),
    type = as.numeric(target_df[[target$event_col]])
  )
  
  source_list <- list(
    list(
      x = as.matrix(source_df_pool[, target$feature_cols, drop = FALSE]),
      stime = as.numeric(source_df_pool[[target$time_col]]),
      type = as.numeric(source_df_pool[[target$event_col]])
    )
  )
  
  all_times <- c(target_list$stime, source_list[[1]]$stime)
  G0 <- fit$basis$g_eval(all_times)
  if (is.null(dim(G0))) G0 <- matrix(G0, ncol = 1)
  p_hat <- ncol(G0)
  
  g_funcs <- lapply(seq_len(p_hat), function(j) {
    force(j)
    function(t) {
      out <- fit$basis$g_eval(t)
      if (is.null(dim(out))) out <- matrix(out, ncol = p_hat)
      as.numeric(out[, j])
    }
  })
  
  L <- ncol(target_list$x)
  d_n <- p_hat + L
  
  lambda_file <- file.path(pilot_dir, tname, "pilot_selected_c_by_column.csv")
  by_column_lambda <- read_by_column_lambda_file(lambda_file, d_n = d_n)
  lambda_delta_vec <- by_column_lambda$lambda
  
  pooled_observed_times <- c(
    target_list$stime[is.finite(target_list$stime)],
    source_list[[1]]$stime[is.finite(source_list[[1]]$stime)]
  )
  
  quant_times <- as.numeric(
    stats::quantile(pooled_observed_times, probs = quant_probs, type = 7, names = FALSE)
  )
  grid_range <- range(pooled_observed_times)
  
  inference <- compute_plugin_inference(
    target = target_list,
    source = source_list,
    beta_hat = fit$beta_hat,
    gamma_hat = fit$gamma_hat,
    g_funcs = g_funcs,
    lambda_delta_vec = lambda_delta_vec,
    quant_times = quant_times,
    grid_range = grid_range,
    quant_probs = quant_probs,
    integration_grid_size = integration_grid_size
  )
  
  boot_out <- compute_multiplier_cumh_intervals(
    inference = inference,
    alpha = 1 - conf_level,
    multiplier_boot = multiplier_boot,
    multiplier = "normal",
    seed = 2026
  )
  
  # -------------------------------------------------------
  # Beta results: estimate / SE / Wald CI / z / p
  # -------------------------------------------------------
  beta_results <- data.frame(
    target_cancer = tname,
    target_label = ifelse(tname %in% names(pretty_name_map), pretty_name_map[[tname]], tname),
    covariate = target$feature_cols,
    estimate = inference$beta_hat_os,
    se = inference$se_beta,
    z = inference$beta_hat_os / inference$se_beta,
    p_value = 2 * pnorm(-abs(inference$beta_hat_os / inference$se_beta)),
    ci_lower = inference$beta_hat_os - z_alpha * inference$se_beta,
    ci_upper = inference$beta_hat_os + z_alpha * inference$se_beta,
    stringsAsFactors = FALSE
  )
  
  # -------------------------------------------------------
  # Cumhaz pointwise Wald CI / z / p against 0
  # -------------------------------------------------------
  cumhaz_wald <- data.frame(
    target_cancer = tname,
    target_label = ifelse(tname %in% names(pretty_name_map), pretty_name_map[[tname]], tname),
    quantile = inference$quant_labels,
    time = inference$quant_times,
    estimate = inference$Lambda_hat_os_quant,
    se = inference$se_Lambda,
    z = inference$Lambda_hat_os_quant / inference$se_Lambda,
    p_value = 2 * pnorm(-abs(inference$Lambda_hat_os_quant / inference$se_Lambda)),
    ci_lower = inference$Lambda_hat_os_quant - z_alpha * inference$se_Lambda,
    ci_upper = inference$Lambda_hat_os_quant + z_alpha * inference$se_Lambda,
    stringsAsFactors = FALSE
  )
  
  # -------------------------------------------------------
  # Cumhaz multiplier pointwise CI
  # -------------------------------------------------------
  cumhaz_boot <- data.frame(
    target_cancer = tname,
    target_label = ifelse(tname %in% names(pretty_name_map), pretty_name_map[[tname]], tname),
    quantile = inference$quant_labels,
    time = inference$quant_times,
    estimate = inference$Lambda_hat_os_quant,
    se_boot = boot_out$pointwise$se,
    ci_lower = boot_out$pointwise$lower,
    ci_upper = boot_out$pointwise$upper,
    stringsAsFactors = FALSE
  )
  
  # -------------------------------------------------------
  # Cumhaz simultaneous confidence band
  # -------------------------------------------------------
  cumhaz_scb <- data.frame(
    target_cancer = tname,
    target_label = ifelse(tname %in% names(pretty_name_map), pretty_name_map[[tname]], tname),
    time = inference$integration_grid,
    estimate = inference$Lambda_hat_os_grid,
    scb_lower = boot_out$scb$lower,
    scb_upper = boot_out$scb$upper,
    scb_crit = boot_out$scb$crit,
    stringsAsFactors = FALSE
  )
  
  write.csv(beta_results,
            file.path(out_dir, paste0("beta_wald_ci_", tname, ".csv")),
            row.names = FALSE)
  
  write.csv(cumhaz_wald,
            file.path(out_dir, paste0("cumhaz_wald_ci_", tname, ".csv")),
            row.names = FALSE)
  
  write.csv(cumhaz_boot,
            file.path(out_dir, paste0("cumhaz_boot_pointwise_ci_", tname, ".csv")),
            row.names = FALSE)
  
  write.csv(cumhaz_scb,
            file.path(out_dir, paste0("cumhaz_scb_", tname, ".csv")),
            row.names = FALSE)
  
  list(
    beta = beta_results,
    cumhaz_wald = cumhaz_wald,
    cumhaz_boot = cumhaz_boot,
    cumhaz_scb = cumhaz_scb,
    inference = inference,
    bootstrap = boot_out
  )
}

# =========================================================
# Run all targets and save combined tables
# =========================================================
run_all_targets_transfer_hessian_inference <- function(target_vec = target_sheets) {
  beta_list <- list()
  ch_wald_list <- list()
  ch_boot_list <- list()
  ch_scb_list <- list()
  
  for (tt in target_vec) {
    message(sprintf("[CI] running target = %s", tt))
    out <- tryCatch(
      run_realdata_transfer_hessian_inference_one(tt),
      error = function(e) {
        message(sprintf("[CI] target=%s failed: %s", tt, conditionMessage(e)))
        NULL
      }
    )
    
    if (!is.null(out)) {
      beta_list[[length(beta_list) + 1]] <- out$beta
      ch_wald_list[[length(ch_wald_list) + 1]] <- out$cumhaz_wald
      ch_boot_list[[length(ch_boot_list) + 1]] <- out$cumhaz_boot
      ch_scb_list[[length(ch_scb_list) + 1]] <- out$cumhaz_scb
    }
  }
  
  beta_all <- bind_rows(beta_list)
  ch_wald_all <- bind_rows(ch_wald_list)
  ch_boot_all <- bind_rows(ch_boot_list)
  ch_scb_all <- bind_rows(ch_scb_list)
  
  write.csv(beta_all,
            file.path(out_dir, "all_targets_beta_wald_ci.csv"),
            row.names = FALSE)
  
  write.csv(ch_wald_all,
            file.path(out_dir, "all_targets_cumhaz_wald_ci.csv"),
            row.names = FALSE)
  
  write.csv(ch_boot_all,
            file.path(out_dir, "all_targets_cumhaz_boot_pointwise_ci.csv"),
            row.names = FALSE)
  
  write.csv(ch_scb_all,
            file.path(out_dir, "all_targets_cumhaz_scb.csv"),
            row.names = FALSE)
  
  list(
    beta = beta_all,
    cumhaz_wald = ch_wald_all,
    cumhaz_boot = ch_boot_all,
    cumhaz_scb = ch_scb_all
  )
}

if (sys.nframe() == 0) {
  run_all_targets_transfer_hessian_inference()
}
