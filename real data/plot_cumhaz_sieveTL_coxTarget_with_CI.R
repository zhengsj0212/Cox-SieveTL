#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(pracma)
  library(survival)
  library(ggplot2)
  library(patchwork)
})

setwd('/Users/yuxisong/Library/CloudStorage/Box-Box/Cox-SieveTL/simulation/results/new_icb_response_result')

source("./codes/cv_source_selection.R")
source("./codes/cv_eval_selected_source_bic.R")

results_dir    <- "./results_allTMB_new_hessian_lasso_new"
data_file      <- "./codes/extracted_cancer_data_by_type.xlsx"
selected_file  <- file.path(results_dir, "selected_sources_summary.csv")
cv_lambda_path <- file.path(results_dir, "cv_metrics_summary_table_with_c_lambda_bic.csv")
pilot_dir      <- file.path(results_dir, "pilot_parameter_by_column_realdata")
plot_dir       <- file.path(results_dir, "cumh_sievetl_vs_coxph_target")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

lambda_zeta_default <- 0.005
lambda_eta_default  <- 0
c_mult_default <- 1
run_targets <- NULL

quant_probs <- seq(0.1, 0.9, by = 0.1)
integration_grid_size <- 500
multiplier_boot <- 500
conf_level <- 0.95
z_alpha <- qnorm(1 - (1 - conf_level) / 2)

cv_tbl <- read.csv("./results_allTMB_new_hessian_lasso_new/cv_metrics_summary_table_with_c_lambda_bic.csv", stringsAsFactors = FALSE)
cv_tbl2 <- cv_tbl %>%
  mutate(
    best_lambda_zeta = ifelse(is.na(best_lambda_zeta) | best_lambda_zeta == 0, lambda_zeta_default, best_lambda_zeta),
    best_lambda_eta  = ifelse(is.na(best_lambda_eta)  | best_lambda_eta  == 0, lambda_eta_default,  best_lambda_eta)
  )


c_map <- cv_tbl$c_mult
all_sheets <- paste0("cancer", c(1:10, 12:16))
names(c_map) <- as.character(all_sheets)

title_map <- c(
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

get_c_mult <- function(tname, default = 1) {
  if (tname %in% names(c_map)) c_map[[tname]] else default
}

parse_sources <- function(s) {
  if (is.null(s) || length(s) == 0 || is.na(s) || trimws(s) == "") return(character(0))
  trimws(strsplit(s, ",")[[1]])
}

pool_sources_df <- function(sel_sources, proc) {
  if (length(sel_sources) == 0) return(NULL)
  bind_rows(lapply(sel_sources, function(sn) proc[[sn]]$df))
}

make_g_funcs_from_fit <- function(fit) {
  p <- length(as.numeric(fit$gamma_hat))
  lapply(seq_len(p), function(j) {
    force(j)
    function(t) {
      gt <- as.matrix(fit$g_eval(t))
      as.numeric(gt[, j])
    }
  })
}

make_baseline_from_gamma <- function(gamma_hat, g_funcs, t_grid) {
  gamma_hat <- as.numeric(gamma_hat)
  t_grid <- sort(unique(as.numeric(t_grid)))
  
  G <- sapply(g_funcs, function(gf) as.numeric(gf(t_grid)))
  G <- as.matrix(G)
  if (nrow(G) != length(t_grid)) G <- t(G)
  
  logh0 <- drop(G %*% gamma_hat)
  h0 <- exp(logh0)
  dt <- diff(t_grid)
  H0 <- c(0, cumsum(dt * (head(h0, -1) + tail(h0, -1)) / 2))
  
  data.frame(time = t_grid, logh0 = logh0, h0 = h0, H0 = H0)
}

.smooth_abs <- function(u, eps) sqrt(u * u + eps * eps) - eps
.smooth_abs_grad <- function(u, eps) u / sqrt(u * u + eps * eps)

combine_domains <- function(target, source_list = NULL) {
  if (is.null(source_list) || length(source_list) == 0L) return(target)
  list(
    x = do.call(rbind, c(list(target$x), lapply(source_list, `[[`, "x"))),
    stime = c(target$stime, unlist(lapply(source_list, `[[`, "stime"))),
    type = c(target$type, unlist(lapply(source_list, `[[`, "type")))
  )
}

regularize_hessian <- function(H, ridge = 1e-6) {
  H <- (H + t(H)) / 2
  H + ridge * diag(ncol(H))
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
    
    exp_eta <- exp(as.vector(gamma_hat %*% G_mat) + sum(beta_hat * x_i))
    Z_mat <- rbind(G_mat, matrix(rep(x_i, times = m), nrow = L))
    
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
      
      exp_eta <- exp(as.vector(gamma_hat %*% G_mat) + sum(beta_hat * x_i))
      Z_mat <- rbind(G_mat, matrix(rep(x_i, times = m), nrow = L))
      
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
  H <- matrix(0, d, d)
  
  for (i in seq_len(n_k)) {
    x_i <- data_k$x[i, ]
    t_i <- data_k$stime[i]
    if (!is.finite(t_i) || t_i <= 0) next
    
    t_grid_i <- time_grid[time_grid <= t_i]
    if (length(t_grid_i) < 2L) next
    
    G_mat <- t(sapply(g_funcs, function(f) f(t_grid_i)))
    m <- length(t_grid_i)
    
    exp_eta_g <- exp(as.vector(gamma_hat %*% G_mat))
    exp_eta_x <- exp(sum(beta_hat * x_i))
    
    arr <- array(0, dim = c(d, d, m))
    for (j in seq_len(m)) {
      z_j <- c(G_mat[, j], x_i)
      arr[, , j] <- tcrossprod(z_j) * exp_eta_g[j]
    }
    
    H <- H + exp_eta_x * apply(arr, c(1, 2), function(v) trapz(t_grid_i, v))
  }
  
  H / n_k
}

compute_pooled_hessian_theta <- function(target, source_list, gamma_hat, beta_hat, g_funcs, time_grid) {
  pooled <- combine_domains(target, source_list)
  compute_target_hessian_theta(pooled, gamma_hat, beta_hat, g_funcs, time_grid)
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
          message(sprintf(
            "[Pilot] solve(H_target, e_%d) failed; using ridge=%g",
            j, ridge
          ))
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
estimate_transfer_inverse_hessian_by_column <- function(H_pooled, H_target, lambda_delta_vec) {
  d <- ncol(H_target)
  Omega <- matrix(0, d, d)
  
  for (j in seq_len(d)) {
    omega_A <- fit_inverse_hessian_initializer_unpenalized(H_pooled, j)
    delta_j <- fit_inverse_hessian_correction(H_target, omega_A, j, lambda_delta_vec[j])
    Omega[, j] <- omega_A + delta_j
  }
  
  Omega
}

read_by_column_lambda_file <- function(path, d_n) {
  if (!file.exists(path)) stop("Missing pilot file: ", path)
  tbl <- read.csv(path, stringsAsFactors = FALSE)
  tbl <- tbl[order(tbl$j), ]
  if (!identical(as.integer(tbl$j), seq_len(d_n))) {
    stop("Column index mismatch in ", path)
  }
  as.numeric(tbl$best_lambda_value)
}

compute_a_matrix <- function(gamma_hat, g_funcs, eval_times) {
  base <- make_baseline_from_gamma(gamma_hat, g_funcs, eval_times)
  p <- length(gamma_hat)
  G <- sapply(g_funcs, function(gf) as.numeric(gf(base$time)))
  G <- as.matrix(G)
  h0 <- base$h0
  dt <- diff(base$time)
  
  a_mat <- matrix(0, nrow = nrow(base), ncol = p)
  for (j in seq_len(p)) {
    integrand <- G[, j] * h0
    a_mat[, j] <- c(0, cumsum(dt * (head(integrand, -1) + tail(integrand, -1)) / 2))
  }
  
  list(time = base$time, a_mat = a_mat, base_curve = base)
}

compute_plugin_inference <- function(target, source, beta_hat, gamma_hat, g_funcs,
                                     lambda_delta_vec,
                                     quant_times,
                                     grid_range,
                                     integration_grid_size = 500) {
  n0 <- nrow(target$x)
  L <- length(beta_hat)
  p <- length(gamma_hat)
  
  time_grid_target <- seq(
    min(target$stime, na.rm = TRUE),
    max(target$stime, na.rm = TRUE),
    length.out = integration_grid_size
  )
  
  pooled_stime <- c(target$stime, unlist(lapply(source, `[[`, "stime")))
  pooled_time_grid <- seq(
    min(pooled_stime, na.rm = TRUE),
    max(pooled_stime, na.rm = TRUE),
    length.out = integration_grid_size
  )
  
  scb_grid <- seq(grid_range[1], grid_range[2], length.out = integration_grid_size)
  
  # -------------------------------------------------------
  # Hessian matrices
  # -------------------------------------------------------
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
  
  # -------------------------------------------------------
  # Transferred inverse Hessian
  # -------------------------------------------------------
  Omega <- estimate_transfer_inverse_hessian_by_column(
    H_pooled = H_pooled,
    H_target = H_target,
    lambda_delta_vec = lambda_delta_vec
  )
  
  # Symmetrize because Omega is estimated column-by-column
  Omega <- (Omega + t(Omega)) / 2
  
  # -------------------------------------------------------
  # IMPORTANT:
  # Use final Cox-SieveTL estimator directly.
  # Do NOT apply theta - Omega * score update here.
  # -------------------------------------------------------
  gamma_os <- as.numeric(gamma_hat)
  beta_os  <- as.numeric(beta_hat)
  
  # -------------------------------------------------------
  # Direct Hessian-based variance:
  # Var(theta_hat) approx Omega / n0
  # -------------------------------------------------------
  Sigma_gamma <- Omega[seq_len(p), seq_len(p), drop = FALSE]
  Sigma_beta  <- Omega[(p + 1):(p + L), (p + 1):(p + L), drop = FALSE]
  
  se_beta <- sqrt(pmax(diag(Sigma_beta) / n0, 0))
  
  # -------------------------------------------------------
  # Cumulative hazard delta-method variance
  # Lambda0(t) = int_0^t exp{gamma^T G(u)} du
  # d Lambda0(t) / d gamma = A(t)
  # Var{Lambda0(t)} approx A(t)^T Sigma_gamma A(t) / n0
  # -------------------------------------------------------
  a_quant <- compute_a_matrix(gamma_os, g_funcs, quant_times)
  A_quant <- a_quant$a_mat
  
  var_Lambda <- diag(A_quant %*% Sigma_gamma %*% t(A_quant))
  se_Lambda <- sqrt(pmax(var_Lambda / n0, 0))
  
  a_grid <- compute_a_matrix(gamma_os, g_funcs, scb_grid)
  A_grid <- a_grid$a_mat
  
  cov_Lambda_grid <- A_grid %*% Sigma_gamma %*% t(A_grid) / n0
  cov_Lambda_grid <- (cov_Lambda_grid + t(cov_Lambda_grid)) / 2
  
  list(
    n0 = n0,
    gamma_os = gamma_os,
    beta_os = beta_os,
    se_beta = se_beta,
    quant_times = quant_times,
    integration_grid = scb_grid,
    se_Lambda = se_Lambda,
    Omega_hat = Omega,
    Sigma_gamma = Sigma_gamma,
    Sigma_beta = Sigma_beta,
    cov_Lambda_grid = cov_Lambda_grid
  )
}

compute_multiplier_cumh_intervals <- function(inference,
                                              alpha = 0.05,
                                              multiplier_boot = 500,
                                              seed = 2026,
                                              sd_floor = 1e-8) {
  set.seed(seed)
  
  Sigma_grid <- inference$cov_Lambda_grid
  Sigma_grid <- (Sigma_grid + t(Sigma_grid)) / 2
  
  # Eigenvalue truncation to avoid numerical non-PSD issues
  eig <- eigen(Sigma_grid, symmetric = TRUE)
  vals <- pmax(eig$values, 0)
  
  sqrt_Sigma <- eig$vectors %*%
    diag(sqrt(vals), nrow = length(vals)) %*%
    t(eig$vectors)
  
  Z <- matrix(
    rnorm(multiplier_boot * nrow(Sigma_grid)),
    nrow = multiplier_boot,
    ncol = nrow(Sigma_grid)
  )
  
  gp_draws <- Z %*% sqrt_Sigma
  
  s_grid <- sqrt(pmax(diag(Sigma_grid), 0))
  s_grid_safe <- pmax(s_grid, sd_floor)
  
  studentized <- abs(sweep(gp_draws, 2, s_grid_safe, "/"))
  sup_stats <- apply(studentized, 1, max, na.rm = TRUE)
  
  crit <- as.numeric(
    quantile(sup_stats, probs = 1 - alpha, na.rm = TRUE, names = FALSE)
  )
  
  list(
    s_grid = s_grid,
    crit = crit
  )
}
compute_sievetl_cumh_real <- function(target_df, source_df, feature_cols, time_col, event_col,
                                      lambda_zeta, lambda_eta, c_mult,
                                      tname, pilot_dir) {
  fit <- fit_algorithm1(
    target_df,
    source_df = source_df,
    feature_cols = feature_cols,
    time_col = time_col,
    event_col = event_col,
    lambda_zeta = lambda_zeta,
    lambda_eta = lambda_eta,
    c_mult = c_mult,
    use_hessian_lasso_update = TRUE,
    ridge_iter = 1e-4
  )
  
  ok <- is.list(fit) &&
    !is.null(fit$conv_transfer) && fit$conv_transfer == 1L &&
    !is.null(fit$conv_debias) && fit$conv_debias == 1L &&
    !is.null(fit$gamma_hat) && !is.null(fit$beta_hat)
  
  if (!ok) return(NULL)
  
  g_funcs <- make_g_funcs_from_fit(fit)
  
  target_list <- list(
    x = as.matrix(target_df[, feature_cols, drop = FALSE]),
    stime = as.numeric(target_df[[time_col]]),
    type = as.numeric(target_df[[event_col]])
  )
  
  source_list <- list()
  if (!is.null(source_df) && nrow(source_df) > 0) {
    source_list[[1]] <- list(
      x = as.matrix(source_df[, feature_cols, drop = FALSE]),
      stime = as.numeric(source_df[[time_col]]),
      type = as.numeric(source_df[[event_col]])
    )
  }
  
  combined_times <- c(target_list$stime, unlist(lapply(source_list, `[[`, "stime")))
  combined_times <- combined_times[is.finite(combined_times)]
  
  p_hat <- length(fit$gamma_hat)
  L <- length(fit$beta_hat)
  d_n <- p_hat + L
  
  lambda_file <- file.path(pilot_dir, tname, "pilot_selected_c_by_column.csv")
  lambda_delta_vec <- read_by_column_lambda_file(lambda_file, d_n)
  
  plot_times <- seq(min(combined_times), max(combined_times), length.out = 400)
  grid_range <- range(combined_times)
  
  inference <- compute_plugin_inference(
    target = target_list,
    source = source_list,
    beta_hat = fit$beta_hat,
    gamma_hat = fit$gamma_hat,
    g_funcs = g_funcs,
    lambda_delta_vec = lambda_delta_vec,
    quant_times = plot_times,
    grid_range = grid_range,
    integration_grid_size = integration_grid_size
  )
  
  boot_out <- compute_multiplier_cumh_intervals(
    inference,
    alpha = 1 - conf_level,
    multiplier_boot = multiplier_boot
  )
  
  base_quant <- make_baseline_from_gamma(inference$gamma_os, g_funcs, plot_times)
  base_grid  <- make_baseline_from_gamma(inference$gamma_os, g_funcs, inference$integration_grid)
  
  curve_wald <- data.frame(
    time = base_quant$time,
    method = "Cox-SieveTL",
    ci_type = "wald",
    cumh = base_quant$H0,
    se = inference$se_Lambda,
    ci_lo = pmax(base_quant$H0 - z_alpha * inference$se_Lambda, 0),
    ci_hi = base_quant$H0 + z_alpha * inference$se_Lambda
  )
  
  scb_halfwidth <- boot_out$crit * pmax(boot_out$s_grid, 1e-8)
  
  scb_curve <- data.frame(
    time = base_grid$time,
    method = "Cox-SieveTL",
    ci_type = "scb",
    cumh = base_grid$H0,
    ci_lo = pmax(base_grid$H0 - scb_halfwidth, 0),
    ci_hi = base_grid$H0 + scb_halfwidth,
    scb_crit = boot_out$crit
  )
  
  list(
    fit = fit,
    inference = inference,
    curve_wald = curve_wald,
    scb_curve = scb_curve,
    combined_time_max = max(combined_times, na.rm = TRUE)
  )
}

compute_coxph_target_cumh_real <- function(target_df, feature_cols, time_col, event_col,
                                           quant_times) {
  xmat <- as.matrix(target_df[, feature_cols, drop = FALSE])
  df <- data.frame(
    stime = as.numeric(target_df[[time_col]]),
    type = as.numeric(target_df[[event_col]]),
    xmat,
    check.names = FALSE
  )
  names(df)[3:ncol(df)] <- paste0("x", seq_len(length(feature_cols)))
  fml <- as.formula(
    paste("survival::Surv(stime, type) ~", paste(names(df)[3:ncol(df)], collapse = " + "))
  )
  
  fit <- tryCatch(
    survival::coxph(fml, data = df, ties = "breslow", x = TRUE),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  
  zero_newdata <- as.data.frame(matrix(0, nrow = 1, ncol = length(feature_cols)))
  names(zero_newdata) <- paste0("x", seq_len(length(feature_cols)))
  
  target_time_max <- max(df$stime[is.finite(df$stime)], na.rm = TRUE)
  quant_times_use <- quant_times[quant_times <= target_time_max]
  if (length(quant_times_use) < 2) return(NULL)
  
  sf <- tryCatch(
    survival::survfit(fit, newdata = zero_newdata, se.fit = TRUE, conf.int = conf_level),
    error = function(e) NULL
  )
  if (is.null(sf)) return(NULL)
  
  sf_sum <- summary(sf, times = quant_times_use, extend = TRUE)
  
  cumhaz <- as.numeric(sf_sum$cumhaz)
  se_ch  <- as.numeric(sf_sum$std.chaz)
  
  curve <- data.frame(
    time = quant_times_use,
    method = "CoxPH (Target)",
    ci_type = "wald",
    cumh = cumhaz,
    se = se_ch,
    ci_lo = pmax(cumhaz - z_alpha * se_ch, 0),
    ci_hi = cumhaz + z_alpha * se_ch
  )
  
  list(
    fit = fit,
    curve = curve,
    target_time_max = target_time_max
  )
}

get_size_map_cumh <- function(proc, selected_df) {
  bind_rows(lapply(names(proc), function(tname) {
    target_info <- proc[[tname]]
    target_df0  <- target_info$df
    
    sel_row <- selected_df[selected_df$target_cancer == tname, , drop = FALSE]
    sel_sources <- character(0)
    if (nrow(sel_row) > 0 &&
        !is.na(sel_row$n_selected[1]) &&
        sel_row$n_selected[1] > 0) {
      sel_sources <- parse_sources(sel_row$selected_sources[1])
      sel_sources <- sel_sources[sel_sources %in% names(proc) & sel_sources != tname]
    }
    
    source_df0 <- pool_sources_df(sel_sources, proc)
    
    apply_list0 <- list()
    if (!is.null(source_df0)) apply_list0$src <- source_df0
    
    pre0 <- preprocess_features(
      target_df0,
      apply_list0,
      target_info$feature_cols,
      target_info$time_col,
      target_info$event_col
    )
    
    tibble(
      target_cancer = tname,
      n_target = nrow(pre0$train_df),
      n_source = if (!is.null(source_df0)) nrow(pre0$apply_df_list$src) else 0
    )
  }))
}

plot_sievetl_vs_coxph_cumh <- function(df_plot, df_scb = NULL,
                                       title_text = NULL,
                                       n_target = NA,
                                       n_source = NA,
                                       x_max = NULL,
                                       y_max = NULL) {
  fill_vals <- c(
    "Cox-SieveTL" = "#E64B35",
    "CoxPH (Target)" = "#4DBBD5",
    "Cox-SieveTL (SCB)" = "grey75"
  )
  
  color_vals <- c(
    "Cox-SieveTL" = "#E64B35",
    "CoxPH (Target)" = "#4DBBD5"
  )
  
  p <- ggplot()
  
  if (!is.null(df_scb) && nrow(df_scb) > 0) {
    df_scb <- df_scb %>%
      mutate(legend_fill = "Cox-SieveTL (SCB)")
    
    p <- p +
      geom_ribbon(
        data = df_scb,
        aes(x = time, ymin = ci_lo, ymax = ci_hi, fill = legend_fill),
        alpha = 0.22,
        color = NA
      )
  }
  
  df_ribbon <- df_plot %>%
    mutate(legend_fill = method) %>%
    filter(method %in% c("Cox-SieveTL", "CoxPH (Target)"))
  
  p <- p +
    geom_ribbon(
      data = df_ribbon,
      aes(x = time, ymin = ci_lo, ymax = ci_hi, fill = legend_fill),
      alpha = 0.18,
      color = NA
    ) +
    geom_line(
      data = df_ribbon,
      aes(x = time, y = cumh, color = method),
      linewidth = 1
    ) +
    scale_fill_manual(
      name = NULL,
      values = fill_vals,
      breaks = c("Cox-SieveTL", "CoxPH (Target)", "Cox-SieveTL (SCB)")
    ) +
    scale_color_manual(
      values = color_vals,
      guide = "none"
    ) +
    labs(
      x = "Time",
      y = expression(hat(Lambda)[0](t)),
      title = title_text
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
      legend.position = "bottom",
      legend.text = element_text(size = 12)
    )
  
  xlim_use <- if (!is.null(x_max) && is.finite(x_max)) c(0, x_max) else NULL
  ylim_use <- if (!is.null(y_max) && is.finite(y_max)) c(0, y_max) else NULL
  
  if (!is.null(xlim_use) || !is.null(ylim_use)) {
    p <- p + coord_cartesian(xlim = xlim_use, ylim = ylim_use)
  }
  
  x_rng <- if (!is.null(xlim_use)) xlim_use else range(df_plot$time, na.rm = TRUE)
  y_all <- c(df_plot$ci_hi, df_plot$cumh)
  if (!is.null(df_scb) && nrow(df_scb) > 0) y_all <- c(y_all, df_scb$ci_hi, df_scb$cumh)
  y_rng <- if (!is.null(ylim_use)) ylim_use else range(y_all, na.rm = TRUE)
  
  label_txt <- paste0("nT=", n_target, "\n", "nS=", n_source)
  
  p +
    annotate(
      "label",
      x = x_rng[1] + 0.015 * diff(x_rng),
      y = y_rng[2] - 0.02 * diff(y_rng),
      label = label_txt,
      hjust = 0,
      vjust = 1,
      size = 5,
      label.size = 0.6,
      fill = "white",
      color = "black"
    )
}

# =========================================================
# LOAD DATA
# =========================================================
all_sheets <- paste0("cancer", c(1:10, 12:16, 18))

cancers <- lapply(all_sheets, function(sh) {
  df <- read_excel(data_file, sheet = sh)
  df <- cap_cancer_df(df)
  df %>% drop_na()
})
names(cancers) <- all_sheets

proc <- lapply(cancers, prep_df)

selected_df <- read.csv(selected_file, stringsAsFactors = FALSE)
cv_tbl <- read.csv(cv_lambda_path, stringsAsFactors = FALSE)

cv_tbl2 <- cv_tbl %>%
  mutate(
    best_lambda_zeta = ifelse(is.na(best_lambda_zeta) | best_lambda_zeta == 0,
                              lambda_zeta_default, best_lambda_zeta),
    best_lambda_eta  = ifelse(is.na(best_lambda_eta) | best_lambda_eta == 0,
                              lambda_eta_default, best_lambda_eta)
  )

lambda_map <- setNames(
  lapply(seq_len(nrow(cv_tbl2)), function(i) {
    list(
      lambda_zeta = as.numeric(cv_tbl2$best_lambda_zeta[i]),
      lambda_eta  = as.numeric(cv_tbl2$best_lambda_eta[i])
    )
  }),
  cv_tbl2$target_cancer
)

target_sheets <- cv_tbl2$target_cancer
target_sheets <- target_sheets[target_sheets %in% all_sheets]

if (!is.null(run_targets)) {
  target_sheets <- intersect(target_sheets, run_targets)
}

size_map_df <- get_size_map_cumh(proc, selected_df)

# =========================================================
# LOOP
# =========================================================
all_curve_rows <- list()
all_scb_rows <- list()
plot_list <- list()

for (tname in target_sheets) {
  message("[CUMH] ", tname)
  
  c_mult <- get_c_mult(tname, c_mult_default)
  target_info <- proc[[tname]]
  target_df0 <- target_info$df
  
  sel_row <- selected_df[selected_df$target_cancer == tname, , drop = FALSE]
  sel_sources <- character(0)
  if (nrow(sel_row) > 0 &&
      !is.na(sel_row$n_selected[1]) &&
      sel_row$n_selected[1] > 0) {
    sel_sources <- parse_sources(sel_row$selected_sources[1])
    sel_sources <- sel_sources[sel_sources %in% names(proc) & sel_sources != tname]
  }
  
  source_df0 <- pool_sources_df(sel_sources, proc)
  
  apply_list0 <- list()
  if (!is.null(source_df0)) apply_list0$src <- source_df0
  
  pre0 <- preprocess_features(
    target_df0,
    apply_list0,
    target_info$feature_cols,
    target_info$time_col,
    target_info$event_col
  )
  
  target_df <- pre0$train_df
  source_df <- if (!is.null(source_df0)) pre0$apply_df_list$src else NULL
  
  lam <- lambda_map[[tname]]
  if (is.null(lam)) {
    lam <- list(lambda_zeta = lambda_zeta_default, lambda_eta = lambda_eta_default)
  }
  
  sievetl_res <- tryCatch(
    compute_sievetl_cumh_real(
      target_df = target_df,
      source_df = source_df,
      feature_cols = target_info$feature_cols,
      time_col = target_info$time_col,
      event_col = target_info$event_col,
      lambda_zeta = lam$lambda_zeta,
      lambda_eta = lam$lambda_eta,
      c_mult = c_mult,
      tname = tname,
      pilot_dir = pilot_dir
    ),
    error = function(e) {
      message("  Cox-SieveTL error: ", conditionMessage(e))
      NULL
    }
  )
  
  if (is.null(sievetl_res)) {
    message("  skipped: Cox-SieveTL failed")
    next
  }
  
  coxph_res <- compute_coxph_target_cumh_real(
    target_df = target_df,
    feature_cols = target_info$feature_cols,
    time_col = target_info$time_col,
    event_col = target_info$event_col,
    quant_times = sievetl_res$curve_wald$time
  )
  
  if (is.null(coxph_res)) {
    message("  skipped: CoxPH failed")
    next
  }
  
  df_plot <- bind_rows(
    sievetl_res$curve_wald,
    coxph_res$curve
  ) %>%
    mutate(
      target_cancer = tname,
      cancer_label = ifelse(tname %in% names(title_map), unname(title_map[tname]), tname),
      c_mult = c_mult
    )
  
  df_scb <- sievetl_res$scb_curve %>%
    mutate(
      target_cancer = tname,
      cancer_label = ifelse(tname %in% names(title_map), unname(title_map[tname]), tname),
      c_mult = c_mult
    )
  
  size_row <- size_map_df %>% filter(target_cancer == tname)
  n_target_use <- if (nrow(size_row) > 0) size_row$n_target[1] else NA
  n_source_use <- if (nrow(size_row) > 0) size_row$n_source[1] else NA
  
  x_max <- if (tname == "cancer15") 40 else sievetl_res$combined_time_max
  y_max <- if (tname == "cancer15") 6 else NULL
  
  p <- plot_sievetl_vs_coxph_cumh(
    df_plot = df_plot,
    df_scb = df_scb,
    title_text = unique(df_plot$cancer_label),
    n_target = n_target_use,
    n_source = n_source_use,
    x_max = x_max,
    y_max = y_max
  )
  
  plot_list[[tname]] <- p
  all_curve_rows[[length(all_curve_rows) + 1]] <- df_plot
  all_scb_rows[[length(all_scb_rows) + 1]] <- df_scb
  
  ggsave(
    file.path(plot_dir, paste0("cumh_sievetl_vs_coxph_", tname, ".pdf")),
    p,
    width = 7,
    height = 5.5
  )
}

if (length(all_curve_rows) > 0) {
  write.csv(
    bind_rows(all_curve_rows),
    file.path(results_dir, "cumh_sievetl_vs_coxph_all_cancers.csv"),
    row.names = FALSE
  )
}

if (length(all_scb_rows) > 0) {
  write.csv(
    bind_rows(all_scb_rows),
    file.path(results_dir, "cumh_sievetl_scb_all_cancers.csv"),
    row.names = FALSE
  )
}

if (length(plot_list) > 0) {
  panel_plot <- wrap_plots(plot_list, ncol = 3, guides = "collect") &
    theme(legend.position = "bottom")
  
  panel_plot <- panel_plot +
    plot_annotation(title = "Baseline cumulative hazard: Cox-SieveTL vs CoxPH (Target)")
  
  ggsave(
    file.path(plot_dir, "cumh_sievetl_vs_coxph_panel.pdf"),
    panel_plot,
    width = 16,
    height = 20
  )
}

