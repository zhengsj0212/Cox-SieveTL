suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(survival)
  library(ggplot2)
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
out_dir   <- file.path(results_dir, "realdata_loris_survival_inference_dense")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# Controls
# =========================================================
lambda_zeta_default <- 0.005
lambda_eta_default  <- 0

quant_probs <- c(0.1, 0.9)
#seq(0.1, 0.9, by = 0.1)
integration_grid_size <- 2000
multiplier_boot <- 1000
conf_level <- 0.95
alpha <- 1 - conf_level
z_alpha <- qnorm(1 - alpha / 2)

# Fixed LORIS values, not estimated from data
loris_values_fixed <- c(0.1, 0.9)
#seq(0.1, 0.9, 0.1)
names(loris_values_fixed) <- paste0("loris_", formatC(loris_values_fixed, format = "f", digits = 1))
loris_col_override <- NULL


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
# Target list / lambda map
# =========================================================
cv_tbl <- read.csv("./results_allTMB_new_hessian_lasso/cv_metrics_summary_table_with_c_lambda_bic.csv", stringsAsFactors = FALSE)
c_map <- cv_tbl$c_mult
all_sheets <- paste0("cancer", c(1:10, 12:16))
names(c_map) <- as.character(all_sheets)
c_default <- 1
cv_tbl2 <- cv_tbl %>%
  mutate(
    best_lambda_zeta = ifelse(is.na(best_lambda_zeta) | best_lambda_zeta == 0,
                              lambda_zeta_default, best_lambda_zeta),
    best_lambda_eta  = ifelse(is.na(best_lambda_eta)  | best_lambda_eta  == 0,
                              lambda_eta_default,  best_lambda_eta)
  )

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

# =========================================================
# Helpers: inverse-Hessian inference
# =========================================================
.smooth_abs <- function(u, eps) sqrt(u * u + eps * eps) - eps
.smooth_abs_grad <- function(u, eps) u / sqrt(u * u + eps * eps)
choose_low_medium_high_profiles <- function(df, value_col = "loris_value") {
  vals <- sort(unique(df[[value_col]]))
  vals <- vals[is.finite(vals)]
  
  if (length(vals) < 3) {
    stop("Need at least 3 unique LORIS values to choose low / medium / high.")
  }
  
  low_val <- min(vals)
  high_val <- max(vals)
  med_val <- vals[which.min(abs(vals - median(vals)))]
  
  keep_vals <- c(low_val, med_val, high_val)
  
  df2 <- df[df[[value_col]] %in% keep_vals, , drop = FALSE]
  
  df2$profile3 <- ifelse(df2[[value_col]] == low_val, "Low",
                         ifelse(df2[[value_col]] == med_val, "Medium", "High"))
  
  df2$profile3 <- factor(df2$profile3, levels = c("Low", "Medium", "High"))
  df2
}

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

fit_inverse_hessian_correction <- function(H_target, omega_A, j, lambda_l1, eps = 1e-4) {
  d <- ncol(H_target)
  e_j <- rep(0, d)
  e_j[j] <- 1
  
  if (lambda_l1 == 0) {
    return(as.vector(solve(H_target, e_j) - omega_A))
  }
  
  obj <- function(delta) {
    omega <- omega_A + delta
    q <- 0.5 * drop(crossprod(omega, H_target %*% omega)) - sum(e_j * omega)
    q + lambda_l1 * sum(.smooth_abs(delta, eps))
  }
  
  grad <- function(delta) {
    omega <- omega_A + delta
    as.vector(H_target %*% omega - e_j + lambda_l1 * .smooth_abs_grad(delta, eps))
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

# =========================================================
# Baseline cumulative hazard helpers
# =========================================================
make_baseline_from_gamma <- function(gamma, g_funcs, t_grid) {
  t_grid <- sort(unique(as.numeric(t_grid)))
  if (length(t_grid) < 2L) stop("Need at least two grid points.")
  
  G <- sapply(g_funcs, function(gf) as.numeric(gf(t_grid)))
  if (is.null(dim(G))) G <- matrix(G, ncol = 1)
  
  logh0 <- drop(G %*% gamma)
  h0 <- exp(logh0)
  
  dt <- diff(t_grid)
  H0 <- c(0, cumsum(dt * (head(h0, -1) + tail(h0, -1)) / 2))
  
  data.frame(
    time = t_grid,
    logh0 = logh0,
    h0 = h0,
    H0 = H0
  )
}

compute_A_matrix <- function(gamma, g_funcs, eval_times) {
  base_curve <- make_baseline_from_gamma(gamma, g_funcs, eval_times)
  p <- length(gamma)
  
  G <- sapply(g_funcs, function(gf) as.numeric(gf(base_curve$time)))
  if (is.null(dim(G))) G <- matrix(G, ncol = 1)
  
  h0 <- as.numeric(base_curve$h0)
  dt <- diff(base_curve$time)
  
  A <- matrix(0, nrow = nrow(base_curve), ncol = p)
  for (j in seq_len(p)) {
    integrand <- G[, j] * h0
    A[, j] <- c(0, cumsum(dt * (head(integrand, -1) + tail(integrand, -1)) / 2))
  }
  
  list(
    time = base_curve$time,
    A = A,
    base_curve = base_curve
  )
}

predict_baseline_cumhaz_inference <- function(gamma_os, IF_gamma_mat, g_funcs, eval_times, n0) {
  Aobj <- compute_A_matrix(gamma_os, g_funcs, eval_times)
  
  IF_Lambda_mat <- IF_gamma_mat %*% t(Aobj$A)
  
  var_Lambda <- colSums(IF_Lambda_mat^2) / n0
  se_Lambda  <- sqrt(pmax(var_Lambda / n0, 0))
  
  list(
    time = Aobj$time,
    estimate = Aobj$base_curve$H0,
    h0 = Aobj$base_curve$h0,
    IF_mat = IF_Lambda_mat,
    var = var_Lambda,
    se = se_Lambda
  )
}

# =========================================================
# Main plugin inference
# =========================================================
compute_plugin_inference <- function(target, source, beta_hat, gamma_hat, g_funcs,
                                     lambda_delta_vec,
                                     quant_times = NULL,
                                     grid_range = NULL,
                                     quant_probs = seq(0.1, 0.9, by = 0.1),
                                     integration_grid_size = 2000) {
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
  
  dense_grid <- seq(grid_start, grid_end, length.out = max(as.integer(integration_grid_size), 2L))
  time_grid_target <- seq(min(stime, na.rm = TRUE), max(stime, na.rm = TRUE),
                          length.out = max(as.integer(integration_grid_size), 2L))
  
  pooled_stime <- c(target$stime, unlist(lapply(source, `[[`, "stime")))
  pooled_time_grid <- seq(min(pooled_stime, na.rm = TRUE), max(pooled_stime, na.rm = TRUE),
                          length.out = max(as.integer(integration_grid_size), 2L))
  
  H_target <- compute_target_hessian_theta(target, gamma_hat, beta_hat, g_funcs, time_grid_target)
  H_pooled <- compute_pooled_hessian_theta(target, source, gamma_hat, beta_hat, g_funcs, pooled_time_grid)
  Omega_hat_by_column <- estimate_transfer_inverse_hessian_by_column(H_pooled, H_target, lambda_delta_vec)
  
  theta_hat <- c(gamma_hat, beta_hat)
  score_vec <- compute_target_score_theta(target, gamma_hat, beta_hat, g_funcs, time_grid_target)
  score_subject_mat <- compute_target_score_theta_subjects(target, gamma_hat, beta_hat, g_funcs, time_grid_target)
  
  IF_theta_hat_mat <- -score_subject_mat %*% t(Omega_hat_by_column)
  IF_gamma_hat_mat <- IF_theta_hat_mat[, seq_len(p), drop = FALSE]
  IF_beta_hat_mat  <- IF_theta_hat_mat[, (p + 1):(p + L), drop = FALSE]
  
  theta_db_by_column <- as.vector(theta_hat - Omega_hat_by_column %*% score_vec)
  gamma_os <- theta_db_by_column[seq_len(p)]
  beta_os  <- theta_db_by_column[(p + 1):(p + L)]
  
  if (is.null(quant_times)) {
    quant_times <- as.numeric(stats::quantile(stime[is.finite(stime)],
                                              probs = quant_probs, type = 7, names = FALSE))
  } else {
    quant_times <- as.numeric(quant_times)
  }
  
  base_quant <- predict_baseline_cumhaz_inference(
    gamma_os = gamma_os,
    IF_gamma_mat = IF_gamma_hat_mat,
    g_funcs = g_funcs,
    eval_times = quant_times,
    n0 = n0
  )
  
  base_dense <- predict_baseline_cumhaz_inference(
    gamma_os = gamma_os,
    IF_gamma_mat = IF_gamma_hat_mat,
    g_funcs = g_funcs,
    eval_times = dense_grid,
    n0 = n0
  )
  
  list(
    n0 = n0,
    gamma_os = gamma_os,
    beta_hat_os = beta_os,
    IF_gamma_hat_mat = IF_gamma_hat_mat,
    IF_beta_hat_mat = IF_beta_hat_mat,
    quant_times = quant_times,
    Lambda_hat_os_quant = base_quant$estimate,
    se_Lambda_quant = base_quant$se,
    IF_Lambda_hat_mat = base_quant$IF_mat,
    integration_grid = dense_grid,
    Lambda_hat_os_grid = base_dense$estimate,
    h0_hat_os_grid = base_dense$h0,
    IF_Lambda_hat_grid = base_dense$IF_mat
  )
}

# =========================================================
# Multiplier CI/SCB for baseline
# =========================================================
compute_multiplier_cumh_intervals <- function(inference,
                                              alpha = 0.05,
                                              multiplier_boot = 1000,
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
# LORIS profile helpers
# =========================================================
detect_loris_col <- function(feature_cols, override = NULL) {
  if (!is.null(override)) {
    if (!(override %in% feature_cols)) {
      stop("Specified loris_col_override not found in feature_cols: ", override)
    }
    return(override)
  }
  
  hits <- feature_cols[grepl("loris", feature_cols, ignore.case = TRUE)]
  if (length(hits) == 0L) {
    stop("Could not detect a LORIS column in feature_cols.")
  }
  if (length(hits) > 1L) {
    stop("Multiple LORIS-like columns found: ", paste(hits, collapse = ", "),
         ". Please set loris_col_override.")
  }
  hits
}

build_reference_profile <- function(target_df, feature_cols, loris_col, loris_value) {
  
  get_mode <- function(x) {
    ux <- unique(x[!is.na(x)])
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  x_ref <- sapply(feature_cols, function(v) {
    xv <- target_df[[v]]
    
    if (is.numeric(xv)) {
      # detect binary (categorical encoded as 0/1)
      ux <- unique(xv[!is.na(xv)])
      
      if (length(ux) <= 2 && all(ux %in% c(0, 1))) {
        # binary → mode
        get_mode(xv)
      } else {
        # continuous → mean
        mean(xv, na.rm = TRUE)
      }
      
    } else {
      # fallback (shouldn't happen after preprocess)
      get_mode(xv)
    }
  })
  
  x_ref <- as.numeric(x_ref)
  names(x_ref) <- feature_cols
  
  # override LORIS
  x_ref[loris_col] <- loris_value
  
  return(x_ref)
}
# =========================================================
# Survival prediction + dense inference
# =========================================================
compute_profile_survival_inference <- function(inference, x_ref,
                                               alpha = 0.05,
                                               multiplier_boot = 1000,
                                               multiplier = c("normal", "rademacher"),
                                               seed = NULL,
                                               sd_floor = 1e-8) {
  multiplier <- match.arg(multiplier)
  n0 <- inference$n0
  
  x_ref <- as.numeric(x_ref)
  stopifnot(length(x_ref) == ncol(inference$IF_beta_hat_mat))
  
  lp_hat <- sum(x_ref * inference$beta_hat_os)
  risk_hat <- exp(lp_hat)
  
  S_quant_hat <- exp(-risk_hat * inference$Lambda_hat_os_quant)
  S_grid_hat  <- exp(-risk_hat * inference$Lambda_hat_os_grid)
  
  IF_lp <- as.numeric(inference$IF_beta_hat_mat %*% x_ref)
  
  IF_S_quant <- - outer(rep(1, n0), S_quant_hat * risk_hat) *
    (inference$IF_Lambda_hat_mat + outer(IF_lp, inference$Lambda_hat_os_quant))
  
  IF_S_grid <- - outer(rep(1, n0), S_grid_hat * risk_hat) *
    (inference$IF_Lambda_hat_grid + outer(IF_lp, inference$Lambda_hat_os_grid))
  
  var_quant <- colSums(IF_S_quant^2) / n0
  se_quant  <- sqrt(pmax(var_quant / n0, 0))
  
  var_grid <- colSums(IF_S_grid^2) / n0
  se_grid  <- sqrt(pmax(var_grid / n0, 0))
  
  z_alpha <- stats::qnorm(1 - alpha / 2)
  
  ci_lower_quant_wald <- pmax(S_quant_hat - z_alpha * se_quant, 0)
  ci_upper_quant_wald <- pmin(S_quant_hat + z_alpha * se_quant, 1)
  
  ci_lower_grid_wald <- pmax(S_grid_hat - z_alpha * se_grid, 0)
  ci_upper_grid_wald <- pmin(S_grid_hat + z_alpha * se_grid, 1)
  
  if (!is.null(seed)) set.seed(seed)
  if (multiplier == "normal") {
    xi <- matrix(rnorm(multiplier_boot * n0), nrow = multiplier_boot, ncol = n0)
  } else {
    xi <- matrix(sample(c(-1, 1), size = multiplier_boot * n0, replace = TRUE),
                 nrow = multiplier_boot, ncol = n0)
  }
  
  z_quant <- (xi %*% IF_S_quant) / sqrt(n0)
  s_quant_boot <- sqrt(pmax(colMeans(z_quant^2), 0))
  se_quant_boot <- s_quant_boot / sqrt(n0)
  ci_lower_quant_boot <- pmax(S_quant_hat - z_alpha * se_quant_boot, 0)
  ci_upper_quant_boot <- pmin(S_quant_hat + z_alpha * se_quant_boot, 1)
  
  z_grid <- (xi %*% IF_S_grid) / sqrt(n0)
  s_grid_boot <- sqrt(pmax(colMeans(z_grid^2), 0))
  se_grid_boot <- s_grid_boot / sqrt(n0)
  
  ci_lower_grid_boot <- pmax(S_grid_hat - z_alpha * se_grid_boot, 0)
  ci_upper_grid_boot <- pmin(S_grid_hat + z_alpha * se_grid_boot, 1)
  
  s_grid_safe <- pmax(s_grid_boot, sd_floor)
  studentized <- abs(sweep(z_grid, 2, s_grid_safe, "/"))
  sup_stats <- apply(studentized, 1, max, na.rm = TRUE)
  scb_crit <- as.numeric(stats::quantile(sup_stats, probs = 1 - alpha, names = FALSE, na.rm = TRUE))
  scb_halfwidth <- scb_crit * s_grid_safe / sqrt(n0)
  
  scb_lower_grid <- pmax(S_grid_hat - scb_halfwidth, 0)
  scb_upper_grid <- pmin(S_grid_hat + scb_halfwidth, 1)
  
  list(
    x_ref = x_ref,
    lp_hat = lp_hat,
    risk_hat = risk_hat,
    S_quant_hat = S_quant_hat,
    S_grid_hat = S_grid_hat,
    IF_S_quant = IF_S_quant,
    IF_S_grid = IF_S_grid,
    se_quant = se_quant,
    se_grid = se_grid,
    pointwise = list(
      lower_wald = ci_lower_quant_wald,
      upper_wald = ci_upper_quant_wald,
      lower_boot = ci_lower_quant_boot,
      upper_boot = ci_upper_quant_boot,
      se_wald = se_quant,
      se_boot = se_quant_boot
    ),
    dense_pointwise = list(
      lower_wald = ci_lower_grid_wald,
      upper_wald = ci_upper_grid_wald,
      lower_boot = ci_lower_grid_boot,
      upper_boot = ci_upper_grid_boot,
      se_wald = se_grid,
      se_boot = se_grid_boot
    ),
    scb = list(
      crit = scb_crit,
      lower = scb_lower_grid,
      upper = scb_upper_grid,
      s_grid = s_grid_boot
    )
  )
}

# =========================================================
# Fit SieveTL on one cancer + transfer inference
# =========================================================
run_fit_and_inference_one <- function(tname) {
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
    alpha = alpha,
    multiplier_boot = multiplier_boot,
    multiplier = "normal",
    seed = 2026
  )
  
  list(
    tname = tname,
    target = target,
    target_df_raw = target_df_raw,
    target_df = target_df,
    fit = fit,
    inference = inference,
    baseline_bootstrap = boot_out
  )
}

# =========================================================
# Build low / medium / high LORIS survival predictions
# =========================================================
run_loris_survival_one <- function(tname,
                                   loris_col_override = NULL,
                                   loris_values_fixed = c(low = 0.25, medium = 0.50, high = 0.75)) {
  obj <- run_fit_and_inference_one(tname)
  
  target <- obj$target
  target_df <- obj$target_df
  inference <- obj$inference
  baseline_boot <- obj$baseline_bootstrap
  
  loris_col <- detect_loris_col(target$feature_cols, override = loris_col_override)
  loris_values <- loris_values_fixed
  
  baseline_pointwise <- data.frame(
    target_cancer = tname,
    target_label = ifelse(tname %in% names(pretty_name_map), pretty_name_map[[tname]], tname),
    time = inference$quant_times,
    Lambda0_estimate = inference$Lambda_hat_os_quant,
    Lambda0_se = inference$se_Lambda_quant,
    Lambda0_ci_lower_wald = inference$Lambda_hat_os_quant - z_alpha * inference$se_Lambda_quant,
    Lambda0_ci_upper_wald = inference$Lambda_hat_os_quant + z_alpha * inference$se_Lambda_quant,
    Lambda0_ci_lower_boot = baseline_boot$pointwise$lower,
    Lambda0_ci_upper_boot = baseline_boot$pointwise$upper,
    stringsAsFactors = FALSE
  )
  
  baseline_dense <- data.frame(
    target_cancer = tname,
    target_label = ifelse(tname %in% names(pretty_name_map), pretty_name_map[[tname]], tname),
    time = inference$integration_grid,
    Lambda0_estimate = inference$Lambda_hat_os_grid,
    h0_estimate = inference$h0_hat_os_grid,
    Lambda0_scb_lower = baseline_boot$scb$lower,
    Lambda0_scb_upper = baseline_boot$scb$upper,
    Lambda0_scb_crit = baseline_boot$scb$crit,
    stringsAsFactors = FALSE
  )
  
  profile_pointwise_list <- list()
  profile_dense_pointwise_list <- list()
  profile_scb_list <- list()
  profile_xref_list <- list()
  
  for (lvl in names(loris_values)) {
    x_ref <- build_reference_profile(
      target_df = target_df,
      feature_cols = target$feature_cols,
      loris_col = loris_col,
      loris_value = loris_values[[lvl]]
    )
    
    prof <- compute_profile_survival_inference(
      inference = inference,
      x_ref = x_ref,
      alpha = alpha,
      multiplier_boot = multiplier_boot,
      multiplier = "normal",
      seed = 2026
    )
    
    profile_xref_list[[lvl]] <- data.frame(
      target_cancer = tname,
      target_label = ifelse(tname %in% names(pretty_name_map), pretty_name_map[[tname]], tname),
      profile = lvl,
      loris_col = loris_col,
      loris_value = loris_values[[lvl]],
      covariate = names(x_ref),
      value = as.numeric(x_ref),
      stringsAsFactors = FALSE
    )
    
    profile_pointwise_list[[lvl]] <- data.frame(
      target_cancer = tname,
      target_label = ifelse(tname %in% names(pretty_name_map), pretty_name_map[[tname]], tname),
      curve_type = "predicted_survival",
      profile = lvl,
      loris_col = loris_col,
      loris_value = loris_values[[lvl]],
      time = inference$quant_times,
      estimate = prof$S_quant_hat,
      se_wald = prof$pointwise$se_wald,
      se_boot = prof$pointwise$se_boot,
      ci_lower_wald = prof$pointwise$lower_wald,
      ci_upper_wald = prof$pointwise$upper_wald,
      ci_lower_boot = prof$pointwise$lower_boot,
      ci_upper_boot = prof$pointwise$upper_boot,
      lp_hat = prof$lp_hat,
      risk_hat = prof$risk_hat,
      stringsAsFactors = FALSE
    )
    
    profile_dense_pointwise_list[[lvl]] <- data.frame(
      target_cancer = tname,
      target_label = ifelse(tname %in% names(pretty_name_map), pretty_name_map[[tname]], tname),
      curve_type = "predicted_survival",
      profile = lvl,
      loris_col = loris_col,
      loris_value = loris_values[[lvl]],
      time = inference$integration_grid,
      estimate = prof$S_grid_hat,
      se_wald = prof$dense_pointwise$se_wald,
      se_boot = prof$dense_pointwise$se_boot,
      ci_lower_wald = prof$dense_pointwise$lower_wald,
      ci_upper_wald = prof$dense_pointwise$upper_wald,
      ci_lower_boot = prof$dense_pointwise$lower_boot,
      ci_upper_boot = prof$dense_pointwise$upper_boot,
      lp_hat = prof$lp_hat,
      risk_hat = prof$risk_hat,
      stringsAsFactors = FALSE
    )
    
    profile_scb_list[[lvl]] <- data.frame(
      target_cancer = tname,
      target_label = ifelse(tname %in% names(pretty_name_map), pretty_name_map[[tname]], tname),
      curve_type = "predicted_survival",
      profile = lvl,
      loris_col = loris_col,
      loris_value = loris_values[[lvl]],
      time = inference$integration_grid,
      estimate = prof$S_grid_hat,
      scb_lower = prof$scb$lower,
      scb_upper = prof$scb$upper,
      scb_crit = prof$scb$crit,
      lp_hat = prof$lp_hat,
      risk_hat = prof$risk_hat,
      stringsAsFactors = FALSE
    )
  }
  
  profile_pointwise <- bind_rows(profile_pointwise_list)
  profile_dense_pointwise <- bind_rows(profile_dense_pointwise_list)
  profile_scb <- bind_rows(profile_scb_list)
  profile_xref <- bind_rows(profile_xref_list)
  
  write.csv(
    baseline_pointwise,
    file.path(out_dir, paste0("baseline_pointwise_", tname, ".csv")),
    row.names = FALSE
  )
  write.csv(
    baseline_dense,
    file.path(out_dir, paste0("baseline_dense_", tname, ".csv")),
    row.names = FALSE
  )
  write.csv(
    profile_pointwise,
    file.path(out_dir, paste0("loris_survival_pointwise_", tname, ".csv")),
    row.names = FALSE
  )
  write.csv(
    profile_dense_pointwise,
    file.path(out_dir, paste0("loris_survival_dense_pointwise_", tname, ".csv")),
    row.names = FALSE
  )
  write.csv(
    profile_scb,
    file.path(out_dir, paste0("loris_survival_scb_", tname, ".csv")),
    row.names = FALSE
  )
  write.csv(
    profile_xref,
    file.path(out_dir, paste0("loris_survival_reference_values_", tname, ".csv")),
    row.names = FALSE
  )
  
  list(
    baseline_pointwise = baseline_pointwise,
    baseline_dense = baseline_dense,
    profile_pointwise = profile_pointwise,
    profile_dense_pointwise = profile_dense_pointwise,
    profile_scb = profile_scb,
    profile_xref = profile_xref,
    loris_col = loris_col,
    loris_values = loris_values
  )
}

# =========================================================
# Run all targets
# =========================================================
run_all_loris_survival <- function(target_vec = target_sheets,
                                   loris_col_override = NULL,
                                   loris_values_fixed = c(low = 0.25, medium = 0.50, high = 0.75)) {
  baseline_pw_list <- list()
  baseline_dense_list <- list()
  profile_pw_list <- list()
  profile_dense_pw_list <- list()
  profile_scb_list <- list()
  xref_list <- list()
  
  for (tt in target_vec) {
    message(sprintf("[LORIS-SURV] running target = %s", tt))
    out <- tryCatch(
      run_loris_survival_one(
        tname = tt,
        loris_col_override = loris_col_override,
        loris_values_fixed = loris_values_fixed
      ),
      error = function(e) {
        message(sprintf("[LORIS-SURV] target=%s failed: %s", tt, conditionMessage(e)))
        NULL
      }
    )
    
    if (!is.null(out)) {
      baseline_pw_list[[length(baseline_pw_list) + 1]] <- out$baseline_pointwise
      baseline_dense_list[[length(baseline_dense_list) + 1]] <- out$baseline_dense
      profile_pw_list[[length(profile_pw_list) + 1]] <- out$profile_pointwise
      profile_dense_pw_list[[length(profile_dense_pw_list) + 1]] <- out$profile_dense_pointwise
      profile_scb_list[[length(profile_scb_list) + 1]] <- out$profile_scb
      xref_list[[length(xref_list) + 1]] <- out$profile_xref
    }
  }
  
  baseline_pw_all <- bind_rows(baseline_pw_list)
  baseline_dense_all <- bind_rows(baseline_dense_list)
  profile_pw_all <- bind_rows(profile_pw_list)
  profile_dense_pw_all <- bind_rows(profile_dense_pw_list)
  profile_scb_all <- bind_rows(profile_scb_list)
  xref_all <- bind_rows(xref_list)
  
  write.csv(
    baseline_pw_all,
    file.path(out_dir, "all_targets_baseline_pointwise.csv"),
    row.names = FALSE
  )
  write.csv(
    baseline_dense_all,
    file.path(out_dir, "all_targets_baseline_dense.csv"),
    row.names = FALSE
  )
  write.csv(
    profile_pw_all,
    file.path(out_dir, "all_targets_loris_survival_pointwise.csv"),
    row.names = FALSE
  )
  write.csv(
    profile_dense_pw_all,
    file.path(out_dir, "all_targets_loris_survival_dense_pointwise.csv"),
    row.names = FALSE
  )
  write.csv(
    profile_scb_all,
    file.path(out_dir, "all_targets_loris_survival_scb.csv"),
    row.names = FALSE
  )
  write.csv(
    xref_all,
    file.path(out_dir, "all_targets_loris_survival_reference_values.csv"),
    row.names = FALSE
  )
  
  list(
    baseline_pointwise = baseline_pw_all,
    baseline_dense = baseline_dense_all,
    profile_pointwise = profile_pw_all,
    profile_dense_pointwise = profile_dense_pw_all,
    profile_scb = profile_scb_all,
    profile_xref = xref_all
  )
}


# =========================================================
# Inputs
# =========================================================
tname <- "cancer12"   # ovarian
ref_file <- file.path(out_dir, "all_targets_loris_survival_reference_values.csv")
results_dir   <- "./results_allTMB_new_hessian_lasso"
# =========================================================
# Read and preprocess ovarian data the same way
# =========================================================
df_raw <- read_excel(data_file, sheet = tname)
df_raw <- cap_cancer_df(df_raw)
df_raw <- na.omit(df_raw)

target_obj <- prep_df(df_raw)

pre <- preprocess_features(
  target_obj$df,
  list(),   # same preprocessing structure
  feature_cols = target_obj$feature_cols,
  time_col = target_obj$time_col,
  event_col = target_obj$event_col
)

target_df <- pre$train_df
feature_cols <- target_obj$feature_cols
time_col <- target_obj$time_col
event_col <- target_obj$event_col

# =========================================================
# Fit CoxPH on the same processed target data
# =========================================================
cox_formula <- as.formula(
  paste0(
    "survival::Surv(", time_col, ", ", event_col, ") ~ ",
    paste(feature_cols, collapse = " + ")
  )
)

fit_cox <- survival::coxph(
  cox_formula,
  data = target_df,
  x = TRUE,
  y = TRUE,
  model = TRUE
)

# =========================================================
# Read the SAME reference covariate values used by Cox-SieveTL
# =========================================================
ref_long <- read.csv(ref_file, stringsAsFactors = FALSE)

ref_long_sub <- ref_long %>%
  filter(target_cancer == tname)

# optional: check available profiles
unique(ref_long_sub$profile)

# =========================================================
# Convert one profile at a time to newdata
# =========================================================
make_newdata_from_reference <- function(ref_df_long, profile_name, feature_cols) {
  x_ref <- ref_df_long %>%
    filter(profile == profile_name) %>%
    dplyr::select(covariate, value)
  
  if (nrow(x_ref) == 0) {
    stop("No reference values found for profile: ", profile_name)
  }
  
  newdat <- x_ref %>%
    tidyr::pivot_wider(names_from = covariate, values_from = value)
  
  # enforce same column order as model features
  missing_cols <- setdiff(feature_cols, names(newdat))
  if (length(missing_cols) > 0) {
    stop("Missing covariates in reference profile: ", paste(missing_cols, collapse = ", "))
  }
  
  newdat <- newdat[, feature_cols, drop = FALSE]
  as.data.frame(newdat)
}

profiles_use <- unique(ref_long_sub$profile)

coxph_pred_list <- lapply(profiles_use, function(pf) {
  newdat <- make_newdata_from_reference(
    ref_df_long = ref_long_sub,
    profile_name = pf,
    feature_cols = feature_cols
  )
  
  sf <- survfit(
    fit_cox,
    newdata = newdat,
    conf.int = 0.95,
    conf.type = "log-log"
  )
  
  data.frame(
    target_cancer = tname,
    target_label = ifelse(tname %in% names(pretty_name_map), pretty_name_map[[tname]], tname),
    method = "coxph",
    profile = pf,
    time = sf$time,
    estimate = sf$surv,
    ci_lower = sf$lower,
    ci_upper = sf$upper,
    stringsAsFactors = FALSE
  )
})

coxph_pred_df <- bind_rows(coxph_pred_list)

write.csv(
  coxph_pred_df,
  file.path(out_dir, paste0("coxph_survival_same_reference_", tname, ".csv")),
  row.names = FALSE
)

# =========================================================
# Plot
# =========================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

tname <- "cancer12"

f_sievetl <- file.path(out_dir, paste0("loris_survival_dense_pointwise_", tname, ".csv"))
f_coxph   <- file.path(out_dir, paste0("coxph_survival_same_reference_", tname, ".csv"))

df_sievetl <- read.csv(f_sievetl, stringsAsFactors = FALSE) %>%
  mutate(
    method = "Cox-SieveTL",
    ci_lower = ci_lower_wald,
    ci_upper = ci_upper_wald
  ) %>%
  dplyr::select(target_cancer, target_label, method, profile, time, estimate, ci_lower, ci_upper) %>%
  filter(profile %in% c('loris_0.1', 'loris_0.9'))

df_coxph <- read.csv(f_coxph, stringsAsFactors = FALSE) %>%
  mutate(
    method = "CoxPH"
  ) %>%
  dplyr::select(target_cancer, target_label, method, profile, time, estimate, ci_lower, ci_upper) %>%
  filter(profile %in% c('loris_0.1', 'loris_0.9'))

# keep only common profiles if needed
common_profiles <- intersect(unique(df_sievetl$profile), unique(df_coxph$profile))

df_plot <- bind_rows(
  df_sievetl %>% filter(profile %in% common_profiles),
  df_coxph %>% filter(profile %in% common_profiles)
)

# optional relabel
# df_plot$profile <- factor(df_plot$profile, levels = sort(unique(df_plot$profile)))
# df_plot$method  <- factor(df_plot$method, levels = c("Cox-SieveTL", "CoxPH"))
# 
# df_plot$curve <- with(df_plot, paste(method, profile, sep = " | "))
# 
# df_plot$curve <- factor(
#   df_plot$curve,
#   levels = c(
#     "Cox-SieveTL | loris_0.1",
#     "Cox-SieveTL | loris_0.9",
#     "CoxPH | loris_0.1",
#     "CoxPH | loris_0.9"
#   ),
#   labels = c(
#     "SieveTL (L=0.1)",
#     "SieveTL (L=0.9)",
#     "CoxPH (L=0.1)",
#     "CoxPH (L=0.9)"
#   )
# )
# 
# p <- ggplot(
#   df_plot,
#   aes(x = time, y = estimate,
#       color = curve, fill = curve, group = curve)
# ) +
#   geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
#               alpha = 0.12, color = NA) +
#   geom_line(linewidth = 1.3) +
#   
#   coord_cartesian(xlim = c(0, 40), ylim = c(0, 1)) +
#   
#   # 🔥 High-contrast, colorblind-safe palette
#   scale_color_manual(values = c(
#     "SieveTL (L=0.1)" = "darkblue",
#     "SieveTL (L=0.9)" = "purple",
#     "CoxPH (L=0.1)"   = "darkred",
#     "CoxPH (L=0.9)"   = "blue"
#   )) +
#   
#   scale_fill_manual(values = c(
#     "SieveTL (L=0.1)" = "darkblue",
#     "SieveTL (L=0.9)" = "purple",
#     "CoxPH (L=0.1)"   = "darkred",
#     "CoxPH (L=0.9)"   = "blue"
#   )) +
#   
#   labs(
#     title = paste0(unique(df_plot$target_label),
#                    ": Cox-SieveTL vs CoxPH"),
#     x = "Time",
#     y = "Predicted survival"
#   ) +
#   
#   theme_bw() +
#   theme(
#     legend.position = "bottom",
#     legend.box = "vertical"
#   )
# 
# print(p)
df_plot$profile <- factor(
  df_plot$profile,
  levels = c("loris_0.1", "loris_0.9"),
  labels = c("LORIS = 0.1", "LORIS = 0.9")
)

df_plot$method <- factor(
  df_plot$method,
  levels = c("Cox-SieveTL", "CoxPH")
)

p <- ggplot(
  df_plot,
  aes(x = time, y = estimate, color = profile, fill = profile, group = profile)
) +
  geom_ribbon(
    aes(ymin = ci_lower, ymax = ci_upper),
    alpha = 0.18,
    color = NA
  ) +
  geom_line(linewidth = 1.3) +
  facet_wrap(~ method, nrow = 1) +
  coord_cartesian(xlim = c(0, 40), ylim = c(0, 1)) +
  scale_color_manual(values = c(
    "LORIS = 0.1" = "#1F77B4",
    "LORIS = 0.9" = "#D62728"
  )) +
  scale_fill_manual(values = c(
    "LORIS = 0.1" = "#9ECAE1",
    "LORIS = 0.9" = "#FC9272"
  )) +
  labs(
    title = paste0(unique(df_plot$target_label), ": Predicted Survival by Method"),
    x = "Time",
    y = "Predicted survival",
    color = NULL,
    fill = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 11)
  )

print(p)
ggsave(
  file.path(out_dir, "ovarian_loris_survival_2panel.pdf"),
  p,
  width = 8,
  height = 4,
  device = cairo_pdf   # better font rendering
)

############# All cancer

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(survival)
})

# =========================================================
# Inputs
# =========================================================
ref_file <- file.path(out_dir, "all_targets_loris_survival_reference_values.csv")

# =========================================================
# Helper: build one newdata row from reference table
# =========================================================
make_newdata_from_reference <- function(ref_df_long, profile_name, feature_cols) {
  x_ref <- ref_df_long %>%
    filter(profile == profile_name) %>%
    dplyr::select(covariate, value)
  
  if (nrow(x_ref) == 0) {
    stop("No reference values found for profile: ", profile_name)
  }
  
  newdat <- x_ref %>%
    tidyr::pivot_wider(names_from = covariate, values_from = value)
  
  missing_cols <- setdiff(feature_cols, names(newdat))
  if (length(missing_cols) > 0) {
    stop("Missing covariates in reference profile: ", paste(missing_cols, collapse = ", "))
  }
  
  newdat <- newdat[, feature_cols, drop = FALSE]
  as.data.frame(newdat)
}

# =========================================================
# Main: generate CoxPH survival predictions for all cancers
# =========================================================
generate_all_cancers_coxph_same_reference <- function(
    target_vec = target_sheets,
    data_file = data_file,
    ref_file = file.path(out_dir, "all_targets_loris_survival_reference_values.csv"),
    out_dir = out_dir
) {
  ref_long <- read.csv(ref_file, stringsAsFactors = FALSE)
  
  all_out <- list()
  
  for (tname in target_vec) {
    message(sprintf("[COXPH] running %s", tname))
    
    out <- tryCatch({
      
      # -----------------------------
      # Read and preprocess cancer data
      # -----------------------------
      df_raw <- read_excel(data_file, sheet = tname)
      df_raw <- cap_cancer_df(df_raw)
      df_raw <- na.omit(df_raw)
      
      target_obj <- prep_df(df_raw)
      
      pre <- preprocess_features(
        target_obj$df,
        list(),
        feature_cols = target_obj$feature_cols,
        time_col = target_obj$time_col,
        event_col = target_obj$event_col
      )
      
      target_df <- pre$train_df
      feature_cols <- target_obj$feature_cols
      time_col <- target_obj$time_col
      event_col <- target_obj$event_col
      
      # -----------------------------
      # Fit CoxPH
      # -----------------------------
      cox_formula <- as.formula(
        paste0(
          "survival::Surv(", time_col, ", ", event_col, ") ~ ",
          paste(feature_cols, collapse = " + ")
        )
      )
      
      fit_cox <- survival::coxph(
        cox_formula,
        data = target_df,
        x = TRUE,
        y = TRUE,
        model = TRUE
      )
      
      # -----------------------------
      # Reference profiles for this cancer
      # -----------------------------
      ref_long_sub <- ref_long %>%
        filter(target_cancer == tname)
      
      if (nrow(ref_long_sub) == 0) {
        stop("No reference values found in reference file for ", tname)
      }
      
      profiles_use <- unique(ref_long_sub$profile)
      
      # -----------------------------
      # Predict one profile at a time
      # -----------------------------
      coxph_pred_list <- lapply(profiles_use, function(pf) {
        newdat <- make_newdata_from_reference(
          ref_df_long = ref_long_sub,
          profile_name = pf,
          feature_cols = feature_cols
        )
        
        sf <- survfit(
          fit_cox,
          newdata = newdat,
          conf.int = 0.95,
          conf.type = "log-log"
        )
        
        data.frame(
          target_cancer = tname,
          target_label = ifelse(
            tname %in% names(pretty_name_map),
            pretty_name_map[[tname]],
            tname
          ),
          method = "coxph",
          profile = pf,
          time = sf$time,
          estimate = sf$surv,
          ci_lower = sf$lower,
          ci_upper = sf$upper,
          stringsAsFactors = FALSE
        )
      })
      
      coxph_pred_df <- bind_rows(coxph_pred_list)
      
      # -----------------------------
      # Save per-cancer file
      # -----------------------------
      out_file <- file.path(out_dir, paste0("coxph_survival_same_reference_", tname, ".csv"))
      write.csv(coxph_pred_df, out_file, row.names = FALSE)
      
      message(sprintf("[COXPH] saved %s", out_file))
      
      coxph_pred_df
      
    }, error = function(e) {
      message(sprintf("[COXPH] %s failed: %s", tname, conditionMessage(e)))
      NULL
    })
    
    if (!is.null(out)) {
      all_out[[length(all_out) + 1]] <- out
    }
  }
  
  all_df <- bind_rows(all_out)
  
  write.csv(
    all_df,
    file.path(out_dir, "all_cancers_coxph_survival_same_reference.csv"),
    row.names = FALSE
  )
  
  invisible(all_df)
}

all_coxph_pred <- generate_all_cancers_coxph_same_reference(
  target_vec = target_sheets,
  data_file = data_file,
  ref_file = file.path(out_dir, "all_targets_loris_survival_reference_values.csv"),
  out_dir = out_dir
)

# Plot

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
})

plot_all_cancers_2panel_pdf <- function(
    target_vec = target_sheets,
    out_dir,
    pretty_name_map,
    pdf_file = file.path(out_dir, "all_cancers_loris_survival_2panel.pdf"),
    profiles_keep = c("loris_0.1", "loris_0.9"),
    x_max = 40,
    plots_per_page = 3   # number of cancer-plots per page
) {
  plot_list <- list()
  
  for (tname in target_vec) {
    f_sievetl <- file.path(out_dir, paste0("loris_survival_dense_pointwise_", tname, ".csv"))
    f_coxph   <- file.path(out_dir, paste0("coxph_survival_same_reference_", tname, ".csv"))
    
    if (!file.exists(f_sievetl) || !file.exists(f_coxph)) {
      message(sprintf("[SKIP] Missing file(s) for %s", tname))
      next
    }
    
    df_sievetl <- read.csv(f_sievetl, stringsAsFactors = FALSE) %>%
      mutate(
        method = "Cox-SieveTL",
        ci_lower = ci_lower_wald,
        ci_upper = ci_upper_wald
      ) %>%
      dplyr::select(
        target_cancer, target_label, method, profile,
        time, estimate, ci_lower, ci_upper
      ) %>%
      filter(profile %in% profiles_keep)
    
    df_coxph <- read.csv(f_coxph, stringsAsFactors = FALSE) %>%
      mutate(method = "CoxPH") %>%
      dplyr::select(
        target_cancer, target_label, method, profile,
        time, estimate, ci_lower, ci_upper
      ) %>%
      filter(profile %in% profiles_keep)
    
    common_profiles <- intersect(unique(df_sievetl$profile), unique(df_coxph$profile))
    if (length(common_profiles) == 0) {
      message(sprintf("[SKIP] No common profiles for %s", tname))
      next
    }
    
    df_plot <- bind_rows(
      df_sievetl %>% filter(profile %in% common_profiles),
      df_coxph %>% filter(profile %in% common_profiles)
    )
    
    if (nrow(df_plot) == 0) {
      message(sprintf("[SKIP] Empty plot data for %s", tname))
      next
    }
    
    df_plot$profile <- factor(
      df_plot$profile,
      levels = profiles_keep,
      labels = paste0("LORIS = ", sub("loris_", "", profiles_keep))
    )
    
    df_plot$method <- factor(
      df_plot$method,
      levels = c("Cox-SieveTL", "CoxPH")
    )
    
    cancer_title <- if (tname %in% names(pretty_name_map)) pretty_name_map[[tname]] else tname
    x_max <- ifelse(tname == 'cancer13', 30, x_max)
    p <- ggplot(
      df_plot,
      aes(x = time, y = estimate, color = profile, fill = profile, group = profile)
    ) +
      geom_ribbon(
        aes(ymin = ci_lower, ymax = ci_upper),
        alpha = 0.18,
        color = NA
      ) +
      geom_line(linewidth = 1.0) +
      facet_wrap(~ method, nrow = 1) +
      coord_cartesian(xlim = c(0, x_max), ylim = c(0, 1)) +
      scale_color_manual(values = c(
        "LORIS = 0.1" = "#1F77B4",
        "LORIS = 0.9" = "#D62728"
      )) +
      scale_fill_manual(values = c(
        "LORIS = 0.1" = "#9ECAE1",
        "LORIS = 0.9" = "#FC9272"
      )) +
      labs(
        title = paste0(cancer_title, ": Predicted Survival by Method"),
        x = "Time",
        y = "Predicted survival",
        color = NULL,
        fill = NULL
      ) +
      theme_bw() +
      theme(
        legend.position = "bottom",
        strip.text = element_text(size = 10),
        plot.title = element_text(size = 12),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 9)
      )
    
    plot_list[[length(plot_list) + 1]] <- p
    message(sprintf("[DONE] %s", tname))
  }
  
  if (length(plot_list) == 0) {
    stop("No plots were created.")
  }
  
  grDevices::pdf(pdf_file, width = 8.5, height = 11, onefile = TRUE)
  
  n_total <- length(plot_list)
  
  for (i in seq(1, n_total, by = plots_per_page)) {
    idx <- i:min(i + plots_per_page - 1, n_total)
    
    gridExtra::grid.arrange(
      grobs = plot_list[idx],
      ncol = 1,   # stack cancer plots vertically
      nrow = length(idx)
    )
  }
  
  grDevices::dev.off()
  message(sprintf("[PDF SAVED] %s", pdf_file))
}

plot_all_cancers_2panel_pdf(
  target_vec = target_sheets,
  out_dir = out_dir,
  pretty_name_map = pretty_name_map,
  pdf_file = file.path(out_dir, "all_cancers_loris_survival_2panel.pdf"),
  profiles_keep = c("loris_0.1", "loris_0.9"),
  x_max = 38,
  plots_per_page = 3
)


########### Reference value table
make_reference_covariates_tex <- function(
    data_file = "./codes/extracted_cancer_data_by_type.xlsx",
    target_vec = target_sheets,
    loris_col_override = NULL,
    out_tex = file.path(out_dir, "reference_covariates.tex"),
    digits = 2
) {
  suppressPackageStartupMessages({
    library(readxl)
    library(dplyr)
    library(tidyr)
    library(knitr)
    library(kableExtra)
  })
  
  rows <- list()
  
  for (tname in target_vec) {
    df_raw <- readxl::read_excel(data_file, sheet = tname)
    df_raw <- cap_cancer_df(df_raw)
    df_raw <- as.data.frame(df_raw)
    
    proc_obj <- prep_df(df_raw)
    feature_cols <- proc_obj$feature_cols
    loris_col <- detect_loris_col(feature_cols, override = loris_col_override)
    
    mean_tbl <- lapply(feature_cols, function(v) {
      if (v == loris_col) return(NULL)
      
      x <- suppressWarnings(as.numeric(df_raw[[v]]))
      
      data.frame(
        covariate = v,
        mean_value = mean(x, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }) %>% bind_rows()
    
    mean_tbl$target_cancer <- tname
    mean_tbl$target_label <- ifelse(
      tname %in% names(pretty_name_map),
      pretty_name_map[[tname]],
      tname
    )
    
    rows[[length(rows) + 1]] <- mean_tbl
  }
  
  df_long <- bind_rows(rows)
  
  df_wide <- df_long %>%
    dplyr::select(target_label, covariate, mean_value) %>%
    tidyr::pivot_wider(names_from = covariate, values_from = mean_value) %>%
    arrange(target_label)
  
  num_cols <- setdiff(names(df_wide), "target_label")
  df_wide[num_cols] <- lapply(df_wide[num_cols], function(x) round(x, digits))
  
  tex_lines <- knitr::kable(
    df_wide,
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    linesep = "",
    escape = TRUE,
    col.names = c("Cancer", setdiff(names(df_wide), "target_label")),
    caption = "Cancer-specific reference covariate values used for prediction. For each cancer type, reference values are computed as the original-scale mean of each non-LORIS covariate in the target dataset. LORIS is omitted because survival curves are generated by varying LORIS separately.",
    label = "tab:reference_covariates"
  ) %>%
    kableExtra::kable_styling(
      latex_options = c("repeat_header", "scale_down"),
      font_size = 8
    )
  
  writeLines(tex_lines, out_tex)
  message("Wrote: ", out_tex)
  
  invisible(df_wide)
}

make_reference_covariates_tex(
  data_file = data_file,
  target_vec = target_sheets,
  loris_col_override = loris_col_override,
  out_tex = file.path(out_dir, "reference_covariates.tex"),
  digits = 0
)

