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

results_dir <- "./results_allTMB_new_hessian_lasso_new"
selected_file <- file.path(results_dir, "selected_sources_summary.csv")
data_file     <- "./codes/extracted_cancer_data_by_type.xlsx"

pilot_out_dir <- file.path(results_dir, "pilot_parameter_by_column_realdata")
dir.create(pilot_out_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# User controls
# =========================================================
c_grid <- c(0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256)
# c_grid <- c(0, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048)
eps_smooth <- 1e-4

lambda_zeta_default <- 0.005
lambda_eta_default  <- 0


# all sheets / target sheets
all_sheets <- paste0("cancer", c(1:10, 12:16))
cv_tbl <- read.csv("./results_allTMB_new_hessian_lasso_new/cv_metrics_summary_table_with_c_lambda_bic.csv",
                   stringsAsFactors = FALSE)
# same c-multiplier map as your old script
# c_map <- c(
#   cancer2  = 0.5,
#   cancer12 = 0.5,
#   cancer13 = 0.3,
#   cancer15 = 0.8
# )
c_map <- cv_tbl$c_mult
names(c_map) <- as.character(all_sheets)
c_default <- 1
target_sheets <- cv_tbl$target_cancer

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

# =========================================================
# Helper functions for inverse-Hessian transfer
# =========================================================
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
# =========================================================
# Build scalar basis functions from fit$basis$g_eval
# =========================================================
build_g_funcs_from_basis <- function(basis_obj, all_times) {
  g_eval <- basis_obj$g_eval
  
  G0 <- g_eval(all_times)
  if (is.null(dim(G0))) {
    G0 <- matrix(G0, ncol = 1)
  }
  p_hat <- ncol(G0)
  
  g_funcs <- lapply(seq_len(p_hat), function(j) {
    force(j)
    function(t) {
      out <- g_eval(t)
      if (is.null(dim(out))) out <- matrix(out, ncol = p_hat)
      as.numeric(out[, j])
    }
  })
  
  list(g_funcs = g_funcs, p_hat = p_hat)
}
regularize_hessian <- function(H, ridge = 1e-6) {
  H <- (H + t(H)) / 2
  H + ridge * diag(ncol(H))
}
# =========================================================
# Run one target cancer
# =========================================================
run_one_target_pilot <- function(tname, proc, selected_df) {
  
  message("===================================================")
  message(sprintf("[Pilot] target = %s", tname))
  
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
    message(sprintf("[Pilot] target=%s skipped: no selected sources", tname))
    return(NULL)
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
  
  # fit full SieveTL
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
  
  if (!conv_ok) {
    message(sprintf("[Pilot] target=%s skipped: full fit did not converge", tname))
    return(NULL)
  }
  
  # build target / source list objects
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
  gf_obj <- build_g_funcs_from_basis(fit$basis, all_times)
  g_funcs <- gf_obj$g_funcs
  p_hat <- gf_obj$p_hat
  L <- ncol(target_list$x)
  d_n <- p_hat + L
  
  time_grid_target <- seq(min(target_list$stime, na.rm = TRUE),
                          max(target_list$stime, na.rm = TRUE),
                          length.out = 500)
  
  pooled_stime <- c(target_list$stime, source_list[[1]]$stime)
  time_grid_pooled <- seq(min(pooled_stime, na.rm = TRUE),
                          max(pooled_stime, na.rm = TRUE),
                          length.out = 500)
  
  H_target <- compute_target_hessian_theta(
    data_k = target_list,
    gamma_hat = fit$gamma_hat,
    beta_hat = fit$beta_hat,
    g_funcs = g_funcs,
    time_grid = time_grid_target
  )
  
  H_pooled <- compute_pooled_hessian_theta(
    target = target_list,
    source_list = source_list,
    gamma_hat = fit$gamma_hat,
    beta_hat = fit$beta_hat,
    g_funcs = g_funcs,
    time_grid = time_grid_pooled
  )
  H_target <- regularize_hessian(H_target, ridge = 1e-6)
  H_pooled <- regularize_hessian(H_pooled, ridge = 1e-6)
  n_target <- nrow(target_list$x)
  
  tie_break_rule <- paste(
    "smallest criterion;",
    "ties broken by smallest residual;",
    "then smaller c"
  )
  
  full_rows <- vector("list", d_n)
  best_rows <- vector("list", d_n)
  
  for (j in seq_len(d_n)) {
    e_j <- rep(0, d_n)
    e_j[j] <- 1
    
    omega_A_j <- fit_inverse_hessian_initializer_unpenalized(H_pooled, j)
    l1_omegaA_j <- sum(abs(omega_A_j))
    
    column_rows <- vector("list", length(c_grid))
    
    tol_df <- 1e-6
    
    for (k in seq_along(c_grid)) {
      c_val <- c_grid[k]
      lambda_val <- c_val * sqrt(log(d_n) / n_target)
      
      delta_j <- fit_inverse_hessian_correction(
        H_target = H_target,
        omega_A = omega_A_j,
        j = j,
        lambda_l1 = lambda_val,
        eps = eps_smooth
      )
      
      omega_j <- omega_A_j + delta_j
      
      residual_j <- sum((H_target %*% omega_j - e_j)^2) / d_n
      
      l1_delta_j <- sum(abs(delta_j))
      l0_delta_j <- sum(abs(delta_j) > tol_df)
      
      l1_omega_j <- sum(abs(omega_j))
      l0_omega_j <- sum(abs(omega_j) > tol_df)
      
      criterion_j <- log((residual_j + 1e-12) / n_target) +
        (log(n_target) / n_target) * l0_delta_j
      
      column_rows[[k]] <- data.frame(
        target_cancer = tname,
        j = j,
        c = c_val,
        lambda = lambda_val,
        d_n = d_n,
        n_target = n_target,
        n_source = nrow(source_list[[1]]$x),
        
        residual = residual_j,
        
        l1_omegaA = l1_omegaA_j,
        l1_delta = l1_delta_j,
        l1_omega = l1_omega_j,
        
        l0_delta = l0_delta_j,
        l0_omega = l0_omega_j,
        
        criterion = criterion_j,
        
        mean_abs_omegaA = mean(abs(omega_A_j)),
        mean_abs_delta = mean(abs(delta_j)),
        max_abs_delta = max(abs(delta_j)),
        
        stringsAsFactors = FALSE
      )
    }
    
    column_table <- bind_rows(column_rows) %>%
      arrange(criterion, residual, c)
    
    full_rows[[j]] <- column_table
    best_rows[[j]] <- data.frame(
      target_cancer = tname,
      j = j,
      best_c = column_table$c[1],
      best_lambda_formula = "c * sqrt(log(d_n) / n_target)",
      best_lambda_value = column_table$lambda[1],
      best_criterion = column_table$criterion[1],
      best_residual = column_table$residual[1],
      best_l1_delta = column_table$l1_delta[1],
      best_l0_delta = column_table$l0_delta[1],
      best_l1_omega = column_table$l1_omega[1],
      best_l0_omega = column_table$l0_omega[1],
      tie_break_rule_used = tie_break_rule,
      stringsAsFactors = FALSE
    )
  }
  
  pilot_bic_tables_by_column <- bind_rows(full_rows)
  pilot_selected_c_by_column <- bind_rows(best_rows)
  
  selected_c_distribution <- pilot_selected_c_by_column %>%
    group_by(best_c) %>%
    summarise(
      freq_selected = n(),
      prop_selected = n() / nrow(pilot_selected_c_by_column),
      mean_selected_criterion = mean(best_criterion),
      mean_selected_residual = mean(best_residual),
      .groups = "drop"
    ) %>%
    arrange(desc(freq_selected), mean_selected_criterion, best_c)
  
  target_dir <- file.path(pilot_out_dir, tname)
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  
  bic_csv  <- file.path(target_dir, "pilot_bic_tables_by_column.csv")
  bic_rds  <- file.path(target_dir, "pilot_bic_tables_by_column.rds")
  best_csv <- file.path(target_dir, "pilot_selected_c_by_column.csv")
  dist_csv <- file.path(target_dir, "selected_c_distribution.csv")
  
  write.csv(pilot_bic_tables_by_column, bic_csv, row.names = FALSE)
  saveRDS(pilot_bic_tables_by_column, bic_rds)
  write.csv(pilot_selected_c_by_column, best_csv, row.names = FALSE)
  write.csv(selected_c_distribution, dist_csv, row.names = FALSE)
  
  message(sprintf("[Pilot] target=%s wrote %s", tname, best_csv))
  
  list(
    target_cancer = tname,
    full = pilot_bic_tables_by_column,
    best = pilot_selected_c_by_column,
    dist = selected_c_distribution
  )
}

# =========================================================
# Main
# =========================================================
main <- function() {
  
  cancers <- lapply(all_sheets, function(sh) {
    df <- read_excel(data_file, sheet = sh)
    df <- cap_cancer_df(df)
    na.omit(df)
  })
  names(cancers) <- all_sheets
  proc <- lapply(cancers, prep_df)
  
  selected_df <- read.csv(selected_file, stringsAsFactors = FALSE)
  
  all_best <- list()
  all_dist <- list()
  
  for (tname in target_sheets) {
    out <- tryCatch(
      run_one_target_pilot(tname, proc, selected_df),
      error = function(e) {
        message(sprintf("[Pilot] target=%s failed: %s", tname, conditionMessage(e)))
        NULL
      }
    )
    
    if (!is.null(out)) {
      all_best[[length(all_best) + 1]] <- out$best
      all_dist[[length(all_dist) + 1]] <- out$dist %>%
        mutate(target_cancer = out$target_cancer, .before = 1)
    }
  }
  
  if (length(all_best) > 0) {
    all_best_df <- bind_rows(all_best)
    write.csv(
      all_best_df,
      file.path(pilot_out_dir, "all_targets_selected_c_by_column.csv"),
      row.names = FALSE
    )
  }
  
  if (length(all_dist) > 0) {
    all_dist_df <- bind_rows(all_dist)
    write.csv(
      all_dist_df,
      file.path(pilot_out_dir, "all_targets_selected_c_distribution.csv"),
      row.names = FALSE
    )
  }
  
  message("[Pilot] done.")
}

if (sys.nframe() == 0) {
  main()
}
