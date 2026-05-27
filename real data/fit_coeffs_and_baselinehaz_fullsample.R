#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(pracma)
  library(reticulate)
  library(Rcpp)
  library(discSurv)
  library(DiscreteKL)
})

setwd('/Users/yuxisong/Library/CloudStorage/Box-Box/Cox-SieveTL/simulation/results/new_icb_response_result')
source("./codes/cv_source_selection.R")
source("/Users/yuxisong/Library/CloudStorage/Box-Box/Cox-SieveTL/simulation/results/new_icb_response_result/codes/cv_eval_selected_source_bic.R")
results_dir <- "./results_allTMB_new_hessian_lasso"
selected_file <- file.path(results_dir, "selected_sources_summary.csv")
data_file <- "./codes/extracted_cancer_data_by_type.xlsx"


# transcox_lr_vec <- c(1e-3, 2e-3, 3e-3, 4e-3)
# transcox_nsteps_vec <- c(200, 800)
# transcox_l1_vec <- c(0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10)
# transcox_l2_vec <- c(0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10)
# transcox_lr_vec <- c(1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 4e-3)
transcox_nsteps_vec <- c(200, 800)
transcox_lr_vec <- c(5e-4, 1e-3, 2e-3, 5e-3)
transcox_l1_vec <- c(5e-4, 1e-3, 0.005, 0.01, 2, 10)
transcox_l2_vec <- c(5e-4, 1e-3, 0.005, 0.01, 2, 10)
transcox_tuning_path <- file.path(results_dir, "transcox_bic_tuning.csv")
transcox_fallback <- list(
  best_lr = 1e-3,
  best_nsteps = 800,
  best_l1 = 0.01,
  best_l2 = 0.01
)

lambda_zeta_default <- 0.005
lambda_eta_default <- 0

cv_tbl <- read.csv("./results_allTMB_new_hessian_lasso/cv_metrics_summary_table_with_c_lambda_bic.csv", stringsAsFactors = FALSE)
cv_tbl2 <- cv_tbl %>%
  mutate(
    best_lambda_zeta = ifelse(is.na(best_lambda_zeta) | best_lambda_zeta == 0, lambda_zeta_default, best_lambda_zeta),
    best_lambda_eta  = ifelse(is.na(best_lambda_eta)  | best_lambda_eta  == 0, lambda_eta_default,  best_lambda_eta)
  )


# cv_tbl2 <- cv_tbl2[-13, ]
lambda_map <- setNames(
  lapply(seq_len(nrow(cv_tbl2)), function(i) {
    list(
      lambda_zeta = cv_tbl2$best_lambda_zeta[i],
      lambda_eta  = cv_tbl2$best_lambda_eta[i]
    )
  }),
  cv_tbl2$target_cancer
)

force_double_df <- function(df) {
  for (nm in names(df)) {
    if (is.numeric(df[[nm]]) || is.integer(df[[nm]]) || is.logical(df[[nm]])) {
      df[[nm]] <- as.numeric(df[[nm]])
    }
  }
  df
}

# All cancers remain candidate sources

all_sheets <- paste0("cancer", c(1:10, 12:16))
  # cv_tbl2$target_cancer
  # paste0("cancer", c(1:16, 18))

# Target cancers (predefined subset of 12)
target_sheets <- cv_tbl2$target_cancer
#   c(
#   "cancer1","cancer2","cancer3","cancer4","cancer5","cancer6","cancer7","cancer8",
#   "cancer9","cancer10","cancer11","cancer12","cancer13","cancer14","cancer15","cancer16","cancer18"
# )

select_transcox_by_bic <- function(target_df, source_df, cov, timevar, statusvar) {
  prim_raw <- target_df[, c(timevar, statusvar, cov), drop = FALSE]
  aux_raw <- source_df[, c(timevar, statusvar, cov), drop = FALSE]
  idx_prim <- complete.cases(prim_raw)
  idx_aux <- complete.cases(aux_raw)
  prim_raw <- prim_raw[idx_prim, , drop = FALSE]
  aux_raw <- aux_raw[idx_aux, , drop = FALSE]
  
  prim_raw[[statusvar]] <- as.numeric(prim_raw[[statusvar]])
  aux_raw[[statusvar]] <- as.numeric(aux_raw[[statusvar]])
  prim_raw[[statusvar]] <- ifelse(prim_raw[[statusvar]] == 1, 1, 2)
  aux_raw[[statusvar]] <- ifelse(aux_raw[[statusvar]] == 1, 1, 2)
  
  p <- length(cov)
  cov_names <- paste0("X", seq_len(p))
  prim_df <- data.frame(
    time = prim_raw[[timevar]],
    status = prim_raw[[statusvar]],
    prim_raw[, cov, drop = FALSE],
    check.names = FALSE
  )
  names(prim_df) <- c("time", "status", cov_names)
  aux_df <- data.frame(
    time = aux_raw[[timevar]],
    status = aux_raw[[statusvar]],
    aux_raw[, cov, drop = FALSE],
    check.names = FALSE
  )
  names(aux_df) <- c("time", "status", cov_names)
  
  message(sprintf("[TransCox BIC] prim nrow=%d time=%d status=%d",
                  nrow(prim_df), length(prim_df$time), length(prim_df$status)))
  message(sprintf("[TransCox BIC] aux nrow=%d time=%d status=%d",
                  nrow(aux_df), length(aux_df$time), length(aux_df$status)))
  if (length(prim_df$time) != length(prim_df$status) ||
      length(aux_df$time) != length(aux_df$status)) {
    stop("BIC data frame length mismatch.")
  }
  
  LRres <- SelLR_By_BIC(
    primData = prim_df,
    auxData = aux_df,
    cov = cov_names,
    statusvar = "status",
    lambda1 = 0.1,
    lambda2 = 0.1,
    learning_rate_vec = transcox_lr_vec,
    nsteps_vec = transcox_nsteps_vec
  )
  
  BICres <- SelParam_By_BIC(
    primData = prim_df,
    auxData = aux_df,
    cov = cov_names,
    statusvar = "status",
    lambda1_vec = transcox_l1_vec,
    lambda2_vec = transcox_l2_vec,
    learning_rate = LRres$best_lr,
    nsteps = LRres$best_nsteps
  )
  
  list(
    best_lr = LRres$best_lr,
    best_nsteps = LRres$best_nsteps,
    best_l1 = BICres$best_la1,
    best_l2 = BICres$best_la2,
    LRres = LRres,
    BICres = BICres
  )
}

use_condaenv("TransCoxEnvi", required = TRUE)
# tf <- import("tensorflow")
source_python(system.file("python", "TransCoxFunction.py", package = "TransCox"))
library(TransCox)

Rcpp::sourceCpp("/Users/yuxisong/Library/CloudStorage/Box-Box/Cox-SieveTL/codes/Discrete_KL_eigen_cloglog.cpp")
if (!exists("contToDisc")) stop("contToDisc not found: please install/load discSurv")

read_optional_csv <- function(path) {
  exists_flag <- file.exists(path)
  message(sprintf("[EXTCSV] exists=%s path=%s", exists_flag, path))
  if (!exists_flag) return(NULL)
  df <- tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df)) {
    message(sprintf("[EXTCSV] read failed path=%s", path))
    return(NULL)
  }
  message(sprintf("[EXTCSV] nrow=%s", nrow(df)))
  message(sprintf("[EXTCSV] colnames=%s", paste(colnames(df), collapse = ",")))
  if ("method" %in% names(df)) {
    message(sprintf("[EXTCSV] unique(method)=%s",
                    paste(unique(df$method), collapse = ",")))
  } else {
    message("[EXTCSV] unique(method)=<missing column>")
  }
  df
}
soft_thresh <- function(z, a) {
  sign(z) * pmax(abs(z) - a, 0)
}

compute_hessian_theta <- function(data, gamma_hat, beta_hat, g_eval, time_grid) {
  n <- nrow(data$x)
  L <- length(beta_hat)
  p <- length(gamma_hat)
  d <- p + L
  
  H <- matrix(0, d, d)
  g_time_mat <- g_eval(time_grid)
  
  for (i in seq_len(n)) {
    x_i <- data$x[i, ]
    t_i <- data$stime[i]
    
    idx <- which(time_grid <= t_i)
    if (length(idx) < 2) next
    
    tg <- time_grid[idx]
    G_i <- g_time_mat[idx, , drop = FALSE]
    m <- length(tg)
    
    exp_eta <- exp(as.vector(G_i %*% gamma_hat) + sum(x_i * beta_hat))
    
    Z_mat <- cbind(
      G_i,
      matrix(rep(x_i, each = m), nrow = m)
    )
    
    for (a in seq_len(d)) {
      for (b in seq_len(d)) {
        integrand <- Z_mat[, a] * Z_mat[, b] * exp_eta
        H[a, b] <- H[a, b] + pracma::trapz(tg, integrand)
      }
    }
  }
  
  H / n
}

run_one_sparse_hessian_update <- function(delta_old,
                                          lambda_zeta, lambda_eta,
                                          L, p,
                                          target,
                                          G_mat_target,
                                          time_grid_target,
                                          g_time_mat_target,
                                          gamma_hat_A,
                                          beta_hat_A,
                                          g_eval,
                                          ridge_iter = 1e-6) {
  gamma_old <- gamma_hat_A + delta_old[1:p]
  beta_old  <- beta_hat_A  + delta_old[(p + 1):(p + L)]
  
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
  
  H <- compute_hessian_theta(
    data = target,
    gamma_hat = gamma_old,
    beta_hat = beta_old,
    g_eval = g_eval,
    time_grid = time_grid_target
  )
  
  d <- length(delta_old)
  H_reg <- H + ridge_iter * diag(d)
  
  lambda_vec <- c(rep(lambda_zeta, p), rep(lambda_eta, L))
  Hdiag <- pmax(diag(H_reg), 1e-8)
  
  z <- delta_old - grad / Hdiag
  delta_new <- soft_thresh(z, lambda_vec / Hdiag)
  
  theta_new <- c(gamma_hat_A, beta_hat_A) + delta_new
  
  list(
    delta_hat = delta_new,
    gamma_hat = theta_new[1:p],
    beta_hat = theta_new[(p + 1):(p + L)],
    H_lasso = H,
    grad_lasso = grad
  )
}
fit_algorithm1 <- function(target_train_df, source_df = NULL,
                           feature_cols, time_col = "time", event_col = "event",
                           lambda_zeta = 0.05, lambda_eta = 0.05,
                           c_mult = 1, maxit = 5000,
                           use_hessian_lasso_update = TRUE,
                           ridge_iter = 1e-6) {
  target_list <- list(
    x = as.matrix(target_train_df[, feature_cols, drop = FALSE]),
    stime = as.numeric(target_train_df[[time_col]]),
    type = as.numeric(target_train_df[[event_col]])
  )
  
  source_list <- NULL
  if (!is.null(source_df)) {
    source_list <- list(
      list(
        x = as.matrix(source_df[, feature_cols, drop = FALSE]),
        stime = as.numeric(source_df[[time_col]]),
        type = as.numeric(source_df[[event_col]])
      )
    )
  }
  
  all_time <- if (is.null(source_list)) {
    target_list$stime
  } else {
    c(target_list$stime, source_list[[1]]$stime)
  }
  
  if (!any(is.finite(all_time))) stop("No finite times for spline basis.")
  
  basis <- build_basis_quantile(all_time, degree = 3, c_mult = c_mult)
  g_funcs <- function(t) basis$g_eval(t)
  
  all.x <- if (is.null(source_list)) {
    target_list$x
  } else {
    rbind(target_list$x, source_list[[1]]$x)
  }
  
  all.stime <- if (is.null(source_list)) {
    target_list$stime
  } else {
    c(target_list$stime, source_list[[1]]$stime)
  }
  
  all.type <- if (is.null(source_list)) {
    target_list$type
  } else {
    c(target_list$type, source_list[[1]]$type)
  }
  
  G_mat <- g_funcs(all.stime)
  p <- ncol(G_mat)
  L <- ncol(all.x)
  
  time_grid <- seq(min(all_time), max(all_time), length.out = 500)
  g_time_mat <- g_funcs(time_grid)
  
  theta0 <- rep(0, L + p)
  
  loss_wrap <- function(th) {
    transfer_loss_fast(
      th, p, L,
      all.x, all.stime, all.type,
      G_mat, time_grid, g_time_mat
    )
  }
  
  grad_wrap <- function(th) {
    transfer_grad_fast(
      th, p, L,
      all.x, all.stime, all.type,
      G_mat, time_grid, g_time_mat
    )
  }
  
  fit <- tryCatch(
    optim(
      par = theta0,
      fn = loss_wrap,
      gr = grad_wrap,
      method = "L-BFGS-B",
      control = list(maxit = maxit, factr = 1e-10, pgtol = 1e-6)
    ),
    error = function(e) e
  )
  
  converged_transfer <- as.integer(
    is.list(fit) &&
      !is.null(fit$convergence) &&
      fit$convergence == 0 &&
      is.finite(fit$value)
  )
  
  if (converged_transfer == 0L) {
    return(list(conv_transfer = 0L, conv_debias = 0L))
  }
  
  theta_hat_A <- fit$par
  gamma_hat_A <- theta_hat_A[1:p]
  beta_hat_A  <- theta_hat_A[(p + 1):(L + p)]
  
  G_mat_target <- g_funcs(target_list$stime)
  time_grid_target <- seq(min(target_list$stime), max(target_list$stime), length.out = 500)
  g_time_mat_target <- g_funcs(time_grid_target)
  
  delta0 <- rep(0, L + p)
  
  deb_loss <- function(d) {
    debias_loss_fast(
      d, lambda_zeta, lambda_eta, L, p,
      target_list, G_mat_target, time_grid_target,
      g_time_mat_target, gamma_hat_A, beta_hat_A
    )
  }
  
  deb_grad <- function(d) {
    debias_grad_fast(
      d, lambda_zeta, lambda_eta, L, p,
      target_list, G_mat_target, time_grid_target,
      g_time_mat_target, gamma_hat_A, beta_hat_A
    )
  }
  
  fit_debias <- tryCatch(
    optim(
      par = delta0,
      fn = deb_loss,
      gr = deb_grad,
      method = "L-BFGS-B",
      control = list(maxit = maxit, factr = 1e-10, pgtol = 1e-6)
    ),
    error = function(e) e
  )
  
  converged_debias <- as.integer(
    is.list(fit_debias) &&
      !is.null(fit_debias$convergence) &&
      fit_debias$convergence == 0 &&
      is.finite(fit_debias$value)
  )
  
  if (converged_debias == 0L) {
    return(list(conv_transfer = converged_transfer, conv_debias = 0L))
  }
  
  delta_hat_optim <- fit_debias$par
  
  if (use_hessian_lasso_update) {
    sparse_fit <- run_one_sparse_hessian_update(
      delta_old = delta_hat_optim,
      lambda_zeta = lambda_zeta,
      lambda_eta = lambda_eta,
      L = L,
      p = p,
      target = target_list,
      G_mat_target = G_mat_target,
      time_grid_target = time_grid_target,
      g_time_mat_target = g_time_mat_target,
      gamma_hat_A = gamma_hat_A,
      beta_hat_A = beta_hat_A,
      g_eval = basis$g_eval,
      ridge_iter = ridge_iter
    )
    
    delta_hat <- sparse_fit$delta_hat
    gamma_hat <- sparse_fit$gamma_hat
    beta_hat  <- sparse_fit$beta_hat
    H_lasso <- sparse_fit$H_lasso
    grad_lasso <- sparse_fit$grad_lasso
  } else {
    delta_hat <- delta_hat_optim
    gamma_hat <- gamma_hat_A + delta_hat[1:p]
    beta_hat  <- beta_hat_A  + delta_hat[(p + 1):(p + L)]
    H_lasso <- NULL
    grad_lasso <- NULL
  }
  
  zeta_hat <- delta_hat[1:p]
  eta_hat  <- delta_hat[(p + 1):(p + L)]
  
  list(
    g_eval = basis$g_eval,
    gamma_hat_A = gamma_hat_A,
    beta_hat_A = beta_hat_A,
    
    delta_hat_optim = delta_hat_optim,
    delta_hat = delta_hat,
    zeta_hat = zeta_hat,
    eta_hat = eta_hat,
    
    gamma_hat = gamma_hat,
    beta_hat = beta_hat,
    
    H_lasso = H_lasso,
    grad_lasso = grad_lasso,
    n_nonzero_delta = sum(abs(delta_hat) > 1e-6),
    n_nonzero_zeta = sum(abs(zeta_hat) > 1e-6),
    n_nonzero_eta = sum(abs(eta_hat) > 1e-6),
    
    time_grid_fit = time_grid,
    g_time_mat_fit = g_time_mat,
    conv_transfer = converged_transfer,
    conv_debias = converged_debias,
    basis = basis,
    use_hessian_lasso_update = use_hessian_lasso_update
  )
}

fit_algorithm1_transfer_only <- function(train_df, source_df, feature_cols, time_col, event_col,
                                         lambda_zeta, lambda_eta,
                                         maxit = 6000, c_mult = 1) {
  fit <- fit_algorithm1(
    train_df,
    source_df = source_df,
    feature_cols = feature_cols,
    time_col = time_col,
    event_col = event_col,
    lambda_zeta = lambda_zeta,
    lambda_eta = lambda_eta,
    maxit = maxit,
    c_mult = c_mult,
    use_hessian_lasso_update = FALSE
  )
  
  if (!is.list(fit) || is.null(fit$conv_transfer) || fit$conv_transfer != 1L) {
    return(list(conv_transfer = 0L, conv_debias = 0L))
  }
  
  list(
    beta_hat = fit$beta_hat_A,
    gamma_hat = fit$gamma_hat_A,
    g_eval = fit$g_eval,
    time_grid_fit = fit$time_grid_fit,
    g_time_mat_fit = fit$g_time_mat_fit,
    basis = fit$basis,
    conv_transfer = fit$conv_transfer,
    conv_debias = 1L
  )
}

compute_baseline <- function(fit) {
  if (is.null(fit$gamma_hat) || is.null(fit$g_eval)) return(NULL)
  if (!is.null(fit$time_grid_fit)) {
    t_grid <- fit$time_grid_fit
    g_time <- fit$g_time_mat_fit
  } else {
    return(NULL)
  }
  h0 <- as.vector(exp(g_time %*% fit$gamma_hat))
  H0 <- pracma::cumtrapz(t_grid, h0)
  S0 <- exp(-H0)
  tibble(t = t_grid, h0 = h0, H0 = H0, S0 = S0)
}

fit_transcox_fullsample <- function(target_df, source_df_pool, feature_cols, time_col, event_col,
                                    l1, l2, learning_rate, nsteps, diag_target = NULL) {
  if (is.null(source_df_pool) || nrow(source_df_pool) == 0) {
    return(list(conv = 0L, reason = "empty_source"))
  }
  p_tc <- target_df
  a_tc <- source_df_pool
  names(p_tc)[names(p_tc) == time_col] <- "time"
  names(p_tc)[names(p_tc) == event_col] <- "status"
  names(a_tc)[names(a_tc) == time_col] <- "time"
  names(a_tc)[names(a_tc) == event_col] <- "status"
  p_tc$status <- as.numeric(p_tc$status)
  a_tc$status <- as.numeric(a_tc$status)
  if (sum(target_df[[event_col]] == 1, na.rm = TRUE) < 2 ||
      sum(source_df_pool[[event_col]] == 1, na.rm = TRUE) < 2) {
    return(list(conv = 0L, reason = "too_few_events"))
  }
  p_tc$status <- ifelse(p_tc$status == 1, 1L, 2L)
  a_tc$status <- ifelse(a_tc$status == 1, 1L, 2L)
  x_names <- paste0("X", seq_along(feature_cols))
  for (i in seq_along(feature_cols)) {
    col <- feature_cols[i]
    if (col %in% names(p_tc)) names(p_tc)[names(p_tc) == col] <- x_names[i]
    if (col %in% names(a_tc)) names(a_tc)[names(a_tc) == col] <- x_names[i]
  }

  out <- tryCatch({
    Cout <- GetAuxSurv(a_tc, cov = x_names)
    if (is.null(Cout$q) || is.null(Cout$q$breakPoints) || is.null(Cout$q$cumHazards)) {
      return(list(conv = 0L, reason = "missing_Cout_q"))
    }
    if (is.null(Cout$estR)) {
      return(list(conv = 0L, reason = "missing_Cout_estR"))
    }
    Pout <- GetPrimaryParam(p_tc, q = Cout$q, estR = Cout$estR)
    if (!is.null(diag_target) && diag_target == "cancer11") {
      diag_cov_stats <- function(df, label) {
        message(sprintf("[TransCox C11] covariate stats (%s)", label))
        for (col in feature_cols) {
          x <- df[[col]]
          if (is.null(x)) next
          qs <- stats::quantile(x, probs = c(0.01, 0.5, 0.99), na.rm = TRUE, names = FALSE, type = 7)
          message(sprintf("[TransCox C11] %s %s: min=%g q01=%g med=%g q99=%g max=%g",
                          label, col, min(x, na.rm = TRUE), qs[1], qs[2], qs[3], max(x, na.rm = TRUE)))
        }
      }
      diag_cov_stats(target_df, "target")
      diag_cov_stats(source_df_pool, "source")
    }
    Tres <- runTransCox_one(Pout, l1 = l1, l2 = l2, learning_rate = learning_rate, nsteps = nsteps,
                            cov = x_names)
    if (!is.null(diag_target) && diag_target == "cancer11") {
      new_beta_dbg <- NULL
      if (!is.null(Tres$new_beta)) {
        new_beta_dbg <- as.numeric(Tres$new_beta)
      } else if (!is.null(Tres$beta)) {
        new_beta_dbg <- as.numeric(Tres$beta)
      }
      if (!is.null(new_beta_dbg)) {
        message(sprintf("[TransCox C11] new_beta head=%s",
                        paste(head(new_beta_dbg), collapse = ",")))
      } else {
        message("[TransCox C11] new_beta head=<missing>")
      }
      if (!is.null(Tres$eta)) {
        message(sprintf("[TransCox C11] eta head=%s",
                        paste(head(as.numeric(Tres$eta)), collapse = ",")))
      } else {
        message("[TransCox C11] eta head=<missing>")
      }
      if (!is.null(Tres$xi)) {
        message(sprintf("[TransCox C11] xi head=%s",
                        paste(head(as.numeric(Tres$xi)), collapse = ",")))
      } else {
        message("[TransCox C11] xi head=<missing>")
      }

      cox_fit <- tryCatch({
        survival::coxph(
          as.formula(paste0("survival::Surv(", time_col, ",", event_col, ") ~ ",
                            paste(feature_cols, collapse = " + "))),
          data = target_df
        )
      }, error = function(e) NULL)
      if (!is.null(cox_fit) && !is.null(new_beta_dbg) &&
          length(new_beta_dbg) == length(feature_cols)) {
        beta_cox <- stats::coef(cox_fit)
        beta_cox <- beta_cox[feature_cols]
        if (all(is.finite(beta_cox))) {
          l2_dist <- sqrt(sum((new_beta_dbg - beta_cox) ^ 2))
          message(sprintf("[TransCox C11] L2(new_beta, coxph_target_beta)=%g", l2_dist))
        } else {
          message("[TransCox C11] L2(new_beta, coxph_target_beta)=<nonfinite>")
        }
      } else {
        message("[TransCox C11] L2(new_beta, coxph_target_beta)=<unavailable>")
      }
    }
    # ===== DEBUG: diagnose new_IntH from TransCox =====
    message(sprintf(
      "[TransCox DEBUG] time range = [%g, %g]",
      min(Tres$time, na.rm = TRUE),
      max(Tres$time, na.rm = TRUE)
    ))

    message(sprintf(
      "[TransCox DEBUG] new_IntH summary: min=%g, max=%g, n_neg=%d / %d",
      min(Tres$new_IntH, na.rm = TRUE),
      max(Tres$new_IntH, na.rm = TRUE),
      sum(Tres$new_IntH < 0, na.rm = TRUE),
      length(Tres$new_IntH)
    ))

    tmp_dbg <- data.frame(
      time = as.numeric(Tres$time),
      x = as.numeric(Tres$new_IntH)
    )
    tmp_dbg <- tmp_dbg[is.finite(tmp_dbg$time) & is.finite(tmp_dbg$x), ]
    tmp_dbg <- tmp_dbg[order(tmp_dbg$time), ]
    tmp_dbg <- tmp_dbg[!duplicated(tmp_dbg$time), ]

    H_A <- tmp_dbg$x
    H_B <- cumsum(tmp_dbg$x)

    message(sprintf(
      "[TransCox DEBUG] H_A (no cumsum): monotone=%s, range=[%g,%g]",
      all(diff(H_A) >= -1e-10),
      min(H_A), max(H_A)
    ))

    message(sprintf(
      "[TransCox DEBUG] H_B (cumsum): monotone=%s, range=[%g,%g]",
      all(diff(H_B) >= -1e-10),
      min(H_B), max(H_B)
    ))
    # ===== END DEBUG =====
    beta_hat <- NULL
    if (!is.null(Tres$new_beta)) {
      beta_hat <- as.numeric(Tres$new_beta)
    } else if (!is.null(Tres$beta)) {
      beta_hat <- as.numeric(Tres$beta)
    } else {
      for (nm in names(Tres)) {
        v <- Tres[[nm]]
        if (is.numeric(v) && length(v) == length(feature_cols)) {
          beta_hat <- as.numeric(v)
          break
        }
      }
    }
    if (is.null(beta_hat) || length(beta_hat) != length(feature_cols)) {
      return(list(conv = 0L, reason = "beta_missing"))
    }
    names(beta_hat) <- feature_cols
    if (is.null(Tres$time) || is.null(Tres$new_IntH)) {
      return(list(conv = 0L, reason = "missing_time_intH"))
    }
    df_dbg <- data.frame(
      time = as.numeric(Tres$time),
      new_IntH_raw = as.numeric(Tres$new_IntH)
    )
    if (any(!is.finite(df_dbg$time)) || any(!is.finite(df_dbg$new_IntH_raw))) {
      return(list(conv = 0L, reason = "nonfinite_time_intH"))
    }
    df_dbg <- df_dbg[order(df_dbg$time), , drop = FALSE]
    df_dbg <- df_dbg[!duplicated(df_dbg$time), , drop = FALSE]
    if (nrow(df_dbg) < 2) {
      return(list(conv = 0L, reason = "H0_insufficient_points"))
    }

    dH <- df_dbg$new_IntH_raw
    message(sprintf("[TransCox] target=%s frac_neg_dH=%g",
                    if (!is.null(diag_target)) diag_target else "<unknown>",
                    mean(dH < 0, na.rm = TRUE)))

    H0_trans_df <- data.frame(
      time = df_dbg$time,
      new_IntH_raw = dH,
      H0_raw_cumsum = cumsum(dH)
    )
    list(conv = 1L, reason = "ok", beta_hat = beta_hat, H0_df = H0_trans_df)
  }, error = function(e) {
    list(conv = 0L, reason = paste0("error: ", conditionMessage(e)))
  })

  out
}

compute_transcox_outputs <- function() {
  cancers <- lapply(all_sheets, function(sh) {
    df <- read_excel(data_file, sheet = sh)
    df <- cap_cancer_df(df)
    na.omit(df)
  })
  names(cancers) <- all_sheets
  proc <- lapply(cancers, prep_df)

  selected_df <- read.csv(selected_file, stringsAsFactors = FALSE)

  tuning_log <- if (file.exists(transcox_tuning_path)) {
    read.csv(transcox_tuning_path, stringsAsFactors = FALSE)
  } else {
    data.frame(
      target_cancer = character(0),
      best_lr = numeric(0),
      best_nsteps = numeric(0),
      best_l1 = numeric(0),
      best_l2 = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  beta_rows <- list()
  base_rows <- list()

  for (tname in target_sheets) {
    target <- proc[[tname]]
    target_df_raw <- target$df

    sel_row <- selected_df[selected_df$target_cancer == tname, , drop = FALSE]
    sel_sources <- character(0)
    if (nrow(sel_row) > 0 &&
        !is.na(sel_row$n_selected[1]) && sel_row$n_selected[1] > 0) {
      sel_sources <- parse_sources(sel_row$selected_sources[1])
      sel_sources <- sel_sources[sel_sources != "" & !is.na(sel_sources)]
      sel_sources <- sel_sources[sel_sources %in% names(proc) & sel_sources != tname]
    }

    source_df_pool <- NULL
    if (length(sel_sources) > 0) {
      source_df_pool <- bind_rows(lapply(sel_sources, function(sn) proc[[sn]]$df))
    }

    apply_list <- list()
    if (!is.null(source_df_pool)) {
      apply_list$src <- source_df_pool
    }
    pre <- preprocess_features(target_df_raw, apply_list, target$feature_cols, target$time_col, target$event_col)
    target_df <- pre$train_df
    if (!is.null(source_df_pool)) {
      source_df_pool <- pre$apply_df_list$src
    }

    bic_sel <- NULL
    if (!is.null(source_df_pool) && nrow(source_df_pool) > 0) {
      cached <- tuning_log[tuning_log$target_cancer == tname, , drop = FALSE]
      if (nrow(cached) == 1 &&
          all(c("best_lr", "best_nsteps", "best_l1", "best_l2") %in% names(cached)) &&
          all(is.finite(cached$best_lr), is.finite(cached$best_nsteps),
              is.finite(cached$best_l1), is.finite(cached$best_l2))) {
        bic_sel <- list(
          best_lr = cached$best_lr,
          best_nsteps = cached$best_nsteps,
          best_l1 = cached$best_l1,
          best_l2 = cached$best_l2
        )
      } else {
        # bic_sel <- select_transcox_by_bic(
        #   target_df = target_df,
        #   source_df = source_df_pool,
        #   cov = target$feature_cols,
        #   timevar = target$time_col,
        #   statusvar = target$event_col
        # )
        bic_sel <- tryCatch(
          select_transcox_by_bic(
            target_df = target_df,
            source_df = source_df_pool,
            cov = target$feature_cols,
            timevar = target$time_col,
            statusvar = target$event_col
          ),
          error = function(e) {
            message(sprintf("[TransCox] BIC selection failed target=%s err=%s", tname, conditionMessage(e)))
            NULL
          }
        )
        if (!is.null(bic_sel)) {
          tuning_log <- rbind(
            tuning_log,
            data.frame(
              target_cancer = tname,
              best_lr = bic_sel$best_lr,
              best_nsteps = bic_sel$best_nsteps,
              best_l1 = bic_sel$best_l1,
              best_l2 = bic_sel$best_l2,
              stringsAsFactors = FALSE
            )
          )
        }
      }
    }

    # if (is.null(bic_sel)) {
    #   message(sprintf("[TransCox] BIC failed -> writing NA (target=%s)", tname))
    #   beta_rows[[length(beta_rows) + 1]] <- tibble(
    #     target_cancer = tname,
    #     method = "transcox"
    #   ) %>% bind_cols(as_tibble(setNames(as.list(rep(NA_real_, length(target$feature_cols))),
    #                                      target$feature_cols)))
    #   base_rows[[length(base_rows) + 1]] <- tibble(
    #     target_cancer = tname,
    #     method = "transcox",
    #     t = NA_real_,
    #     h0 = NA_real_,
    #     H0 = NA_real_,
    #     S0 = NA_real_
    #   )
    #   next
    # }
    if (is.null(bic_sel)) {
      message(sprintf(
        "[TransCox] BIC failed -> using fallback values for target=%s: lr=%g, nsteps=%d, l1=%g, l2=%g",
        tname,
        transcox_fallback$best_lr,
        transcox_fallback$best_nsteps,
        transcox_fallback$best_l1,
        transcox_fallback$best_l2
      ))
      bic_sel <- transcox_fallback
    }

    fit <- fit_transcox_fullsample(target_df, source_df_pool,
                                   feature_cols = target$feature_cols,
                                   time_col = target$time_col,
                                   event_col = target$event_col,
                                   l1 = bic_sel$best_l1,
                                   l2 = bic_sel$best_l2,
                                   learning_rate = bic_sel$best_lr,
                                   nsteps = bic_sel$best_nsteps,
                                   diag_target = tname)
    # fit_test <- fit_transcox_fullsample(
    #   target_df, source_df_pool,
    #   feature_cols= target$feature_cols,
    #   time_col= target$time_col,
    #   event_col= target$event_col,
    #   l1 = 1e2,
    #   l2 = 1e2,
    #   learning_rate = 1e-6,
    #   nsteps = 800
    # )
    message(sprintf("[TransCox] target=%s conv=%d reason=%s", tname, fit$conv, fit$reason))

    if (is.list(fit) && !is.null(fit$conv) && fit$conv == 1L) {
      beta_vals <- rep(NA_real_, length(target$feature_cols))
      names(beta_vals) <- target$feature_cols
      if (!is.null(fit$beta_hat)) {
        beta_vals[names(fit$beta_hat)] <- as.numeric(fit$beta_hat)
      }
      beta_rows[[length(beta_rows) + 1]] <- tibble(
        target_cancer = tname,
        method = "transcox"
      ) %>% bind_cols(as_tibble(as.list(beta_vals)))

      df_native <- data.frame(
        time = fit$H0_df$time,
        new_IntH_raw = fit$H0_df$new_IntH_raw,
        H0_raw_cumsum = fit$H0_df$H0_raw_cumsum
      )
      df_native <- df_native[is.finite(df_native$time) & is.finite(df_native$H0_raw_cumsum), , drop = FALSE]
      df_native <- df_native[order(df_native$time), , drop = FALSE]
      if (nrow(df_native) >= 2) {
        H0 <- df_native$H0_raw_cumsum
        S0 <- exp(-H0)
        base_rows[[length(base_rows) + 1]] <- tibble(
          target_cancer = tname,
          method = "transcox",
          t = df_native$time,
          h0 = NA_real_,
          H0 = H0,
          S0 = S0,
          new_IntH_raw = df_native$new_IntH_raw,
          H0_raw_cumsum = df_native$H0_raw_cumsum
        )
      }
    } else {
      beta_rows[[length(beta_rows) + 1]] <- tibble(
        target_cancer = tname,
        method = "transcox"
      ) %>% bind_cols(as_tibble(setNames(as.list(rep(NA_real_, length(target$feature_cols))),
                                         target$feature_cols)))
    }
  }

  beta_df <- bind_rows(beta_rows)
  base_df_out <- bind_rows(base_rows)

  write.csv(beta_df, file.path(results_dir, "transcox_beta_estimates.csv"), row.names = FALSE)
  write.csv(base_df_out, file.path(results_dir, "transcox_baseline_hazard.csv"), row.names = FALSE)
  if (nrow(tuning_log) > 0) {
    tuning_log <- tuning_log[!duplicated(tuning_log$target_cancer, fromLast = TRUE), , drop = FALSE]
    write.csv(tuning_log, transcox_tuning_path, row.names = FALSE)
  }
  message(sprintf("[TransCox] wrote %s", file.path(results_dir, "transcox_beta_estimates.csv")))
  message(sprintf("[TransCox] wrote %s", file.path(results_dir, "transcox_baseline_hazard.csv")))
}

fit_discretekl_fullsample <- function(target_df, source_df_pool, feature_cols, time_col, event_col) {
  if (is.null(source_df_pool) || nrow(source_df_pool) == 0) {
    return(list(conv = 0L, reason = "empty_source"))
  }

  train <- target_df
  aData <- source_df_pool
  if (time_col != "time") {
    names(train)[names(train) == time_col] <- "time"
    names(aData)[names(aData) == time_col] <- "time"
  }
  if (event_col != "status") {
    names(train)[names(train) == event_col] <- "status"
    names(aData)[names(aData) == event_col] <- "status"
  }

  train$status <- ifelse(as.numeric(train$status) == 1, 1L, 0L)
  aData$status <- ifelse(as.numeric(aData$status) == 1, 1L, 0L)

  time_grid_kl <- as.numeric(stats::quantile(train$time, probs = seq(0, 1, length.out = 20), na.rm = TRUE))
  time_grid_kl <- sort(unique(time_grid_kl[is.finite(time_grid_kl)]))
  if (length(time_grid_kl) < 2) {
    return(list(conv = 0L, reason = "time_grid_insufficient"))
  }

  out <- tryCatch({
    t_disc_source <- contToDisc(
      dataShort = aData,
      timeColumn = "time",
      intervalLimits = time_grid_kl
    )$timeDisc
    X_source <- dplyr::select(aData, -c("time", "status"))
    keep_source <- !is.na(t_disc_source)
    if (sum(keep_source) < 2) {
      return(list(conv = 0L, reason = "source_disc_insufficient"))
    }
    prior_fit <- discSurv_cloglog(
      t_disc_source[keep_source],
      X_source[keep_source, , drop = FALSE],
      aData$status[keep_source]
    )

    t_disc_target <- contToDisc(
      dataShort = train,
      timeColumn = "time",
      intervalLimits = time_grid_kl
    )$timeDisc
    X_target <- dplyr::select(train, -c("time", "status"))
    df_input <- cbind(X_target, time = t_disc_target, status = train$status)
    df_input <- df_input[complete.cases(df_input), , drop = FALSE]
    if (nrow(df_input) < 2) {
      return(list(conv = 0L, reason = "target_input_insufficient"))
    }

    eta_min <- 1
    eta_max <- 1
    eta_interval <- 0.01
    kl_result <- kl_cloglog(
      day_prior = prior_fit$beta_t,
      beta_prior = prior_fit$beta_v,
      eta_min = eta_min,
      eta_max = eta_max,
      eta_interval = eta_interval,
      df_input = df_input
    )
    beta_v <- kl_result$model$beta_v
    beta_t <- as.numeric(kl_result$model$beta_t)
    if (is.null(beta_v) || length(beta_v) != length(feature_cols)) {
      return(list(conv = 0L, reason = "beta_missing"))
    }
    names(beta_v) <- feature_cols

    beta_t_cap <- pmin(pmax(beta_t, -50), 50)
    haz <- -expm1(-exp(beta_t_cap))
    haz <- pmin(pmax(haz, 0), 1)
    cumhaz <- cumsum(haz)
    basehaz_df <- data.frame(time = time_grid_kl, hazard = cumhaz)
    basehaz_df <- basehaz_df[is.finite(basehaz_df$time) & is.finite(basehaz_df$hazard), , drop = FALSE]
    basehaz_df <- basehaz_df[order(basehaz_df$time), , drop = FALSE]
    basehaz_df <- basehaz_df[!duplicated(basehaz_df$time), , drop = FALSE]
    if (nrow(basehaz_df) < 2) {
      return(list(conv = 0L, reason = "H0_insufficient_points"))
    }

    list(conv = 1L, reason = "ok", beta_hat = beta_v, H0_df = basehaz_df)
  }, error = function(e) {
    list(conv = 0L, reason = paste0("error: ", conditionMessage(e)))
  })

  out
}

compute_discretekl_outputs <- function() {
  cancers <- lapply(all_sheets, function(sh) {
    df <- read_excel(data_file, sheet = sh)
    df <- cap_cancer_df(df)
    na.omit(df)
  })
  names(cancers) <- all_sheets
  proc <- lapply(cancers, prep_df)

  selected_df <- read.csv(selected_file, stringsAsFactors = FALSE)

  beta_rows <- list()
  base_rows <- list()

  for (tname in target_sheets) {
    target <- proc[[tname]]
    target_df <- target$df
    sel_row <- selected_df[selected_df$target_cancer == tname, , drop = FALSE]
    sel_sources <- character(0)
    if (nrow(sel_row) > 0 &&
        !is.na(sel_row$n_selected[1]) && sel_row$n_selected[1] > 0) {
      sel_sources <- parse_sources(sel_row$selected_sources[1])
      sel_sources <- sel_sources[sel_sources != "" & !is.na(sel_sources)]
      sel_sources <- sel_sources[sel_sources %in% names(proc) & sel_sources != tname]
    }

    source_df_pool <- NULL
    if (length(sel_sources) > 0) {
      source_df_pool <- bind_rows(lapply(sel_sources, function(sn) proc[[sn]]$df))
    }

    apply_list <- list()
    if (!is.null(source_df_pool)) {
      apply_list$src <- source_df_pool
    }
    pre <- preprocess_features(target_df, apply_list, target$feature_cols, target$time_col, target$event_col)
    target_df <- pre$train_df
    if (!is.null(source_df_pool)) {
      source_df_pool <- pre$apply_df_list$src
    }

    fit <- fit_discretekl_fullsample(target_df, source_df_pool,
                                     feature_cols = target$feature_cols,
                                     time_col = target$time_col,
                                     event_col = target$event_col)
    message(sprintf("[DiscreteKL] target=%s conv=%d reason=%s", tname, fit$conv, fit$reason))

    if (is.list(fit) && !is.null(fit$conv) && fit$conv == 1L) {
      beta_vals <- as.numeric(fit$beta_hat)
      names(beta_vals) <- target$feature_cols
      beta_rows[[length(beta_rows) + 1]] <- tibble(
        target_cancer = tname,
        method = "discretekl"
      ) %>% bind_cols(as_tibble(as.list(beta_vals)))

      df_native <- data.frame(time = fit$H0_df$time, hazard = fit$H0_df$hazard)
      df_native <- df_native[is.finite(df_native$time) & is.finite(df_native$hazard), , drop = FALSE]
      df_native <- df_native[order(df_native$time), , drop = FALSE]
      if (nrow(df_native) >= 2) {
        H0 <- df_native$hazard
        S0 <- exp(-H0)
        base_rows[[length(base_rows) + 1]] <- tibble(
          target_cancer = tname,
          method = "discretekl",
          t = df_native$time,
          h0 = NA_real_,
          H0 = H0,
          S0 = S0
        )
      }
    } else {
      beta_rows[[length(beta_rows) + 1]] <- tibble(
        target_cancer = tname,
        method = "discretekl"
      ) %>% bind_cols(as_tibble(setNames(as.list(rep(NA_real_, length(target$feature_cols))),
                                         target$feature_cols)))
    }
  }

  beta_df <- bind_rows(beta_rows)
  base_df_out <- bind_rows(base_rows)

  write.csv(beta_df, file.path(results_dir, "discretekl_beta_estimates.csv"), row.names = FALSE)
  write.csv(base_df_out, file.path(results_dir, "discretekl_baseline_hazard.csv"), row.names = FALSE)
  message(sprintf("[DiscreteKL] wrote %s", file.path(results_dir, "discretekl_beta_estimates.csv")))
  message(sprintf("[DiscreteKL] wrote %s", file.path(results_dir, "discretekl_baseline_hazard.csv")))
}

collect_external_methods <- function(results_dir) {
  transcox_beta_path <- file.path(results_dir, "transcox_beta_estimates.csv")
  transcox_base_path <- file.path(results_dir, "transcox_baseline_hazard.csv")
  discretekl_beta_path <- file.path(results_dir, "discretekl_beta_estimates.csv")
  discretekl_base_path <- file.path(results_dir, "discretekl_baseline_hazard.csv")

  beta_list <- list(
    read_optional_csv(transcox_beta_path),
    read_optional_csv(discretekl_beta_path)
  )
  base_list <- list(
    read_optional_csv(transcox_base_path),
    read_optional_csv(discretekl_base_path)
  )
  beta_list <- beta_list[!vapply(beta_list, is.null, logical(1))]
  base_list <- base_list[!vapply(base_list, is.null, logical(1))]

  list(
    beta = if (length(beta_list) > 0) bind_rows(beta_list) else NULL,
    base = if (length(base_list) > 0) bind_rows(base_list) else NULL
  )
}

write_all_methods_outputs <- function() {
  full_beta_path <- file.path(results_dir, "fullsample_beta_estimates.csv")
  full_base_path <- file.path(results_dir, "fullsample_baseline_hazard.csv")
  cox_beta_path <- file.path(results_dir, "coxph_beta_estimates.csv")
  cox_base_path <- file.path(results_dir, "coxph_baseline_hazard.csv")

  full_beta <- read_optional_csv(full_beta_path)
  full_base <- read_optional_csv(full_base_path)
  cox_beta <- read_optional_csv(cox_beta_path)
  cox_base <- read_optional_csv(cox_base_path)
  external <- collect_external_methods(results_dir)

  beta_list <- list(full_beta, cox_beta, external$beta)
  beta_list <- beta_list[!vapply(beta_list, is.null, logical(1))]
  base_list <- list(full_base, cox_base, external$base)
  base_list <- base_list[!vapply(base_list, is.null, logical(1))]

  if (length(beta_list) > 0) {
    all_beta <- bind_rows(beta_list)
    write.csv(all_beta, file.path(results_dir, "all_methods_beta_estimates.csv"), row.names = FALSE)
  }
  if (length(base_list) > 0) {
    all_base <- bind_rows(base_list)
    write.csv(all_base, file.path(results_dir, "all_methods_baseline_hazard.csv"), row.names = FALSE)
  }
}
# c_map <- c(
#   cancer2 = 0.5,
#   cancer12 = 0.5,
#   cancer13 = 0.3,
#   cancer15 = 0.8
# )
# c_default <- 1
c_map <- cv_tbl$c_mult
names(c_map) <- as.character(all_sheets)
main <- function() {
  cancers <- lapply(all_sheets, function(sh) {
    df <- read_excel(data_file, sheet = sh)
    df <- cap_cancer_df(df)
    na.omit(df)
  })
  names(cancers) <- all_sheets
  proc <- lapply(cancers, prep_df)
  
  selected_df <- read.csv(selected_file, stringsAsFactors = FALSE)
  
  beta_rows <- list()
  gamma_rows <- list()
  base_rows <- list()
  
  for (tname in target_sheets) {
    target <- proc[[tname]]
    target_df <- target$df
    lam <- lambda_map[[tname]]
    if (is.null(lam)) {
      lam <- list(lambda_zeta = lambda_zeta_default, lambda_eta = lambda_eta_default)
    }
    c_mult <- if (tname %in% names(c_map)) c_map[[tname]] else c_default
    lambda_zeta_t <- lam$lambda_zeta
    lambda_eta_t <- lam$lambda_eta
    stopifnot(is.numeric(lambda_zeta_t), length(lambda_zeta_t) == 1, is.finite(lambda_zeta_t))
    stopifnot(is.numeric(lambda_eta_t), length(lambda_eta_t) == 1, is.finite(lambda_eta_t))
    message(sprintf("[LAM] target=%s lambda_zeta=%g lambda_eta=%g",
                    tname, lambda_zeta_t, lambda_eta_t))
    
    sel_row <- selected_df[selected_df$target_cancer == tname, , drop = FALSE]
    sel_sources <- character(0)
    if (nrow(sel_row) > 0 &&
        !is.na(sel_row$n_selected[1]) && sel_row$n_selected[1] > 0) {
      sel_sources <- parse_sources(sel_row$selected_sources[1])
      sel_sources <- sel_sources[sel_sources != "" & !is.na(sel_sources)]
      sel_sources <- sel_sources[sel_sources %in% names(proc) & sel_sources != tname]
    }
    
    source_df_pool <- NULL
    if (length(sel_sources) > 0) {
      source_df_pool <- bind_rows(lapply(sel_sources, function(sn) proc[[sn]]$df))
    }
    
    apply_list <- list()
    if (!is.null(source_df_pool)) {
      apply_list$src <- source_df_pool
    }
    pre <- preprocess_features(target_df, apply_list, target$feature_cols, target$time_col, target$event_col)
    target_df <- pre$train_df
    if (!is.null(source_df_pool)) {
      source_df_pool <- pre$apply_df_list$src
    }
    
    fit_target <- fit_algorithm1_transfer_only(target_df, source_df = NULL,
                                               feature_cols = target$feature_cols,
                                               time_col = target$time_col,
                                               event_col = target$event_col,
                                               lambda_zeta = lambda_zeta_t,
                                               lambda_eta = lambda_eta_t,
                                               c_mult = c_mult)
    fit_transfer <- fit_algorithm1_transfer_only(target_df, source_df_pool,
                                                 feature_cols = target$feature_cols,
                                                 time_col = target$time_col,
                                                 event_col = target$event_col,
                                                 lambda_zeta = lambda_zeta_t,
                                                 lambda_eta = lambda_eta_t,
                                                 c_mult = c_mult)
    fit_full <- fit_algorithm1(target_df, source_df = source_df_pool,
                               feature_cols = target$feature_cols,
                               time_col = target$time_col,
                               event_col = target$event_col,
                               lambda_zeta = lambda_zeta_t,
                               lambda_eta = lambda_eta_t,
                               c_mult = c_mult, use_hessian_lasso_update = FALSE)
    
    fit_list <- list(
      target_only = fit_target,
      combined_transfer = fit_transfer,
      full_method = fit_full
    )
    
    for (mname in names(fit_list)) {
      fit <- fit_list[[mname]]
      conv_ok <- is.list(fit) && !is.null(fit$conv_transfer) &&
        fit$conv_transfer == 1L &&
        (!is.null(fit$beta_hat) && !is.null(fit$gamma_hat))
      if (mname == "full_method") {
        conv_ok <- conv_ok && !is.null(fit$conv_debias) && fit$conv_debias == 1L
      }
      if (!conv_ok) {
        beta_rows[[length(beta_rows) + 1]] <- tibble(
          target_cancer = tname,
          method = mname
        ) %>% bind_cols(as_tibble(setNames(as.list(rep(NA_real_, length(target$feature_cols))),
                                           target$feature_cols)))
        gamma_rows[[length(gamma_rows) + 1]] <- tibble(
          target_cancer = tname,
          method = mname,
          basis_term = NA_character_,
          gamma_value = NA_real_
        )
        base_rows[[length(base_rows) + 1]] <- tibble(
          target_cancer = tname,
          method = mname,
          t = NA_real_,
          h0 = NA_real_,
          H0 = NA_real_,
          S0 = NA_real_
        )
        next
      }
      
      beta_vals <- as.numeric(fit$beta_hat)
      names(beta_vals) <- target$feature_cols
      beta_rows[[length(beta_rows) + 1]] <- tibble(
        target_cancer = tname,
        method = mname
      ) %>% bind_cols(as_tibble(as.list(beta_vals)))
      
      gamma_vals <- as.numeric(fit$gamma_hat)
      gamma_rows[[length(gamma_rows) + 1]] <- tibble(
        target_cancer = tname,
        method = mname,
        basis_term = paste0("g", seq_along(gamma_vals)),
        gamma_value = gamma_vals
      )
      
      base_df <- compute_baseline(fit)
      if (!is.null(base_df) && all(c("t", "H0") %in% names(base_df))) {
        base_df <- base_df[is.finite(base_df$t) & is.finite(base_df$H0), , drop = FALSE]
        base_df <- base_df[order(base_df$t), , drop = FALSE]
        if (nrow(base_df) >= 2) {
          base_rows[[length(base_rows) + 1]] <- base_df %>%
            mutate(target_cancer = tname, method = mname) %>%
            dplyr::select(target_cancer, method, t, h0, H0, S0)
        }
      }
    }
  }
  
  beta_df <- bind_rows(beta_rows)
  gamma_df <- bind_rows(gamma_rows)
  base_df <- bind_rows(base_rows)
  
  write.csv(beta_df, file.path(results_dir, "fullsample_beta_estimates.csv"), row.names = FALSE)
  write.csv(gamma_df, file.path(results_dir, "fullsample_gamma_estimates.csv"), row.names = FALSE)
  write.csv(base_df, file.path(results_dir, "fullsample_baseline_hazard.csv"), row.names = FALSE)
}

plot_baseline_hazards <- function() {
  suppressPackageStartupMessages({
    library(ggplot2)
  })
  base_path <- file.path(results_dir, "fullsample_baseline_hazard.csv")
  if (!file.exists(base_path)) return(invisible(NULL))
  base_df <- read.csv(base_path, stringsAsFactors = FALSE)
  base_df <- base_df[base_df$method %in% c("target_only", "combined_transfer", "full_method"), ]
  if (nrow(base_df) == 0) return(invisible(NULL))

  clean_curve <- function(df, method_name, target_name) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    sub <- df[df$target_cancer == target_name & df$method == method_name, , drop = FALSE]
    if (nrow(sub) == 0 || !all(c("t", "H0") %in% names(sub))) return(NULL)
    sub <- sub[is.finite(sub$t) & is.finite(sub$H0), , drop = FALSE]
    if (nrow(sub) < 2) return(NULL)
    sub <- sub[order(sub$t), , drop = FALSE]
    sub$time <- sub$t
    sub$cumhaz <- sub$H0
    sub$method <- method_name
    sub[, c("time", "cumhaz", "method"), drop = FALSE]
  }

  method_levels <- c("target_only", "full_method")
  color_map <- c(
    #combined_transfer = "#F8766D",
    target_only = "#B79F00",
    full_method = "#B3B3B3"
  )
  label_map <- c(
    #combined_transfer = "Combined",
    target_only = "TargetOnly",
    full_method = "Cox-SieveTL"
  )
  for (cc in unique(base_df$target_cancer)) {
    curves <- list(
      clean_curve(base_df, "full_method", cc),
      # clean_curve(base_df, "combined_transfer", cc),
      clean_curve(base_df, "target_only", cc)
    )
    df_plot <- bind_rows(curves)
    if (nrow(df_plot) == 0) next
    df_plot$method <- factor(df_plot$method, levels = method_levels)
    y_ref <- df_plot$cumhaz[is.finite(df_plot$cumhaz)]
    ylim_vals <- if (length(y_ref) >= 2) range(y_ref) else NULL
    p <- ggplot(df_plot, aes(x = time, y = cumhaz, color = method, group = method)) +
      geom_line(linewidth = 1, na.rm = TRUE) +
      scale_color_manual(values = color_map, breaks = method_levels, labels = label_map, drop = FALSE) +
      labs(title = paste0("Cumulative baseline hazard H0(t): ", cc),
           x = "Time", y = "H0(t)", color = "method") +
      theme_bw()
    if (!is.null(ylim_vals)) {
      p <- p + coord_cartesian(ylim = ylim_vals)
    }
    print(p)
  }
}


compute_coxph_outputs <- function() {
  suppressPackageStartupMessages({
    library(survival)
    library(readxl)
    library(dplyr)
    library(tibble)
  })
  
  cancers <- lapply(all_sheets, function(sh) {
    df <- read_excel(data_file, sheet = sh)
    df <- cap_cancer_df(df)
    df %>% drop_na()
  })
  names(cancers) <- all_sheets
  proc <- lapply(cancers, prep_df)
  
  selected_df <- read.csv(selected_file, stringsAsFactors = FALSE)
  
  beta_rows <- list()
  base_rows <- list()
  
  safe_coxph_refit_drop_na_coef <- function(dat, feature_cols, time_col, event_col) {
    if (is.null(dat) || nrow(dat) < 3) {
      return(list(
        fit_full = NULL,
        fit_final = NULL,
        dropped_vars = character(0),
        beta_full = setNames(rep(NA_real_, length(feature_cols)), feature_cols),
        beta_out  = setNames(rep(NA_real_, length(feature_cols)), feature_cols)
      ))
    }
    
    fml_full <- as.formula(
      paste0(
        "survival::Surv(", time_col, ",", event_col, ") ~ ",
        paste(feature_cols, collapse = " + ")
      )
    )
    
    fit_full <- tryCatch(
      survival::coxph(fml_full, data = dat, x = TRUE),
      error = function(e) NULL
    )
    
    beta_template <- setNames(rep(NA_real_, length(feature_cols)), feature_cols)
    
    if (is.null(fit_full)) {
      return(list(
        fit_full = NULL,
        fit_final = NULL,
        dropped_vars = character(0),
        beta_full = beta_template,
        beta_out  = beta_template
      ))
    }
    
    coef_full <- stats::coef(fit_full)
    coef_named <- beta_template
    if (!is.null(coef_full)) {
      nm_intersect <- intersect(names(coef_full), feature_cols)
      coef_named[nm_intersect] <- coef_full[nm_intersect]
    }
    
    na_vars <- names(coef_named)[is.na(coef_named)]
    
    if (length(na_vars) == 0) {
      return(list(
        fit_full = fit_full,
        fit_final = fit_full,
        dropped_vars = character(0),
        beta_full = coef_named,
        beta_out  = coef_named
      ))
    }
    
    keep_vars <- setdiff(feature_cols, na_vars)
    
    if (length(keep_vars) == 0) {
      return(list(
        fit_full = fit_full,
        fit_final = NULL,
        dropped_vars = na_vars,
        beta_full = coef_named,
        beta_out  = beta_template
      ))
    }
    
    fml_refit <- as.formula(
      paste0(
        "survival::Surv(", time_col, ",", event_col, ") ~ ",
        paste(keep_vars, collapse = " + ")
      )
    )
    
    fit_refit <- tryCatch(
      survival::coxph(fml_refit, data = dat, x = TRUE),
      error = function(e) NULL
    )
    
    beta_out <- beta_template
    if (!is.null(fit_refit)) {
      coef_refit <- stats::coef(fit_refit)
      nm_intersect <- intersect(names(coef_refit), feature_cols)
      beta_out[nm_intersect] <- coef_refit[nm_intersect]
    }
    
    list(
      fit_full = fit_full,
      fit_final = fit_refit,
      dropped_vars = na_vars,
      beta_full = coef_named,
      beta_out  = beta_out
    )
  }
  
  for (tname in target_sheets) {
    target <- proc[[tname]]
    target_df <- target$df
    
    sel_row <- selected_df[selected_df$target_cancer == tname, , drop = FALSE]
    sel_sources <- character(0)
    if (nrow(sel_row) > 0 &&
        !is.na(sel_row$n_selected[1]) && sel_row$n_selected[1] > 0) {
      sel_sources <- parse_sources(sel_row$selected_sources[1])
      sel_sources <- sel_sources[sel_sources != "" & !is.na(sel_sources)]
      sel_sources <- sel_sources[sel_sources %in% names(proc) & sel_sources != tname]
    }
    
    source_df_pool <- NULL
    if (length(sel_sources) > 0) {
      source_df_pool <- bind_rows(lapply(sel_sources, function(sn) proc[[sn]]$df))
    }
    
    apply_list <- list()
    if (!is.null(source_df_pool)) {
      apply_list$src <- source_df_pool
    }
    
    pre <- preprocess_features(
      target_df,
      apply_list,
      target$feature_cols,
      target$time_col,
      target$event_col
    )
    
    target_df <- pre$train_df
    if (!is.null(source_df_pool)) {
      source_df_pool <- pre$apply_df_list$src
    }
    
    # -------------------------
    # target-only CoxPH
    # -------------------------
    target_res <- safe_coxph_refit_drop_na_coef(
      dat = target_df,
      feature_cols = target$feature_cols,
      time_col = target$time_col,
      event_col = target$event_col
    )
    
    message(sprintf(
      "[coxph_target] target=%s dropped_vars=%s",
      tname,
      ifelse(length(target_res$dropped_vars) == 0,
             "<none>",
             paste(target_res$dropped_vars, collapse = ", "))
    ))
    
    beta_rows[[length(beta_rows) + 1]] <- tibble(
      target_cancer = tname,
      method = "coxph_target"
    ) %>% bind_cols(as_tibble(as.list(target_res$beta_out)))
    
    if (!is.null(target_res$fit_final)) {
      bh <- tryCatch(
        survival::basehaz(target_res$fit_final, centered = FALSE),
        error = function(e) NULL
      )
      
      if (!is.null(bh)) {
        df_native <- data.frame(time = bh$time, hazard = bh$hazard)
        df_native <- df_native[is.finite(df_native$time) & is.finite(df_native$hazard), , drop = FALSE]
        df_native <- df_native[order(df_native$time), , drop = FALSE]
        if (nrow(df_native) >= 2) {
          H0 <- df_native$hazard
          S0 <- exp(-H0)
          base_rows[[length(base_rows) + 1]] <- tibble(
            target_cancer = tname,
            method = "coxph_target",
            t = df_native$time,
            h0 = NA_real_,
            H0 = H0,
            S0 = S0
          )
        }
      }
    }
    
    # -------------------------
    # combined target + selected sources CoxPH
    # -------------------------
    combined_df <- if (is.null(source_df_pool) || nrow(source_df_pool) == 0) {
      target_df
    } else {
      bind_rows(target_df, source_df_pool)
    }
    
    combined_res <- safe_coxph_refit_drop_na_coef(
      dat = combined_df,
      feature_cols = target$feature_cols,
      time_col = target$time_col,
      event_col = target$event_col
    )
    
    message(sprintf(
      "[coxph_combined] target=%s dropped_vars=%s",
      tname,
      ifelse(length(combined_res$dropped_vars) == 0,
             "<none>",
             paste(combined_res$dropped_vars, collapse = ", "))
    ))
    
    beta_rows[[length(beta_rows) + 1]] <- tibble(
      target_cancer = tname,
      method = "coxph_combined"
    ) %>% bind_cols(as_tibble(as.list(combined_res$beta_out)))
    
    if (!is.null(combined_res$fit_final)) {
      bh <- tryCatch(
        survival::basehaz(combined_res$fit_final, centered = FALSE),
        error = function(e) NULL
      )
      
      if (!is.null(bh)) {
        df_native <- data.frame(time = bh$time, hazard = bh$hazard)
        df_native <- df_native[is.finite(df_native$time) & is.finite(df_native$hazard), , drop = FALSE]
        df_native <- df_native[order(df_native$time), , drop = FALSE]
        if (nrow(df_native) >= 2) {
          H0 <- df_native$hazard
          S0 <- exp(-H0)
          base_rows[[length(base_rows) + 1]] <- tibble(
            target_cancer = tname,
            method = "coxph_combined",
            t = df_native$time,
            h0 = NA_real_,
            H0 = H0,
            S0 = S0
          )
        }
      }
    }
  }
  
  beta_df <- bind_rows(beta_rows)
  base_df_out <- bind_rows(base_rows)
  
  write.csv(beta_df, file.path(results_dir, "coxph_beta_estimates.csv"), row.names = FALSE)
  write.csv(base_df_out, file.path(results_dir, "coxph_baseline_hazard.csv"), row.names = FALSE)
}

plot_baseline_hazards_with_coxph <- function(exclude_targets = character(0)) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(cowplot)
    library(readxl)
    library(dplyr)
  })
  
  # ---- compute n_target / n_source for each target cancer ----
  cancers <- lapply(all_sheets, function(sh) {
    df <- read_excel(data_file, sheet = sh)
    df <- cap_cancer_df(df)
    na.omit(df)
  })
  names(cancers) <- all_sheets
  proc <- lapply(cancers, prep_df)
  
  selected_df <- read.csv(selected_file, stringsAsFactors = FALSE)
  
  size_map <- list()
  for (tname in target_sheets) {
    target <- proc[[tname]]
    target_df_raw <- target$df
    
    sel_row <- selected_df[selected_df$target_cancer == tname, , drop = FALSE]
    sel_sources <- character(0)
    if (nrow(sel_row) > 0 &&
        !is.na(sel_row$n_selected[1]) && sel_row$n_selected[1] > 0) {
      sel_sources <- parse_sources(sel_row$selected_sources[1])
      sel_sources <- sel_sources[sel_sources != "" & !is.na(sel_sources)]
      sel_sources <- sel_sources[sel_sources %in% names(proc) & sel_sources != tname]
    }
    
    source_df_pool <- NULL
    if (length(sel_sources) > 0) {
      source_df_pool <- bind_rows(lapply(sel_sources, function(sn) proc[[sn]]$df))
    }
    
    apply_list <- list()
    if (!is.null(source_df_pool)) {
      apply_list$src <- source_df_pool
    }
    
    pre <- preprocess_features(
      target_df_raw,
      apply_list,
      target$feature_cols,
      target$time_col,
      target$event_col
    )
    
    target_df <- pre$train_df
    if (!is.null(source_df_pool)) {
      source_df_pool <- pre$apply_df_list$src
    }
    
    size_map[[tname]] <- list(
      n_target = nrow(target_df),
      n_source = if (is.null(source_df_pool)) 0L else nrow(source_df_pool)
    )
  }
  
  base_path <- file.path(results_dir, "fullsample_baseline_hazard.csv")
  cox_path  <- file.path(results_dir, "coxph_baseline_hazard.csv")
  if (!file.exists(base_path) || !file.exists(cox_path)) return(invisible(NULL))
  
  base_df <- read.csv(base_path, stringsAsFactors = FALSE)
  cox_df  <- read.csv(cox_path, stringsAsFactors = FALSE)
  external <- collect_external_methods(results_dir)
  all_df <- bind_rows(base_df, cox_df, external$base)
  all_df <- all_df[!all_df$method %in% c('target_only', 'CoxPH Source'),]
  message(sprintf("[PLOT] all_df nrow=%s", nrow(all_df)))
  if ("method" %in% names(all_df)) {
    message("[PLOT] all_df method table:")
    print(table(all_df$method, useNA = "ifany"))
  } else {
    message("[PLOT] all_df method table: <missing column>")
  }
  message("[PLOT] all_df t summary:")
  print(summary(all_df$t))
  message("[PLOT] all_df H0 summary:")
  print(summary(all_df$H0))
  
  all_df$method <- tolower(all_df$method)
  message(sprintf("[PLOT] unique methods before filter=%s",
                  paste(unique(all_df$method), collapse = ",")))
  
  all_df <- all_df[all_df$method %in% c("target_only", "combined_transfer", "full_method",
                                        "coxph_target", "coxph_combined", "transcox",
                                        "discretekl"), ]
  if (nrow(all_df) == 0) return(invisible(NULL))
  
  clean_curve <- function(df, method_name, target_name) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    sub <- df[df$target_cancer == target_name & df$method == method_name, , drop = FALSE]
    if (nrow(sub) == 0 || !all(c("t", "H0") %in% names(sub))) return(NULL)
    sub <- sub[is.finite(sub$t) & is.finite(sub$H0), , drop = FALSE]
    if (nrow(sub) < 2) return(NULL)
    sub <- sub[order(sub$t), , drop = FALSE]
    sub$time <- sub$t
    sub$cumhaz <- sub$H0
    sub$method <- method_name
    sub[, c("time", "cumhaz", "method"), drop = FALSE]
  }
  
  method_levels <- c("full_method",
                     "coxph_target", "coxph_combined", "transcox",
                     "discretekl")
  color_map <- c(
    target_only = "#B79F00",
    transcox = "#619CFF",
    discretekl = "#C77CFF",
    full_method = "#F8766D",
    coxph_target = "#000000",
    coxph_combined = "#009E73"
  )
  label_map <- c(
    #combined_transfer = "SieveTL Combined",
    #target_only = "TargetOnly",
    transcox = "TransCox",
    discretekl = "DiscreteKL",
    full_method = "Cox-SieveTL",
    coxph_target = "CoxPH Target",
    coxph_combined = "CoxPH Combined"
  )
  linetype_map <- c(
    # target_only = "solid",
    #combined_transfer = "longdash",
    full_method = "solid",
    transcox = "dotdash",
    discretekl = "dashed",
    coxph_target = "solid",
    coxph_combined = "dotted"
  )
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
  
  plot_list <- list()
  legend_df <- NULL
  ordered_targets <- target_sheets[target_sheets %in% unique(all_df$target_cancer)]
  ordered_targets <- ordered_targets[!ordered_targets %in% c(exclude_targets)]
  
  for (cc in ordered_targets) {
    curves <- list(
      clean_curve(all_df, "full_method", cc),
      #clean_curve(all_df, "combined_transfer", cc),
      #clean_curve(all_df, "target_only", cc),
      clean_curve(all_df, "coxph_target", cc),
      clean_curve(all_df, "coxph_combined", cc),
      clean_curve(all_df, "transcox", cc),
      clean_curve(all_df, "discretekl", cc)
    )
    df_plot <- bind_rows(curves)
    if (nrow(df_plot) == 0) next
    if (is.null(legend_df)) {
      legend_df <- df_plot
    }
    
    df_plot$method <- factor(df_plot$method, levels = method_levels)
    y_ref <- df_plot$cumhaz[is.finite(df_plot$cumhaz) & df_plot$method != "transcox"]
    ylim_vals <- if (length(y_ref) >= 2) range(y_ref) else NULL
    
    n_target_cc <- size_map[[cc]]$n_target
    n_source_cc <- size_map[[cc]]$n_source
    
    note_txt <- paste0("nT=", n_target_cc, "\n", "nS=", n_source_cc)
    
    p <- ggplot(df_plot, aes(x = time, y = cumhaz, color = method, linetype = method, group = method)) +
      geom_line(linewidth = 1, na.rm = TRUE) +
      scale_color_manual(values = color_map, breaks = method_levels, labels = label_map, drop = FALSE) +
      scale_linetype_manual(values = linetype_map, breaks = method_levels, labels = label_map, drop = FALSE) +
      labs(title = title_map[[cc]], color = NULL, linetype = NULL) +
      annotate(
        "label",
        x = -Inf, y = Inf,
        label = note_txt,
        hjust = -0.1, vjust = 1.1,
        size = 3,
        linewidth = 0.25
      ) +
      theme_bw() +
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none"
      )
    
    if (!is.null(ylim_vals)) {
      if (abs(max(ylim_vals) - quantile(df_plot$cumhaz, 0.95, na.rm = TRUE)) > 10) {
        p <- p + coord_cartesian(ylim = c(0, 2))
      } else {
        p <- p + coord_cartesian(ylim = ylim_vals)
      }
    }
    
    plot_list[[length(plot_list) + 1]] <- p
  }
  
  if (length(plot_list) == 0) return(invisible(NULL))
  if (is.null(legend_df)) return(invisible(NULL))
  
  legend_dummy <- data.frame(
    time = rep(1:2, times = length(method_levels)),
    cumhaz = rep(c(0, 1), times = length(method_levels)),
    method = factor(rep(method_levels, each = 2), levels = method_levels)
  )
  
  legend_plot <- ggplot(legend_dummy, aes(x = time, y = cumhaz, color = method, linetype = method)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = color_map, breaks = method_levels, labels = label_map, drop = FALSE) +
    scale_linetype_manual(values = linetype_map, breaks = method_levels, labels = label_map, drop = FALSE) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.key.width = grid::unit(2.2, "cm"),
      legend.title = element_blank()
    )
  
  shared_legend <- cowplot::get_legend(legend_plot)
  grid_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = 5, nrow = 3)
  combined_plot <- cowplot::plot_grid(grid_plot, shared_legend, ncol = 1, rel_heights = c(1, 0.08))
  
  out_path <- file.path(results_dir, "cumbasehaz_3x3_sharedlegend_new.eps")
  ggsave(out_path, combined_plot, width = 14, height = 14, units = "in", device = "eps")
  
  return(combined_plot)
}

if (sys.nframe() == 0) {
  main()
  compute_coxph_outputs()
  compute_transcox_outputs()
  compute_discretekl_outputs()
  write_all_methods_outputs()
  plot_baseline_hazards_with_coxph(exclude_targets = c("cancer11"))
}

