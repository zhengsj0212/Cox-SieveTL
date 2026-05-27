#!/usr/bin/env Rscript
# Use BIC to select lambda_zeta and lambda_eta for SieveTL,
# while still calculating CV prediction metrics.
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(timeROC)
  library(survival)
  library(survAUC)
  library(pracma)
  library(stringr)
})
setwd('/Users/yuxisong/Library/CloudStorage/Box-Box/Cox-SieveTL/simulation/results/new_icb_response_result')
source("./codes/cv_source_selection.R")
results_dir <- "./results_allTMB_new_hessian_lasso"
data_file <- "./codes/extracted_cancer_data_by_type.xlsx"
selected_file <- file.path(results_dir, "selected_sources_summary.csv")

# lambda_zeta_grid <- c(1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2)
# lambda_eta_grid  <- c(0, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2)


lambda_grid <- list(
  lambda_zeta = c(1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2),
  lambda_eta  = c(1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2)
)
lambda_grid_special <- lambda_grid

# lambda_grid_special <- list(
#   cancer2 = list(   # Breast
#     lambda_zeta = c(1e-3, 2e-3, 5e-3, 1e-2),
#     lambda_eta  = c(0, 1e-6, 1e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2)
#   ),
#   cancer10 = list(  # Mesothelioma
#     lambda_zeta = c(1e-3, 2e-3, 5e-3, 1e-2),
#     lambda_eta  = c(0, 1e-6, 1e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2)
#   ),
#   cancer13 = list(  # Pancreatic
#     lambda_zeta = c(1e-3, 2e-3, 5e-3, 1e-2),
#     lambda_eta  = c(0, 1e-6, 1e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2)
#   ),
#   cancer18 = list(  # CNS
#     lambda_zeta = c(1e-3, 2e-3, 5e-3, 1e-2),
#     lambda_eta  = c(0, 1e-6, 1e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2)
#   )
# )

get_lambda_grid <- function(tname) {
  if (tname %in% names(lambda_grid_special)) {
    lambda_grid_special[[tname]]
  } else {
    lambda_grid
  }
}

lambda_zeta_default <- 0.005
lambda_eta_default <- 0

parse_sources <- function(s) {
  if (is.null(s) || length(s) == 0 || is.na(s) || trimws(s) == "") return(character(0))
  trimws(strsplit(s, ",")[[1]])
}

sd_safe <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) NA_real_ else sd(x)
}

mean_safe <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 1) NA_real_ else mean(x)
}

count_nonzero <- function(x, tol = 1e-8) {
  sum(abs(x) > tol, na.rm = TRUE)
}

# ------------------------------------------------------------
# Metric evaluation
# ------------------------------------------------------------

eval_fold_metrics <- function(fit, train_df, val_df, feature_cols, time_col, event_col) {
  if (is.null(fit$beta_hat) || is.null(fit$gamma_hat) || is.null(fit$g_eval)) {
    return(list(
      C_harrell = NA_real_,
      C_uno = NA_real_,
      auc_vec = NA_real_,
      mean_auc = NA_real_,
      brier = NA_real_
    ))
  }
  
  lp <- as.vector(as.matrix(val_df[, feature_cols, drop = FALSE]) %*% fit$beta_hat)
  
  val_time  <- as.numeric(val_df[[time_col]])
  val_event <- as.numeric(val_df[[event_col]])
  train_time  <- as.numeric(train_df[[time_col]])
  train_event <- as.numeric(train_df[[event_col]])
  
  # -------------------
  # Harrell C
  # -------------------
  C_harrell <- tryCatch({
    survival::concordance(Surv(val_time, val_event) ~ lp, reverse = TRUE)$concordance
  }, error = function(e) {
    message("[HarrellC ERROR] ", conditionMessage(e))
    NA_real_
  })
  
  # -------------------
  # Uno C
  # -------------------
  C_uno <- tryCatch({
    if (sum(train_event == 1, na.rm = TRUE) < 2 || sum(val_event == 1, na.rm = TRUE) < 2) {
      return(NA_real_)
    }
    
    uno_time <- as.numeric(stats::quantile(train_time, 0.8, na.rm = TRUE))
    uno_time <- min(uno_time, max(val_time, na.rm = TRUE) - 1e-8)
    uno_time <- max(uno_time, min(val_time, na.rm = TRUE) + 1e-8)
    
    get_uno <- function(score) {
      as.numeric(survAUC::UnoC(
        survival::Surv(train_time, train_event),
        survival::Surv(val_time, val_event),
        lpnew = as.numeric(score),
        time  = uno_time
      ))
    }
    
    C1 <- tryCatch(get_uno(lp), error = function(e) NA_real_)
    C2 <- tryCatch(get_uno(-lp), error = function(e) NA_real_)
    if (is.finite(C1) && (is.na(C2) || C1 >= C2)) C1 else C2
  }, error = function(e) {
    message("[UnoC ERROR] ", conditionMessage(e))
    NA_real_
  })
  
  # -------------------
  # time-dependent AUC
  # -------------------
  auc_vec <- NA_real_
  mean_auc <- NA_real_
  
  if (sum(val_event == 1, na.rm = TRUE) >= 2 &&
      sum(val_event == 0, na.rm = TRUE) >= 2 &&
      length(unique(lp[is.finite(lp)])) >= 2) {
    
    event_times <- sort(unique(val_time[val_event == 1 & is.finite(val_time)]))
    
    if (length(event_times) >= 3) {
      # use interior event times only
      probs <- c(0.2, 0.4, 0.6, 0.8)
      times_auc <- unique(as.numeric(stats::quantile(event_times, probs = probs, na.rm = TRUE, type = 7)))
      times_auc <- times_auc[is.finite(times_auc)]
      
      tmin <- min(val_time, na.rm = TRUE) + 1e-8
      tmax <- max(val_time, na.rm = TRUE) - 1e-8
      times_auc <- times_auc[times_auc > tmin & times_auc < tmax]
      times_auc <- sort(unique(times_auc))
      
      if (length(times_auc) >= 1) {
        auc_vec <- tryCatch({
          roc_obj <- timeROC::timeROC(
            T = val_time,
            delta = val_event,
            marker = lp,
            cause = 1,
            times = times_auc,
            iid = FALSE
          )
          as.numeric(roc_obj$AUC)
        }, error = function(e) {
          message("[timeROC ERROR] ", conditionMessage(e))
          NA_real_
        })
      }
    }
  }
  
  mean_auc <- if (all(is.na(auc_vec))) NA_real_ else mean(auc_vec, na.rm = TRUE)
  
  # -------------------
  # Brier
  # -------------------
  t_grid <- sort(unique(train_df[[time_col]]))
  if (length(t_grid) < 2 || all(is.na(auc_vec))) {
    brier <- NA_real_
  } else {
    brier <- tryCatch({
      times_brier <- times_auc
      g_time <- fit$g_eval(t_grid)
      exp_g <- exp(g_time %*% fit$gamma_hat)
      H0 <- pracma::cumtrapz(t_grid, exp_g)
      H0_t <- approx(t_grid, H0, xout = times_brier, rule = 2)$y
      S_mat <- exp(-exp(lp) %o% H0_t)
      
      cens_fit <- survival::survfit(Surv(train_df[[time_col]], 1 - train_df[[event_col]]) ~ 1)
      G_t  <- summary(cens_fit, times = times_brier, extend = TRUE)$surv
      G_Ti <- summary(cens_fit, times = pmax(val_time - 1e-8, 0), extend = TRUE)$surv
      
      epsG <- 1e-6
      G_t  <- pmax(G_t, epsG)
      G_Ti <- pmax(G_Ti, epsG)
      
      y_mat <- outer(val_time, times_brier, function(ti, tt) as.numeric(ti > tt))
      w_mat <- matrix(0, nrow = length(val_time), ncol = length(times_brier))
      for (j in seq_along(times_brier)) {
        t_j <- times_brier[j]
        w_mat[val_time <= t_j & val_event == 1, j] <- 1 / G_Ti[val_time <= t_j & val_event == 1]
        w_mat[val_time > t_j, j] <- 1 / G_t[j]
      }
      brier_vec <- colMeans(w_mat * (y_mat - S_mat)^2, na.rm = TRUE)
      mean(brier_vec, na.rm = TRUE)
    }, error = function(e) {
      message("[Brier ERROR] ", conditionMessage(e))
      NA_real_
    })
  }
  
  list(
    C_harrell = C_harrell,
    C_uno = C_uno,
    auc_vec = auc_vec,
    mean_auc = mean_auc,
    brier = brier
  )
}

# ------------------------------------------------------------
# BIC for fitted SieveTL
# ------------------------------------------------------------

compute_fit_bic <- function(fit, train_df, feature_cols, time_col, event_col, tol = 1e-6) {
  if (!is.list(fit) ||
      is.null(fit$conv_transfer) || fit$conv_transfer != 1L ||
      is.null(fit$conv_debias)   || fit$conv_debias   != 1L ||
      is.null(fit$beta_hat) || is.null(fit$gamma_hat) || is.null(fit$g_eval)) {
    return(Inf)
  }
  
  X <- as.matrix(train_df[, feature_cols, drop = FALSE])
  time  <- as.numeric(train_df[[time_col]])
  event <- as.numeric(train_df[[event_col]])
  
  keep <- is.finite(time) & is.finite(event) & complete.cases(X)
  X <- X[keep, , drop = FALSE]
  time <- time[keep]
  event <- event[keep]
  
  n <- length(time)
  if (n < 5) return(Inf)
  
  lp <- as.vector(X %*% fit$beta_hat)
  
  t_grid <- sort(unique(time))
  if (length(t_grid) < 2) return(Inf)
  
  g_time <- fit$g_eval(t_grid)
  if (is.null(dim(g_time))) return(Inf)
  
  h0_grid <- as.vector(exp(g_time %*% fit$gamma_hat))
  H0_grid <- pracma::cumtrapz(t_grid, h0_grid)
  
  h0_i <- approx(t_grid, h0_grid, xout = time, rule = 2)$y
  H0_i <- approx(t_grid, H0_grid, xout = time, rule = 2)$y
  
  if (any(!is.finite(h0_i)) || any(!is.finite(H0_i))) return(Inf)
  
  loglik <- sum(event * (log(pmax(h0_i, 1e-12)) + lp) - exp(lp) * H0_i)
  if (!is.finite(loglik)) return(Inf)
  
  if (!is.null(fit$zeta_hat) && !is.null(fit$eta_hat)) {
    k_eff <- count_nonzero(fit$zeta_hat, tol = tol) +
      count_nonzero(fit$eta_hat, tol = tol)
  } else if (!is.null(fit$delta_hat)) {
    k_eff <- count_nonzero(fit$delta_hat, tol = tol)
  } else {
    k_eff <- count_nonzero(fit$eta_hat, tol = tol) +
      count_nonzero(fit$zeta_hat, tol = tol)
  }
  
  bic <- -2 * loglik + log(n) * k_eff
  bic
}

build_basis_quantile <- function(times,
                                 degree = 3,
                                 n_internal = NULL,
                                 c_mult = 1,
                                 boundary = range(times, na.rm = TRUE)) {
  times <- times[is.finite(times)]
  if (length(times) == 0) stop("No finite times for spline basis.")
  
  if (is.null(n_internal)) {
    n_internal <- max(1, floor(c_mult * length(unique(times))^(1/3)))
  }
  
  probs <- seq(0, 1, length.out = n_internal + 2)
  knots_all <- as.numeric(quantile(times, probs = probs, na.rm = TRUE, type = 7))
  boundary <- c(knots_all[1], knots_all[length(knots_all)])
  knots <- sort(unique(knots_all[-c(1, length(knots_all))]))
  
  g_eval <- function(tvec) {
    splines::bs(
      tvec,
      degree = degree,
      knots = knots,
      Boundary.knots = boundary,
      intercept = TRUE
    )
  }
  
  list(
    knots = knots,
    boundary = boundary,
    degree = degree,
    g_eval = g_eval,
    n_internal = n_internal,
    c_mult = c_mult
  )
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
    x_i <- as.numeric(data$x[i, ])
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

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------
c_map <- c(
  cancer2 = 0.5,
  cancer12 = 0.5,
  cancer13 = 0.3,
  cancer15 = 0.8
)
#   c(
#   cancer13 = 0.5
# )

c_default <- 1

main <- function(
    K = 5,
    seed = 2026,
    run_targets = NULL
) {
  sheets <- paste0("cancer", 1:18)
  
  cancers <- lapply(sheets, function(sh) {
    df <- read_excel(data_file, sheet = sh)
    df <- cap_cancer_df(df)
    df %>% drop_na()
  })
  names(cancers) <- sheets
  
  proc <- lapply(cancers, prep_df)
  for (nm in names(proc)) proc[[nm]]$name <- nm
  
  selected_df <- read.csv(selected_file, stringsAsFactors = FALSE)
  
  if (!is.null(run_targets)) {
    selected_df <- selected_df %>%
      filter(target_cancer %in% run_targets)
  }
  
  all_summaries <- list()
  tuned_rows <- list()
  grid_rows <- list()
  all_fold_rows <- list()
  
  for (i in seq_len(nrow(selected_df))) {
    
    tname <- selected_df$target_cancer[i]
    if (!tname %in% names(proc)) next
    if(tname == 'cancer13') {
      c_grid = c(0.3, 0.5, 0.8, 1.0)
    } else {
      c_grid = c(0.5, 0.8, 1.0)
    }
    message("\n==============================")
    message("Target cancer: ", tname)
    message("==============================")
    
    target <- proc[[tname]]
    
    lam_grid <- get_lambda_grid(tname)
    lambda_zeta_grid_base <- lam_grid$lambda_zeta
    lambda_eta_grid_base  <- lam_grid$lambda_eta
    
    sel_sources <- character(0)
    if (!is.na(selected_df$n_selected[i]) && selected_df$n_selected[i] > 0) {
      sel_sources <- parse_sources(selected_df$selected_sources[i])
      sel_sources <- sel_sources[sel_sources %in% names(proc) & sel_sources != tname]
    }
    
    # ============================================================
    # Step 1: Full-sample preprocessing
    # ============================================================
    
    target_full_df <- target$df
    
    source_full_df <- NULL
    if (length(sel_sources) > 0) {
      source_full_df <- bind_rows(lapply(sel_sources, function(sn) proc[[sn]]$df))
    }
    
    apply_list_full <- list()
    if (!is.null(source_full_df)) {
      apply_list_full$src <- source_full_df
    }
    
    pre_full <- preprocess_features(
      target_full_df,
      apply_list_full,
      target$feature_cols,
      target$time_col,
      target$event_col
    )
    
    target_full_df <- pre_full$train_df
    if (!is.null(source_full_df)) {
      source_full_df <- pre_full$apply_df_list$src
    }
    
    all_time <- if (is.null(source_full_df)) {
      target_full_df[[target$time_col]]
    } else {
      c(
        target_full_df[[target$time_col]],
        source_full_df[[target$time_col]]
      )
    }
    
    L <- length(target$feature_cols)
    
    # ============================================================
    # Step 2: Joint BIC tuning over c, lambda_zeta, lambda_eta
    # ============================================================
    
    grid_metrics <- list()
    
    for (c_mult in c_grid) {
      
      n_internal <- max(
        1,
        floor(c_mult * length(unique(all_time))^(1/3))
      )
      
      pn <- n_internal + 4
      
      fix_le <- sqrt(log(pn + L) / nrow(target_full_df))
      fix_lz <- sqrt(log(pn + L) / nrow(target_full_df))
      
      lambda_eta_grid  <- fix_le * lambda_eta_grid_base
      lambda_zeta_grid <- fix_lz * lambda_zeta_grid_base
      
      for (lz in lambda_zeta_grid) {
        for (le in lambda_eta_grid) {
          
          fit_full <- fit_algorithm1(
            target_full_df,
            source_df = source_full_df,
            feature_cols = target$feature_cols,
            time_col = target$time_col,
            event_col = target$event_col,
            lambda_zeta = lz,
            lambda_eta = le,
            c_mult = c_mult,
            use_hessian_lasso_update = TRUE
          )
          
          conv_full <- as.integer(
            is.list(fit_full) &&
              !is.null(fit_full$conv_transfer) &&
              !is.null(fit_full$conv_debias) &&
              fit_full$conv_transfer == 1L &&
              fit_full$conv_debias == 1L
          )
          
          bic_val <- if (conv_full == 1L) {
            compute_fit_bic(
              fit = fit_full,
              train_df = target_full_df,
              feature_cols = target$feature_cols,
              time_col = target$time_col,
              event_col = target$event_col
            )
          } else {
            Inf
          }
          
          message(sprintf(
            "[FULL-BIC] target=%s c=%g lambda_zeta=%g lambda_eta=%g conv=%d BIC=%s",
            tname, c_mult, lz, le, conv_full,
            ifelse(is.finite(bic_val), sprintf("%.3f", bic_val), "Inf")
          ))
          
          grid_metrics[[length(grid_metrics) + 1]] <- tibble(
            target_cancer = tname,
            c_mult = c_mult,
            n_internal = n_internal,
            pn = pn,
            lambda_zeta = lz,
            lambda_eta = le,
            BIC = bic_val,
            conv_full = conv_full
          )
        }
      }
    }
    
    grid_df <- bind_rows(grid_metrics) %>%
      arrange(BIC, desc(conv_full), c_mult, lambda_zeta, lambda_eta)
    
    grid_rows[[length(grid_rows) + 1]] <- grid_df
    
    best_row <- grid_df[1, , drop = FALSE]
    
    message(sprintf(
      "[SELECTED] target=%s best_c=%g best_lambda_zeta=%g best_lambda_eta=%g BIC=%s",
      tname,
      best_row$c_mult,
      best_row$lambda_zeta,
      best_row$lambda_eta,
      ifelse(is.finite(best_row$BIC), sprintf("%.3f", best_row$BIC), "Inf")
    ))
    
    fit_full <- fit_algorithm1(
      target_full_df,
      source_df = source_full_df,
      feature_cols = target$feature_cols,
      time_col = target$time_col,
      event_col = target$event_col,
      lambda_zeta = best_row$lambda_zeta,
      lambda_eta = best_row$lambda_eta,
      c_mult = best_row$c_mult,
      use_hessian_lasso_update = TRUE
    )
    
    tuned_rows[[length(tuned_rows) + 1]] <- tibble(
      target_cancer = tname,
      best_c_mult = best_row$c_mult,
      best_n_internal = best_row$n_internal,
      best_pn = best_row$pn,
      best_lambda_zeta = best_row$lambda_zeta,
      best_lambda_eta = best_row$lambda_eta,
      BIC = best_row$BIC,
      conv_full = best_row$conv_full
    )
    
    # ============================================================
    # Step 3: CV metrics using selected c/lambda only
    # ============================================================
    
    folds <- make_folds(target$df[[target$event_col]], k = K, seed = seed)
    
    C_harrell_vec <- numeric(0)
    C_uno_vec <- numeric(0)
    auc_fold_vec <- numeric(0)
    brier_vec <- numeric(0)
    conv_vec <- integer(0)
    
    for (fold_idx in seq_len(K)) {
      
      val_idx <- folds[[fold_idx]]
      train_idx <- setdiff(seq_len(nrow(target$df)), val_idx)
      
      train_df <- target$df[train_idx, , drop = FALSE]
      val_df   <- target$df[val_idx, , drop = FALSE]
      
      source_df_pool <- NULL
      if (length(sel_sources) > 0) {
        source_df_pool <- bind_rows(lapply(sel_sources, function(sn) proc[[sn]]$df))
      }
      
      apply_list <- list(val = val_df)
      if (!is.null(source_df_pool)) {
        apply_list$src <- source_df_pool
      }
      
      pre <- preprocess_features(
        train_df,
        apply_list,
        target$feature_cols,
        target$time_col,
        target$event_col
      )
      
      train_df <- pre$train_df
      val_df <- pre$apply_df_list$val
      
      if (!is.null(source_df_pool)) {
        source_df_pool <- pre$apply_df_list$src
      }
      
      fit <- fit_algorithm1(
        train_df,
        source_df = source_df_pool,
        feature_cols = target$feature_cols,
        time_col = target$time_col,
        event_col = target$event_col,
        lambda_zeta = best_row$lambda_zeta,
        lambda_eta = best_row$lambda_eta,
        c_mult = best_row$c_mult,
        use_hessian_lasso_update = TRUE
      )
      
      conv <- as.integer(
        is.list(fit) &&
          !is.null(fit$conv_transfer) &&
          !is.null(fit$conv_debias) &&
          fit$conv_transfer == 1L &&
          fit$conv_debias == 1L
      )
      
      message(sprintf(
        "[CV] target=%s c=%g lambda_zeta=%g lambda_eta=%g fold=%d conv=%d",
        tname,
        best_row$c_mult,
        best_row$lambda_zeta,
        best_row$lambda_eta,
        fold_idx,
        conv
      ))
      
      if (conv == 1L) {
        metrics <- eval_fold_metrics(
          fit,
          train_df,
          val_df,
          target$feature_cols,
          target$time_col,
          target$event_col
        )
      } else {
        metrics <- list(
          C_harrell = NA_real_,
          C_uno = NA_real_,
          auc_vec = NA_real_,
          mean_auc = NA_real_,
          brier = NA_real_
        )
      }
      
      C_harrell_vec <- c(C_harrell_vec, metrics$C_harrell)
      C_uno_vec <- c(C_uno_vec, metrics$C_uno)
      auc_fold_vec <- c(auc_fold_vec, metrics$mean_auc)
      brier_vec <- c(brier_vec, metrics$brier)
      conv_vec <- c(conv_vec, conv)
      
      all_fold_rows[[length(all_fold_rows) + 1]] <- tibble(
        target_cancer = tname,
        method = "full_method",
        fold = fold_idx,
        c_mult = best_row$c_mult,
        lambda_zeta = best_row$lambda_zeta,
        lambda_eta = best_row$lambda_eta,
        conv = conv,
        C_harrell = metrics$C_harrell,
        C_uno = metrics$C_uno,
        auc_fold = metrics$mean_auc,
        brier = metrics$brier
      )
    }
    
    summary_df <- tibble(
      target_cancer = tname,
      c_mult = best_row$c_mult,
      BIC = if (is.finite(best_row$BIC)) best_row$BIC else NA_real_,
      Harrell_C = mean_safe(C_harrell_vec),
      Uno_C = mean_safe(C_uno_vec),
      AUC = mean_safe(auc_fold_vec),
      Brier = mean_safe(brier_vec),
      ConvRate = if (length(conv_vec) > 0) mean(conv_vec) else NA_real_
    ) %>%
      mutate(
        BIC = sprintf("%.3f", BIC),
        Harrell_C = sprintf("%.3f", Harrell_C),
        Uno_C = sprintf("%.3f", Uno_C),
        AUC = sprintf("%.3f", AUC),
        Brier = sprintf("%.3f", Brier),
        ConvRate = sprintf("%.2f", ConvRate)
      )
    
    all_summaries[[tname]] <- summary_df
  }
  
  tuning_df <- bind_rows(tuned_rows)
  # tuning_df <- read.csv(file.path(results_dir, "c_lambda_tuning_by_cancer_bic.csv"))
  write.csv(
    tuning_df,
    file.path(results_dir, "c_lambda_tuning_by_cancer_bic.csv"),
    row.names = FALSE
  )
  
  grid_df_all <- bind_rows(grid_rows)
  write.csv(
    grid_df_all,
    file.path(results_dir, "c_lambda_tuning_grid_results_bic.csv"),
    row.names = FALSE
  )
  
  all_summary_df <- bind_rows(all_summaries)
  write.csv(
    all_summary_df,
    file.path(results_dir, "cv_metrics_summary_table_bic.csv"),
    row.names = FALSE
  )
  # all_summary_df <- read.csv(file.path(results_dir, "cv_metrics_summary_table_bic.csv"))
  summary_with_lambda <- all_summary_df %>%
    left_join(
      tuning_df %>%
        dplyr::select(
          target_cancer,
          best_c_mult,
          best_n_internal,
          best_pn,
          best_lambda_zeta,
          best_lambda_eta
        ),
      by = "target_cancer"
    )
  
  write.csv(
    summary_with_lambda,
    file.path(results_dir, "cv_metrics_summary_table_with_c_lambda_bic.csv"),
    row.names = FALSE
  )
  
  all_folds_df <- bind_rows(all_fold_rows)
  
  write.csv(
    all_folds_df,
    file.path(results_dir, "cv_metrics_allfolds_with_c_lambda_bic.csv"),
    row.names = FALSE
  )
  
  invisible(list(
    tuning = tuning_df,
    grid = grid_df_all,
    summary = summary_with_lambda,
    folds = all_folds_df
  ))
}

plot_km_target_source_combined <- function(
    out_file = file.path(results_dir, "km_target_source_combined.pdf"),
    ncol = 3
) {
  suppressPackageStartupMessages({
    library(readxl)
    library(dplyr)
    library(ggplot2)
    library(survival)
  })
  
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
  
  sheets <- paste0("cancer", c(1:16, 18))
  
  cancers <- lapply(sheets, function(sh) {
    df <- read_excel(data_file, sheet = sh)
    df <- cap_cancer_df(df)
    df %>% drop_na()
  })
  names(cancers) <- sheets
  
  proc <- lapply(cancers, prep_df)
  for (nm in names(proc)) proc[[nm]]$name <- nm
  
  selected_df <- read.csv(selected_file, stringsAsFactors = FALSE)
  
  make_km_df <- function(df, time_col, event_col, group_label, target_name) {
    if (is.null(df) || nrow(df) < 2) return(NULL)
    
    fit <- tryCatch(
      survfit(Surv(df[[time_col]], df[[event_col]]) ~ 1),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    
    tibble(
      time = c(0, fit$time),
      surv = c(1, fit$surv),
      group = group_label,
      target_cancer = target_name
    )
  }
  
  km_rows <- list()
  
  for (i in seq_len(nrow(selected_df))) {
    tname <- selected_df$target_cancer[i]
    if (!tname %in% names(proc)) next
    
    target <- proc[[tname]]
    target_df <- target$df
    
    sel_sources <- character(0)
    if (!is.na(selected_df$n_selected[i]) && selected_df$n_selected[i] > 0) {
      sel_sources <- parse_sources(selected_df$selected_sources[i])
      sel_sources <- sel_sources[sel_sources %in% names(proc) & sel_sources != tname]
    }
    
    # target curve
    df_target <- make_km_df(
      df = target_df,
      time_col = target$time_col,
      event_col = target$event_col,
      group_label = "Target",
      target_name = tname
    )
    km_rows <- append(km_rows, list(df_target))
    
    # each source separately
    if (length(sel_sources) > 0) {
      for (sn in sel_sources) {
        source_df <- proc[[sn]]$df
        
        source_label <- if (sn %in% names(title_map)) {
          paste0("Source: ", unname(title_map[sn]))
        } else {
          paste0("Source: ", sn)
        }
        
        df_source <- make_km_df(
          df = source_df,
          time_col = proc[[sn]]$time_col,
          event_col = proc[[sn]]$event_col,
          group_label = source_label,
          target_name = tname
        )
        
        km_rows <- append(km_rows, list(df_source))
      }
    }
  }
  
  km_df <- bind_rows(km_rows)
  
  if (nrow(km_df) == 0) {
    message("No KM data generated.")
    return(invisible(NULL))
  }
  
  km_df <- km_df %>%
    mutate(
      cancer_label = ifelse(
        target_cancer %in% names(title_map),
        unname(title_map[target_cancer]),
        target_cancer
      )
    )
  
  cancer_levels <- unname(title_map[names(title_map) %in% unique(km_df$target_cancer)])
  km_df$cancer_label <- factor(km_df$cancer_label, levels = cancer_levels)
  
  # keep Target first, then sources
  # keep Target first, then sources
  source_levels <- sort(unique(km_df$group[km_df$group != "Target"]))
  group_levels <- c("Target", source_levels)
  km_df$group <- factor(km_df$group, levels = group_levels)
  
  # generate enough colors
  n_groups <- length(group_levels)
  color_vals <- setNames(grDevices::hcl.colors(n_groups, palette = "Dark 3"), group_levels)
  color_vals["Target"] <- "black"
  
  # line types: target solid, sources recycled
  lty_pool <- c("dashed", "dotted", "dotdash", "longdash", "twodash")
  lty_vals <- setNames(rep(lty_pool, length.out = n_groups), group_levels)
  lty_vals["Target"] <- "solid"
  
  # remove bad rows before plotting
  km_df_plot <- km_df %>%
    filter(is.finite(time), is.finite(surv), !is.na(group), !is.na(cancer_label))
  
  p <- ggplot(km_df_plot, aes(x = time, y = surv, color = group, linetype = group)) +
    geom_step(linewidth = 1) +
    facet_wrap(~ cancer_label, scales = "free_x", ncol = ncol) +
    scale_color_manual(values = color_vals) +
    scale_linetype_manual(values = lty_vals) +
    guides(
      color = guide_legend(ncol = 8, byrow = TRUE),
      linetype = guide_legend(ncol = 8, byrow = TRUE)
    ) +
    labs(
      x = "Time",
      y = "Survival probability",
      color = NULL,
      linetype = NULL,
      title = NULL,
        #"Kaplan-Meier curves: target and selected sources"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 9)
    )
  
  print(p)
  
  if (!is.null(out_file)) {
    ggsave(out_file, p, width = 14, height = 14)
  }
  
  return(list(
    plot = p,
    data = km_df
  ))
}

if (sys.nframe() == 0) {
  main(run_targets = NULL)
  plot_km_target_source_combined()
}

