#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(timeROC)
  library(splines)
  library(pracma)
})

setwd('/Users/yuxisong/Library/CloudStorage/Box-Box/Cox-SieveTL/simulation/results/new_icb_response_result')

data_file <- "./codes/extracted_cancer_data_by_type.xlsx"
results_dir <- "./results_allTMB_new_hessian_lasso"

if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

log_conv <- function(...) message(sprintf(...))

detect_columns <- function(df) {
  cn <- names(df)
  time_col <- cn[cn %in% c("OS_Months", "time", "Time", "time_months")]
  event_col <- cn[cn %in% c("OS_Event", "event", "delta", "status")]
  if (length(time_col) == 0 || length(event_col) == 0) {
    stop("Could not detect time/event columns.")
  }
  list(time = time_col[1], event = event_col[1])
}

prep_df <- function(df) {
  cols <- detect_columns(df)
  df <- df %>% dplyr::select(-cancer) %>% drop_na()
  
  time <- as.numeric(df[[cols$time]])
  event <- as.numeric(df[[cols$event]])
  feat <- df %>% dplyr::select(-all_of(c(cols$time, cols$event)))
  
  X <- model.matrix(~ . - 1, data = feat)
  covar_cols <- colnames(X)
  
  df_clean <- data.frame(time = time, event = event, X, check.names = FALSE)
  
  list(
    df = df_clean,
    feature_cols = covar_cols,
    time_col = "time",
    event_col = "event"
  )
}

cap_cancer_df <- function(df) {
  if ("Age" %in% names(df)) df$Age <- pmin(df$Age, 85)
  if ("NLR" %in% names(df)) df$NLR <- pmin(df$NLR, 25)
  if ("TMB" %in% names(df)) df$TMB <- pmin(df$TMB, 50)
  df
}

preprocess_features <- function(train_df, apply_df_list, feature_cols, time_col, event_col) {
  train_X <- train_df[, feature_cols, drop = FALSE]
  num_cols <- feature_cols[vapply(train_X, is.numeric, logical(1))]
  
  if (length(num_cols) == 0) {
    return(list(train_df = train_df, apply_df_list = apply_df_list,
                center = numeric(0), scale = numeric(0)))
  }
  
  is_binary <- vapply(train_X[, num_cols, drop = FALSE], function(x) {
    ux <- unique(na.omit(x))
    length(ux) > 0 && all(ux %in% c(0, 1))
  }, logical(1))
  
  scale_cols <- num_cols[!is_binary]
  
  if (length(scale_cols) == 0) {
    return(list(train_df = train_df, apply_df_list = apply_df_list,
                center = numeric(0), scale = numeric(0)))
  }
  
  means <- vapply(train_X[, scale_cols, drop = FALSE], mean, numeric(1), na.rm = TRUE)
  sds <- vapply(train_X[, scale_cols, drop = FALSE], sd, numeric(1), na.rm = TRUE)
  
  adj_means <- ifelse(is.finite(means), means, 0)
  adj_sds <- ifelse(is.finite(sds) & sds > 0, sds, 1)
  
  scale_df <- function(df) {
    df2 <- df
    for (col in scale_cols) {
      df2[[col]] <- (df2[[col]] - adj_means[col]) / adj_sds[col]
    }
    df2
  }
  
  list(
    train_df = scale_df(train_df),
    apply_df_list = lapply(apply_df_list, scale_df),
    center = adj_means,
    scale = adj_sds
  )
}

build_basis_quantile <- function(times,
                                 degree = 3,
                                 n_internal = NULL,
                                 c_mult = 1) {
  times <- times[is.finite(times)]
  if (length(times) == 0) stop("No finite times for spline basis.")
  
  if (is.null(n_internal)) {
    n_internal <- max(1, floor(c_mult * length(unique(times))^(1 / 3)))
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

.smooth_abs <- function(u, eps) sqrt(u * u + eps * eps) - eps
.smooth_abs_grad <- function(u, eps) u / sqrt(u * u + eps * eps)

soft_thresh <- function(z, a) {
  sign(z) * pmax(abs(z) - a, 0)
}

transfer_loss_fast <- function(theta, p, L,
                               all.x, all.stime, all.type,
                               G_mat, time_grid, g_time_mat) {
  gamma <- theta[1:p]
  beta  <- theta[(p + 1):(L + p)]
  
  xb <- as.vector(all.x %*% beta)
  g_gamma <- if (p == 1) as.vector(G_mat * gamma) else as.vector(G_mat %*% gamma)
  
  term1 <- all.type * (xb + g_gamma)
  
  exp_gamma_g <- as.vector(exp(g_time_mat %*% gamma))
  idx_end <- findInterval(all.stime, time_grid)
  cumI <- pracma::cumtrapz(time_grid, exp_gamma_g)
  
  int_vec <- numeric(length(all.stime))
  valid <- idx_end >= 2
  int_vec[valid] <- cumI[idx_end[valid]]
  
  term2 <- exp(xb) * int_vec
  
  -mean(term1 - term2, na.rm = TRUE)
}

transfer_grad_fast <- function(theta, p, L,
                               all.x, all.stime, all.type,
                               G_mat, time_grid, g_time_mat) {
  gamma <- theta[1:p]
  beta  <- theta[(p + 1):(L + p)]
  
  xb <- as.vector(all.x %*% beta)
  exp_xb <- exp(xb)
  
  eg_time <- as.vector(g_time_mat %*% gamma)
  exp_eg <- exp(eg_time)
  
  idx_end <- findInterval(all.stime, time_grid)
  cumI <- pracma::cumtrapz(time_grid, exp_eg)
  
  int_vec <- numeric(length(all.stime))
  valid <- idx_end >= 2
  int_vec[valid] <- cumI[idx_end[valid]]
  
  term_mat <- matrix(0, nrow = length(all.stime), ncol = p)
  for (j in seq_len(p)) {
    cumIj <- pracma::cumtrapz(time_grid, g_time_mat[, j] * exp_eg)
    term_mat[valid, j] <- cumIj[idx_end[valid]]
  }
  
  grad_gamma <- -colMeans(G_mat * all.type, na.rm = TRUE) +
    colMeans(exp_xb * term_mat, na.rm = TRUE)
  
  grad_beta <- -colMeans(all.x * all.type, na.rm = TRUE) +
    colMeans(all.x * (exp_xb * int_vec), na.rm = TRUE)
  
  c(grad_gamma, grad_beta)
}

debias_loss_fast <- function(delta, lambda_zeta, lambda_eta, L, p,
                             target, G_mat_target, time_grid_target,
                             g_time_mat_target, gamma_hat_A, beta_hat_A,
                             eps = 1e-4) {
  zeta <- delta[1:p]
  eta <- delta[(p + 1):(L + p)]
  
  gamma <- gamma_hat_A + zeta
  beta <- beta_hat_A + eta
  
  xb <- as.vector(target$x %*% beta)
  g_gamma <- if (p == 1) {
    as.vector(G_mat_target * gamma)
  } else {
    as.vector(G_mat_target %*% gamma)
  }
  
  term1 <- target$type * (xb + g_gamma)
  
  eg_time <- as.vector(g_time_mat_target %*% gamma)
  exp_eg <- exp(eg_time)
  
  idx_end <- findInterval(target$stime, time_grid_target)
  cumI <- pracma::cumtrapz(time_grid_target, exp_eg)
  
  int_vec <- numeric(length(target$stime))
  valid <- idx_end >= 2
  int_vec[valid] <- cumI[idx_end[valid]]
  
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
  
  idx_end <- findInterval(target$stime, time_grid_target)
  cumI <- pracma::cumtrapz(time_grid_target, exp_eg)
  
  int_vec <- numeric(length(target$stime))
  valid <- idx_end >= 2
  int_vec[valid] <- cumI[idx_end[valid]]
  
  term_mat <- matrix(0, nrow = length(target$stime), ncol = p)
  for (j in seq_len(p)) {
    cumIj <- pracma::cumtrapz(time_grid_target, g_time_mat_target[, j] * exp_eg)
    term_mat[valid, j] <- cumIj[idx_end[valid]]
  }
  
  grad_zeta <- -colMeans(G_mat_target * target$type, na.rm = TRUE) +
    colMeans(exp_xb * term_mat, na.rm = TRUE)
  
  grad_eta <- -colMeans(target$x * target$type, na.rm = TRUE) +
    colMeans(target$x * (exp_xb * int_vec), na.rm = TRUE)
  
  grad_zeta <- grad_zeta + lambda_zeta * .smooth_abs_grad(zeta, eps)
  grad_eta <- grad_eta + lambda_eta * .smooth_abs_grad(eta, eps)
  
  c(grad_zeta, grad_eta)
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
    
    exp_eta <- exp(as.vector(G_i %*% gamma_hat) + sum(x_i * beta_hat))
    
    for (a in seq_along(tg)) {
      z_a <- c(G_i[a, ], x_i)
      H <- H + tcrossprod(z_a) * exp_eta[a]
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
    H = H,
    H_reg = H_reg,
    grad = grad
  )
}

fit_algorithm1 <- function(target_train_df, source_df = NULL,
                           feature_cols,
                           time_col = "time",
                           event_col = "event",
                           lambda_zeta = 0.05,
                           lambda_eta = 0.05,
                           c_mult = 1,
                           maxit = 5000,
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
  
  basis <- build_basis_quantile(all_time, degree = 3, c_mult = c_mult)
  g_eval <- basis$g_eval
  
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
  
  G_mat <- g_eval(all.stime)
  p <- ncol(G_mat)
  L <- ncol(all.x)
  
  time_grid <- seq(min(all_time), max(all_time), length.out = 500)
  g_time_mat <- g_eval(time_grid)
  
  theta0 <- rep(0, L + p)
  
  loss_wrap <- function(th) {
    transfer_loss_fast(
      theta = th,
      p = p,
      L = L,
      all.x = all.x,
      all.stime = all.stime,
      all.type = all.type,
      G_mat = G_mat,
      time_grid = time_grid,
      g_time_mat = g_time_mat
    )
  }
  
  grad_wrap <- function(th) {
    transfer_grad_fast(
      theta = th,
      p = p,
      L = L,
      all.x = all.x,
      all.stime = all.stime,
      all.type = all.type,
      G_mat = G_mat,
      time_grid = time_grid,
      g_time_mat = g_time_mat
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
  beta_hat_A <- theta_hat_A[(p + 1):(L + p)]
  
  G_mat_target <- g_eval(target_list$stime)
  time_grid_target <- seq(min(target_list$stime), max(target_list$stime), length.out = 500)
  g_time_mat_target <- g_eval(time_grid_target)
  
  delta0 <- rep(0, L + p)
  
  deb_loss <- function(d) {
    debias_loss_fast(
      delta = d,
      lambda_zeta = lambda_zeta,
      lambda_eta = lambda_eta,
      L = L,
      p = p,
      target = target_list,
      G_mat_target = G_mat_target,
      time_grid_target = time_grid_target,
      g_time_mat_target = g_time_mat_target,
      gamma_hat_A = gamma_hat_A,
      beta_hat_A = beta_hat_A
    )
  }
  
  deb_grad <- function(d) {
    debias_grad_fast(
      delta = d,
      lambda_zeta = lambda_zeta,
      lambda_eta = lambda_eta,
      L = L,
      p = p,
      target = target_list,
      G_mat_target = G_mat_target,
      time_grid_target = time_grid_target,
      g_time_mat_target = g_time_mat_target,
      gamma_hat_A = gamma_hat_A,
      beta_hat_A = beta_hat_A
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
      g_eval = g_eval,
      ridge_iter = ridge_iter
    )
    
    delta_hat <- sparse_fit$delta_hat
    gamma_hat <- sparse_fit$gamma_hat
    beta_hat <- sparse_fit$beta_hat
    H_lasso <- sparse_fit$H
    grad_lasso <- sparse_fit$grad
  } else {
    delta_hat <- delta_hat_optim
    gamma_hat <- gamma_hat_A + delta_hat[1:p]
    beta_hat <- beta_hat_A + delta_hat[(p + 1):(L + p)]
    H_lasso <- NULL
    grad_lasso <- NULL
  }
  
  zeta_hat <- delta_hat[1:p]
  eta_hat <- delta_hat[(p + 1):(L + p)]
  
  list(
    g_eval = g_eval,
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
    use_hessian_lasso_update = use_hessian_lasso_update,
    basis = basis
  )
}

neg_loglik_sievetl <- function(fit, val_df, feature_cols, time_col, event_col) {
  X_val <- as.matrix(val_df[, feature_cols, drop = FALSE])
  T_val <- as.numeric(val_df[[time_col]])
  delta_val <- as.numeric(val_df[[event_col]])
  
  beta_hat <- fit$beta_hat
  gamma_hat <- fit$gamma_hat
  g_eval <- fit$g_eval
  
  G_val <- g_eval(T_val)
  xb <- as.vector(X_val %*% beta_hat)
  eta <- xb + as.vector(G_val %*% gamma_hat)
  
  time_grid <- fit$time_grid_fit
  g_time <- fit$g_time_mat_fit
  
  exp_g <- as.vector(exp(g_time %*% gamma_hat))
  idx_end <- findInterval(T_val, time_grid)
  cumI <- pracma::cumtrapz(time_grid, exp_g)
  
  int_vec <- numeric(length(T_val))
  valid <- idx_end >= 2
  int_vec[valid] <- cumI[idx_end[valid]]
  
  -mean(delta_val * eta - exp(xb) * int_vec, na.rm = TRUE)
}

make_folds <- function(event, k = 3, seed = 123) {
  set.seed(seed)
  idx1 <- sample(which(event == 1))
  idx0 <- sample(which(event == 0))
  
  folds <- vector("list", k)
  
  for (i in seq_along(idx1)) {
    folds[[ (i - 1) %% k + 1 ]] <- c(folds[[ (i - 1) %% k + 1 ]], idx1[i])
  }
  
  for (i in seq_along(idx0)) {
    folds[[ (i - 1) %% k + 1 ]] <- c(folds[[ (i - 1) %% k + 1 ]], idx0[i])
  }
  
  lapply(folds, sort)
}

cv_screen_target <- function(target_list, sources, target_name,
                             k = 3, seed = 123, C0 = 1.0,
                             lambda_zeta = 0.05,
                             lambda_eta = 0.05,
                             c_mult = 1,
                             use_hessian_lasso_update = TRUE) {
  
  folds <- make_folds(target_list$df[[target_list$event_col]], k = k, seed = seed)
  
  L0_vec <- numeric(k)
  Ls_mat <- matrix(Inf, nrow = length(sources), ncol = k,
                   dimnames = list(names(sources), NULL))
  
  conv0_tr <- integer(k)
  conv0_db <- integer(k)
  
  convS_tr <- matrix(0L, nrow = length(sources), ncol = k,
                     dimnames = list(names(sources), NULL))
  convS_db <- matrix(0L, nrow = length(sources), ncol = k,
                     dimnames = list(names(sources), NULL))
  
  source_ok <- setNames(rep(TRUE, length(sources)), names(sources))
  
  for (fold_idx in seq_len(k)) {
    val_idx <- folds[[fold_idx]]
    train_idx <- setdiff(seq_len(nrow(target_list$df)), val_idx)
    
    train_df <- target_list$df[train_idx, , drop = FALSE]
    val_df <- target_list$df[val_idx, , drop = FALSE]
    
    apply_list <- c(list(val = val_df), lapply(sources, function(s) s$df))
    names(apply_list) <- c("val", names(sources))
    
    pre <- preprocess_features(
      train_df = train_df,
      apply_df_list = apply_list,
      feature_cols = target_list$feature_cols,
      time_col = target_list$time_col,
      event_col = target_list$event_col
    )
    
    train_df <- pre$train_df
    val_df <- pre$apply_df_list$val
    sources_scaled <- pre$apply_df_list[names(sources)]
    
    fit0 <- fit_algorithm1(
      target_train_df = train_df,
      source_df = NULL,
      feature_cols = target_list$feature_cols,
      time_col = target_list$time_col,
      event_col = target_list$event_col,
      lambda_zeta = lambda_zeta,
      lambda_eta = lambda_eta,
      c_mult = c_mult,
      use_hessian_lasso_update = use_hessian_lasso_update
    )
    
    conv0_tr[fold_idx] <- fit0$conv_transfer
    conv0_db[fold_idx] <- fit0$conv_debias
    
    log_conv(
      "[target=%s fold=%d] target-only conv_transfer=%d conv_debias=%d nonzero_delta=%s",
      target_name, fold_idx,
      fit0$conv_transfer, fit0$conv_debias,
      ifelse(is.null(fit0$n_nonzero_delta), NA, fit0$n_nonzero_delta)
    )
    
    if (!(fit0$conv_transfer == 1L && fit0$conv_debias == 1L)) {
      log_conv("[target=%s] SKIP target: target-only nonconverged at fold=%d",
               target_name, fold_idx)
      return(NULL)
    }
    
    L0_vec[fold_idx] <- neg_loglik_sievetl(
      fit0, val_df,
      target_list$feature_cols,
      target_list$time_col,
      target_list$event_col
    )
    
    for (sname in names(sources)) {
      if (!source_ok[[sname]]) {
        Ls_mat[sname, fold_idx] <- Inf
        next
      }
      
      s_df <- sources_scaled[[sname]]
      
      fit_s <- fit_algorithm1(
        target_train_df = train_df,
        source_df = s_df,
        feature_cols = target_list$feature_cols,
        time_col = target_list$time_col,
        event_col = target_list$event_col,
        lambda_zeta = lambda_zeta,
        lambda_eta = lambda_eta,
        c_mult = c_mult,
        use_hessian_lasso_update = use_hessian_lasso_update
      )
      
      convS_tr[sname, fold_idx] <- fit_s$conv_transfer
      convS_db[sname, fold_idx] <- fit_s$conv_debias
      
      log_conv(
        "[target=%s fold=%d source=%s] conv_transfer=%d conv_debias=%d nonzero_delta=%s",
        target_name, fold_idx, sname,
        fit_s$conv_transfer, fit_s$conv_debias,
        ifelse(is.null(fit_s$n_nonzero_delta), NA, fit_s$n_nonzero_delta)
      )
      
      if (!(fit_s$conv_transfer == 1L && fit_s$conv_debias == 1L)) {
        source_ok[[sname]] <- FALSE
        Ls_mat[sname, fold_idx] <- Inf
        log_conv("[target=%s source=%s] DROP source: nonconverged at fold=%d",
                 target_name, sname, fold_idx)
        next
      }
      
      Ls_mat[sname, fold_idx] <- neg_loglik_sievetl(
        fit_s, val_df,
        target_list$feature_cols,
        target_list$time_col,
        target_list$event_col
      )
    }
  }
  
  L0_bar <- mean(L0_vec, na.rm = TRUE)
  sigma_hat <- sqrt(0.5 * sum((L0_vec - L0_bar)^2, na.rm = TRUE))
  threshold <- C0 * max(sigma_hat, 0.01)
  
  Ls_bar <- rowMeans(Ls_mat)
  delta_L <- Ls_bar - L0_bar
  pass <- (delta_L <= threshold) & source_ok
  
  n_sources <- length(Ls_bar)
  
  conv_target_tr_df <- as.data.frame(
    as.list(setNames(conv0_tr, paste0("conv_target_tr_fold", seq_len(k)))),
    check.names = FALSE
  )
  conv_target_db_df <- as.data.frame(
    as.list(setNames(conv0_db, paste0("conv_target_db_fold", seq_len(k)))),
    check.names = FALSE
  )
  
  conv_target_tr_df <- conv_target_tr_df[rep(1, n_sources), , drop = FALSE]
  conv_target_db_df <- conv_target_db_df[rep(1, n_sources), , drop = FALSE]
  
  conv_source_tr_df <- as.data.frame(convS_tr, check.names = FALSE)
  colnames(conv_source_tr_df) <- paste0("conv_source_tr_fold", seq_len(k))
  conv_source_tr_df <- conv_source_tr_df[names(Ls_bar), , drop = FALSE]
  
  conv_source_db_df <- as.data.frame(convS_db, check.names = FALSE)
  colnames(conv_source_db_df) <- paste0("conv_source_db_fold", seq_len(k))
  conv_source_db_df <- conv_source_db_df[names(Ls_bar), , drop = FALSE]
  
  tibble(
    source_cancer = names(Ls_bar),
    L_source = as.numeric(Ls_bar),
    L_targetonly = L0_bar,
    diff = delta_L,
    pass = pass,
    threshold = threshold
  ) %>%
    bind_cols(conv_target_tr_df, conv_target_db_df,
              conv_source_tr_df, conv_source_db_df) %>%
    arrange(diff)
}

load_cancer <- function(sheet) {
  df <- read_excel(data_file, sheet = sheet)
  n_raw <- nrow(df)
  
  df <- cap_cancer_df(df)
  n_after_cap <- nrow(df)
  
  df <- df %>% drop_na()
  n_after_cap_dropna <- nrow(df)
  
  message(
    sheet,
    " rows raw: ", n_raw,
    "; rows after cap: ", n_after_cap,
    "; rows after cap+drop_na: ", n_after_cap_dropna
  )
  
  df
}

main <- function() {
  sheets <- paste0("cancer", c(1:10, 12:16))
  
  cancers <- lapply(sheets, load_cancer)
  names(cancers) <- sheets
  
  surv_rate <- readxl::read_excel(
    '/Users/yuxisong/Library/CloudStorage/Box-Box/Cox-SieveTL/real data/5yearsurvival.xlsx'
  )
  
  surv_rate <- surv_rate[c('Cancer Type', 'Index', '5-year relative  survival')]
  surv_rate$`5-year relative  survival` <- round(
    as.numeric(surv_rate$`5-year relative  survival`),
    4
  )
  
  proc <- lapply(cancers, prep_df)
  
  n_min <- 0
  surv_bound <- 0.2
  
  target_sheets <- names(proc)[sapply(proc, function(x) nrow(x$df) >= n_min)]
  
  all_results <- list()
  sel_summary <- list()
  
  for (tname in target_sheets) {
    target <- proc[[tname]]
    target_sr <- surv_rate$`5-year relative  survival`[surv_rate$Index == tname]
    
    sources <- proc[
      names(proc) != tname &
        sapply(names(proc), function(nm) {
          cond_n <- nrow(proc[[nm]]$df) >= n_min
          
          sr_src <- surv_rate$`5-year relative  survival`[surv_rate$Index == nm]
          sr_tgt <- target_sr
          
          cond_sr <- if (length(sr_src) == 0L || is.na(sr_src) || is.na(sr_tgt)) {
            TRUE
          } else {
            abs(sr_src - sr_tgt) <= surv_bound
          }
          
          cond_n & cond_sr
        })
    ]
    
    res <- cv_screen_target(
      target_list = target,
      sources = sources,
      target_name = tname,
      k = 3,
      seed = 2025,
      C0 = 1.0,
      lambda_zeta = 0.05,
      lambda_eta = 0.05,
      c_mult = 1,
      use_hessian_lasso_update = TRUE
    )
    
    if (is.null(res) || nrow(res) == 0) next
    
    res <- res %>%
      mutate(target_cancer = tname) %>%
      select(target_cancer, everything())
    
    all_results[[tname]] <- res
    
    write.csv(
      res,
      file.path(results_dir, paste0("source_selection_", tname, ".csv")),
      row.names = FALSE
    )
    
    selected <- res %>% filter(pass) %>% pull(source_cancer)
    
    sel_summary[[tname]] <- tibble(
      target_cancer = tname,
      selected_sources = paste(selected, collapse = ", "),
      n_selected = length(selected)
    )
  }
  
  all_res_df <- bind_rows(all_results)
  write.csv(
    all_res_df,
    file.path(results_dir, "source_selection_all.csv"),
    row.names = FALSE
  )
  
  sel_summary_df <- bind_rows(sel_summary)
  write.csv(
    sel_summary_df,
    file.path(results_dir, "selected_sources_summary.csv"),
    row.names = FALSE
  )
}

if (sys.nframe() == 0) {
  main()
}