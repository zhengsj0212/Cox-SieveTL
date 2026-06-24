#!/usr/bin/env Rscript

# Column-specific pilot tuning for the transfer inverse-Hessian step.
# Paper-matching rule:
# 1. pooled initializer is unpenalized;
# 2. only the target correction penalty c_j is tuned by column;
# 3. selection criterion uses 0.5 * residual + (log n0 / n0) * ||delta_j||_0.

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  normalizePath(getwd())
}

script_dir <- get_script_dir()
source(file.path(script_dir, "run_PLC_CR_spline_hoffman2.R"))

hessian_c_grid <- c(0, 0.1, 0.5, 1, 2, 4, 8)
npilot <- as.integer(Sys.getenv("NPILOT", "20"))
pilot_seed_base <- as.integer(Sys.getenv("PILOT_SEED_BASE", "1"))
pilot_task_id <- as.integer(Sys.getenv("SGE_TASK_ID", "1"))
pilot_base_out <- path.expand(Sys.getenv("PILOT_OUT_BASE", unset = getwd()))

if (!is.finite(pilot_task_id) || pilot_task_id < 1L || pilot_task_id > npilot) {
  stop("Invalid pilot task id: ", pilot_task_id, " with NPILOT=", npilot)
}

pilot_seed <- pilot_seed_base + pilot_task_id
set.seed(pilot_seed)

pilot_design_dir <- default_pilot_design_dir(
  beta2 = beta2, kappa = kappa, c1 = c1,
  n0 = n_target, nA = n_source, L = L,
  shift = covariate_shift, base_dir = pilot_base_out
)
dir.create(pilot_design_dir, recursive = TRUE, showWarnings = FALSE)

build_basis_objects <- function(target, source, cspline) {
  pooled_stime <- c(target$stime, source[[1]]$stime)
  N <- length(unique(pooled_stime[c(target$type, source[[1]]$type) == 1]))
  J <- floor(cspline * N^(1/3))
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
  g_funcs <- lapply(seq_len(ncol(B0)), function(j) {
    function(t) {
      Bt <- splines::bs(t, degree = deg, knots = knots,
                        Boundary.knots = bkn, intercept = TRUE)
      as.numeric(Bt[, j])
    }
  })
  list(
    pooled_stime = pooled_stime,
    g_funcs = g_funcs
  )
}

D <- simulate_once_all(
  n_target = n_target, n_source = n_source,
  beta2_source = beta2, kappa_source = kappa,
  c1 = c1, covariate_shift = covariate_shift, L = L
)

tuning_file <- default_sievetl_tuning_path(
  beta2 = beta2, kappa = kappa, c1 = c1,
  n0 = n_target, nA = n_source, L = L,
  shift = covariate_shift, base_dir = pilot_base_out
)
tuning_row <- read_sievetl_tuning_file(tuning_file)

cspline_hat <- as.numeric(tuning_row$c_hat[1])
cz_hat <- as.numeric(tuning_row$c_zeta_hat[1])
ce_hat <- as.numeric(tuning_row$c_eta_hat[1])

result <- sievetl_approx(
  target = D$target, source = D$source,
  L = L,
  c_zeta = cz_hat,
  c_eta = ce_hat,
  n0_penalty = n_target,
  c = cspline_hat
)

basis_obj <- build_basis_objects(D$target, D$source, cspline_hat)
time_grid <- seq(min(D$target$stime, na.rm = TRUE),
                 max(D$target$stime, na.rm = TRUE), length.out = 500)
pooled_time_grid <- seq(min(basis_obj$pooled_stime, na.rm = TRUE),
                        max(basis_obj$pooled_stime, na.rm = TRUE), length.out = 500)
H_target <- compute_target_hessian_theta(
  D$target, result$gamma.hat, result$beta.hat,
  g_funcs = basis_obj$g_funcs, time_grid = time_grid
)
H_pooled <- compute_pooled_hessian_theta(
  D$target, D$source, result$gamma.hat, result$beta.hat,
  g_funcs = basis_obj$g_funcs, time_grid = pooled_time_grid
)

d_n <- ncol(H_target)
tie_break_rule <- paste(
  "smallest BIC_j;",
  "ties broken by smallest residual;",
  "then smaller c"
)

full_rows <- vector("list", d_n)
best_rows <- vector("list", d_n)

for (j in seq_len(d_n)) {
  e_j <- rep(0, d_n)
  e_j[j] <- 1
  omega_A_j <- fit_inverse_hessian_initializer_unpenalized(H_pooled, j)

  column_rows <- vector("list", length(hessian_c_grid))

  for (k in seq_along(hessian_c_grid)) {
    c_val <- hessian_c_grid[k]
    lambda_val <- c_val * sqrt(log(d_n) / n_target)

    delta_j <- fit_inverse_hessian_correction(H_target, omega_A_j, j, lambda_val, eps = 1e-4)
    omega_j <- omega_A_j + delta_j

    residual_j <- sum((H_target %*% omega_j - e_j)^2)
    l0_delta_j <- sum(abs(delta_j) > 1e-4)
    criterion_j <- 0.5 * residual_j +
      (log(n_target) / n_target) * l0_delta_j

    column_rows[[k]] <- data.frame(
      j = j,
      c = c_val,
      lambda = lambda_val,
      d_n = d_n,
      n_target = n_target,
      residual = residual_j,
      l0_delta = l0_delta_j,
      criterion = criterion_j,
      stringsAsFactors = FALSE
    )
  }

  column_table <- do.call(rbind, column_rows)
  column_table <- column_table[order(
    column_table$criterion,
    column_table$residual,
    column_table$c
  ), ]

  full_rows[[j]] <- column_table
  best_rows[[j]] <- data.frame(
    j = j,
    c_hat = cspline_hat,
    c_zeta_hat = cz_hat,
    c_eta_hat = ce_hat,
    lambda_zeta_hat = result$lambda_zeta,
    lambda_eta_hat = result$lambda_eta,
    d_hat = d_n,
    best_c = column_table$c[1],
    best_lambda_formula = "c * sqrt(log(d_n) / n_target)",
    best_lambda_value = column_table$lambda[1],
    best_criterion = column_table$criterion[1],
    best_residual = column_table$residual[1],
    best_l0_delta = column_table$l0_delta[1],
    tie_break_rule_used = tie_break_rule,
    stringsAsFactors = FALSE
  )
}

pilot_bic_tables_by_column <- do.call(rbind, full_rows)
pilot_selected_c_by_column <- do.call(rbind, best_rows)
selected_c_distribution <- pilot_selected_c_by_column |>
  dplyr::group_by(best_c) |>
  dplyr::summarise(
    freq_selected = dplyr::n(),
    prop_selected = dplyr::n() / nrow(pilot_selected_c_by_column),
    mean_selected_criterion = mean(best_criterion),
    mean_selected_residual = mean(best_residual),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(freq_selected), mean_selected_criterion, best_c)

bic_csv <- file.path(pilot_design_dir, sprintf("pilot_bic_tables_by_column_seed%04d.csv", pilot_seed))
bic_rds <- file.path(pilot_design_dir, sprintf("pilot_bic_tables_by_column_seed%04d.rds", pilot_seed))
best_csv <- file.path(pilot_design_dir, sprintf("pilot_selected_c_by_column_seed%04d.csv", pilot_seed))
dist_csv <- file.path(pilot_design_dir, sprintf("selected_c_distribution_seed%04d.csv", pilot_seed))

write.csv(pilot_bic_tables_by_column, bic_csv, row.names = FALSE)
saveRDS(pilot_bic_tables_by_column, bic_rds)
write.csv(pilot_selected_c_by_column, best_csv, row.names = FALSE)
write.csv(selected_c_distribution, dist_csv, row.names = FALSE)

cat("Saved column-specific pilot tuning outputs:\n")
cat(bic_csv, "\n")
cat(bic_rds, "\n")
cat(best_csv, "\n")
cat(dist_csv, "\n")
