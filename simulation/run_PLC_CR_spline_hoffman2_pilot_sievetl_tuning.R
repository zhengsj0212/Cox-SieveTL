#!/usr/bin/env Rscript

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

c_sieve_grid <- c(0.5, 1)
c_zeta_grid <- c(0, 0.05, 0.1, 0.2)
c_eta_grid <- c(0, 0.05, 0.1, 0.2)
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

D <- simulate_once_all(
  n_target = n_target, n_source = n_source,
  beta2_source = beta2, kappa_source = kappa,
  c1 = c1, covariate_shift = covariate_shift, L = L
)

pilot_tuning <- sievetl_ic_approx(
  target1 = D$target, source1 = D$source,
  L = L,
  c_grid = c_sieve_grid,
  c_zeta_grid = c_zeta_grid,
  c_eta_grid = c_eta_grid,
  lambda_fixed_stage1 = 0
)

selected_row <- data.frame(
  seed = pilot_seed,
  c_hat = pilot_tuning$best_c,
  c_zeta_hat = pilot_tuning$best_c_zeta,
  c_eta_hat = pilot_tuning$best_c_eta,
  lambda_zeta_hat = pilot_tuning$best_lambda_zeta,
  lambda_eta_hat = pilot_tuning$best_lambda_eta,
  p_hat = pilot_tuning$best_p,
  d_hat = pilot_tuning$best_d,
  stringsAsFactors = FALSE
)

bic_csv <- file.path(pilot_design_dir, sprintf("pilot_bic_tables_sievetl_seed%04d.csv", pilot_seed))
sel_csv <- file.path(pilot_design_dir, sprintf("pilot_selected_sievetl_tuning_seed%04d.csv", pilot_seed))

write.csv(pilot_tuning$bic_table, bic_csv, row.names = FALSE)
write.csv(selected_row, sel_csv, row.names = FALSE)

cat("Saved SieveTL pilot tuning outputs:\n")
cat(bic_csv, "\n")
cat(sel_csv, "\n")
