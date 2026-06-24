#!/usr/bin/env Rscript

# Hoffman2 array runner for plugin-Wald coverage of the SieveTL transfer estimator.
#
# Each array task handles one design row and one chunk of Monte Carlo replications.
# The heavy lifting lives in coverage_plugin_wald_sievetl.R; this file only maps
# SGE_TASK_ID to a design/chunk, runs those seeds, and saves one .rds file per seed.
#
# Beta coverage uses the one-step center with influence-function standard error
# from the joint theta influence construction. The cumulative-hazard center is
# Lambda0_os based on gamma_os, and conditional survival is evaluated at
# x0 = rep(0.5, L).

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  flag <- "--file="
  match_idx <- grep(flag, args)
  if (length(match_idx) > 0) {
    return(dirname(normalizePath(sub(flag, "", args[match_idx[1]]))))
  }
  normalizePath(getwd())
}

script_dir <- get_script_dir()
source(file.path(script_dir, "coverage_plugin_wald_sievetl.R"))

format_kappa_tag <- function(kappa) {
  if (isTRUE(all.equal(kappa, round(kappa), tolerance = 1e-8))) {
    return(sprintf("%d", as.integer(round(kappa))))
  }
  formatC(kappa, format = "f", digits = 1)
}

task_id <- as.integer(Sys.getenv("SGE_TASK_ID", "1"))
nsim <- as.integer(Sys.getenv("NSIM", "100"))
reps_per_task <- as.integer(Sys.getenv("REPS_PER_TASK", "5"))
seed_base <- as.integer(Sys.getenv("SEED_BASE", "100000"))
conf_level <- as.numeric(Sys.getenv("CONF_LEVEL", "0.95"))
multiplier_boot <- as.integer(Sys.getenv("MULTIPLIER_BOOT", "1000"))
L <- as.integer(Sys.getenv("LVAL", "5"))
x0 <- rep(0.5, L)
debug_q10 <- Sys.getenv("DEBUG_Q10", "0") %in% c("1", "TRUE", "true", "True")
beta2 <- as.numeric(Sys.getenv("BETA2", "0.5"))
kappa <- as.numeric(Sys.getenv("KAPPA", "2"))
shift_env <- Sys.getenv("SHIFT", "0")
c1 <- as.numeric(Sys.getenv("C1", "1"))
n_target <- as.integer(Sys.getenv("N_TARGET", "100"))
n_source <- as.integer(Sys.getenv("N_SOURCE", "1000"))

if (!is.finite(task_id) || task_id < 1L) stop("Invalid SGE_TASK_ID.")
if (!is.finite(nsim) || nsim < 1L) stop("NSIM must be positive.")
if (!is.finite(reps_per_task) || reps_per_task < 1L) stop("REPS_PER_TASK must be positive.")
if (!is.finite(multiplier_boot) || multiplier_boot < 50L) stop("MULTIPLIER_BOOT must be at least 50.")
if (!is.finite(beta2)) stop("BETA2 must be numeric.")
if (!is.finite(kappa)) stop("KAPPA must be numeric.")
if (!is.finite(c1)) stop("C1 must be numeric.")
if (!is.finite(n_target) || n_target < 1L) stop("N_TARGET must be positive.")
if (!is.finite(n_source) || n_source < 1L) stop("N_SOURCE must be positive.")

shift_flag <- shift_env %in% c("1", "TRUE", "true", "True")
design <- default_sievetl_coverage_design(
  beta2 = beta2,
  kappa = kappa,
  c = c1,
  shift = shift_flag,
  n0 = n_target,
  nA = n_source
)
n_design <- nrow(design)
n_chunks <- ceiling(nsim / reps_per_task)
total_jobs <- n_design * n_chunks

if (task_id > total_jobs) {
  message(sprintf("SGE_TASK_ID=%d exceeds total_jobs=%d. Nothing to run.", task_id, total_jobs))
  quit(save = "no", status = 0)
}

design_idx <- ((task_id - 1L) %% n_design) + 1L
chunk_idx <- ((task_id - 1L) %/% n_design) + 1L
rep_start <- (chunk_idx - 1L) * reps_per_task + 1L
rep_end <- min(chunk_idx * reps_per_task, nsim)
chunk_nsim <- rep_end - rep_start + 1L

row <- design[design_idx, , drop = FALSE]
row_seed_base <- seed_base + 1000000L * (design_idx - 1L) + 10000L * (chunk_idx - 1L)
tuning_row <- resolve_sievetl_tuning_row(
  beta2 = row$beta2[1], kappa = row$kappa[1], c1 = row$c[1],
  n0 = row$n0[1], nA = row$nA[1], L = L,
  shift = isTRUE(row$shift[1]), base_dir = script_dir
)
c_sieve_use <- as.numeric(tuning_row$c_hat[1])

base_root <- path.expand(Sys.getenv("OUT_BASE", unset = getwd()))
if (!dir.exists(base_root)) dir.create(base_root, recursive = TRUE, showWarnings = FALSE)

base_out <- file.path(base_root, sprintf("coverage_n%dN%d", row$n0[1], row$nA[1]))
if (!dir.exists(base_out)) dir.create(base_out, recursive = TRUE, showWarnings = FALSE)

out_dir <- file.path(
  base_out,
  sprintf(
    "beta%.1f_kappa%s_L%d_c%.2f_C%.2f_shift%s",
    row$beta2[1], format_kappa_tag(row$kappa[1]), L, row$c[1], c_sieve_use, ifelse(isTRUE(row$shift[1]), "TRUE", "FALSE")
  )
)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message(sprintf(
  "Task %d/%d: %s, chunk %d/%d, reps %d-%d",
  task_id, total_jobs, basename(out_dir), chunk_idx, n_chunks, rep_start, rep_end
))

for (global_rep_id in rep_start:rep_end) {
  out_file <- file.path(out_dir, sprintf("seed_%04d.rds", global_rep_id))
  rep_seed_base <- seed_base + 1000000L * (design_idx - 1L) + global_rep_id

  ans <- run_sievetl_plugin_coverage_study(
    nsim = 1,
    design = row,
  L = L,
  c_sieve = c_sieve_use,
  lambda_base_dir = script_dir,
  x0 = x0,
  conf_level = conf_level,
    multiplier_boot = multiplier_boot,
    seed_base = rep_seed_base,
    verbose = TRUE,
    debug_q10 = debug_q10
  )

  saveRDS(
    list(
      seed_id = global_rep_id,
      design_row = row,
      summary = ans$summary,
      records = ans$records,
      replications = ans$replications
    ),
    out_file
  )
  message("Saved: ", out_file)
}
