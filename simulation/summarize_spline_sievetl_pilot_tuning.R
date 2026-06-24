#!/usr/bin/env Rscript

pilot_base_out <- path.expand(Sys.getenv("PILOT_OUT_BASE", unset = getwd()))
pilot_search_root <- file.path(pilot_base_out, "pilot_parameter_by_column")
if (!dir.exists(pilot_search_root)) {
  pilot_search_root <- pilot_base_out
}
candidate_files <- list.files(
  pilot_search_root,
  pattern = "^pilot_selected_sievetl_tuning_seed[0-9]+\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
candidate_dirs <- sort(unique(dirname(candidate_files)))

if (length(candidate_dirs) == 0L) {
  stop("No SieveTL pilot tuning files found under ", pilot_search_root)
}

for (pilot_design_dir in candidate_dirs) {
  selected_files <- sort(Sys.glob(file.path(pilot_design_dir, "pilot_selected_sievetl_tuning_seed*.csv")))
  bic_files <- sort(Sys.glob(file.path(pilot_design_dir, "pilot_bic_tables_sievetl_seed*.csv")))

  selected_rows <- lapply(selected_files, read.csv, stringsAsFactors = FALSE)
  selected_df <- do.call(rbind, selected_rows)
  bic_rows <- lapply(bic_files, read.csv, stringsAsFactors = FALSE)
  bic_df <- if (length(bic_rows) > 0L) do.call(rbind, bic_rows) else data.frame()

  final_tuning <- data.frame(
    c_hat = stats::median(selected_df$c_hat),
    c_zeta_hat = stats::median(selected_df$c_zeta_hat),
    c_eta_hat = stats::median(selected_df$c_eta_hat),
    lambda_zeta_hat = stats::median(selected_df$lambda_zeta_hat),
    lambda_eta_hat = stats::median(selected_df$lambda_eta_hat),
    p_hat = stats::median(selected_df$p_hat),
    d_hat = stats::median(selected_df$d_hat),
    n_pilot = nrow(selected_df)
  )

  final_file <- file.path(pilot_design_dir, "final_selected_tuning.csv")
  selected_all_file <- file.path(pilot_design_dir, "pilot_selected_sievetl_tuning_all.csv")
  bic_all_file <- file.path(pilot_design_dir, "pilot_bic_tables_sievetl_all.csv")

  write.csv(final_tuning, final_file, row.names = FALSE)
  write.csv(selected_df, selected_all_file, row.names = FALSE)
  if (nrow(bic_df) > 0L) {
    write.csv(bic_df, bic_all_file, row.names = FALSE)
  }

  cat("Wrote:\n")
  cat(final_file, "\n")
  cat(selected_all_file, "\n")
  if (nrow(bic_df) > 0L) cat(bic_all_file, "\n")
}
