#!/usr/bin/env Rscript

pilot_base_out <- path.expand(Sys.getenv("PILOT_OUT_BASE", unset = getwd()))
pilot_search_root <- file.path(pilot_base_out, "pilot_parameter_by_column")
if (!dir.exists(pilot_search_root)) {
  pilot_search_root <- pilot_base_out
}
candidate_files <- list.files(
  pilot_search_root,
  pattern = "^pilot_selected_c_by_column_seed[0-9]+\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
candidate_dirs <- sort(unique(dirname(candidate_files)))

if (length(candidate_dirs) == 0L) {
  stop("No by-column pilot files found under ", pilot_base_out)
}

for (pilot_design_dir in candidate_dirs) {
  selected_files <- sort(Sys.glob(file.path(pilot_design_dir, "pilot_selected_c_by_column_seed*.csv")))
  bic_files <- sort(Sys.glob(file.path(pilot_design_dir, "pilot_bic_tables_by_column_seed*.csv")))

  tuning_file <- file.path(pilot_design_dir, "final_selected_tuning.csv")
  if (!file.exists(tuning_file)) {
    stop("Stage 1 final_selected_tuning.csv not found in ", pilot_design_dir)
  }
  final_tuning <- read.csv(tuning_file, stringsAsFactors = FALSE)

  selected_rows <- lapply(selected_files, read.csv, stringsAsFactors = FALSE)
  selected_df <- do.call(rbind, selected_rows)

  full_rows <- lapply(bic_files, read.csv, stringsAsFactors = FALSE)
  full_df <- do.call(rbind, full_rows)
  n_match_pilot <- nrow(selected_df) / length(unique(selected_df$j))

  final_by_column <- selected_df |>
    dplyr::group_by(j, best_c) |>
    dplyr::summarise(
      freq_selected = dplyr::n(),
      mean_best_criterion = mean(best_criterion),
      mean_best_residual = mean(best_residual),
      mean_best_l0_delta = mean(best_l0_delta),
      best_lambda_value = mean(best_lambda_value),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      j,
      dplyr::desc(freq_selected),
      best_c
    ) |>
    dplyr::group_by(j) |>
    dplyr::slice(1L) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      prop_selected = freq_selected / n_match_pilot,
      best_lambda_formula = "c * sqrt(log(d_n) / n_target)",
      c_hat = final_tuning$c_hat[1],
      c_zeta_hat = final_tuning$c_zeta_hat[1],
      c_eta_hat = final_tuning$c_eta_hat[1],
      lambda_zeta_hat = final_tuning$lambda_zeta_hat[1],
      lambda_eta_hat = final_tuning$lambda_eta_hat[1],
      d_hat = final_tuning$d_hat[1],
      tie_break_rule_used = paste(
        "highest frequency across pilot seeds;",
        "ties broken by smaller c"
      )
    ) |>
    dplyr::select(
      c_hat, c_zeta_hat, c_eta_hat, lambda_zeta_hat, lambda_eta_hat, d_hat,
      j, best_c, best_lambda_formula, best_lambda_value,
      freq_selected, prop_selected, mean_best_criterion, mean_best_residual,
      mean_best_l0_delta, tie_break_rule_used
    )

  selection_distribution <- final_by_column |>
    dplyr::group_by(best_c) |>
    dplyr::summarise(
      freq_selected = dplyr::n(),
      prop_selected = dplyr::n() / nrow(final_by_column),
      mean_selected_criterion = mean(mean_best_criterion),
      mean_selected_residual = mean(mean_best_residual),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(freq_selected), best_c)

  final_file <- file.path(pilot_design_dir, "final_selected_c_by_column.csv")
  selected_all_file <- file.path(pilot_design_dir, "pilot_selected_c_by_column_all.csv")
  full_all_file <- file.path(pilot_design_dir, "pilot_bic_tables_by_column_all.csv")
  dist_file <- file.path(pilot_design_dir, "selected_c_distribution_across_columns.csv")

  write.csv(final_by_column, final_file, row.names = FALSE)
  write.csv(selected_df, selected_all_file, row.names = FALSE)
  write.csv(full_df, full_all_file, row.names = FALSE)
  write.csv(selection_distribution, dist_file, row.names = FALSE)

  cat("Wrote:\n")
  cat(final_file, "\n")
  cat(selected_all_file, "\n")
  cat(full_all_file, "\n")
  cat(dist_file, "\n")
  cat(tuning_file, "\n")
}
