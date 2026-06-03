#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readxl)
  library(patchwork)
})

# ------------------------------------------------------------
# paths
# ------------------------------------------------------------
setwd('/Users/yuxisong/Library/CloudStorage/Box-Box/Cox-SieveTL/simulation/results/new_icb_response_result')

source("./codes/cv_source_selection.R")

results_dir   <- "./results_allTMB_new_hessian_lasso_new"
grid_file     <- file.path(results_dir, "c_lambda_tuning_grid_results_bic.csv")
data_file     <- "./codes/extracted_cancer_data_by_type.xlsx"
selected_file <- file.path(results_dir, "selected_sources_summary.csv")

stopifnot(file.exists(grid_file))
stopifnot(file.exists(data_file))
stopifnot(file.exists(selected_file))

# ------------------------------------------------------------
# read grid results
# ------------------------------------------------------------
all_sheets <- paste0("cancer", c(1:10, 12:16))
cv_tbl <- read.csv("./results_allTMB_new_hessian_lasso_new/cv_metrics_summary_table_with_c_lambda_bic.csv", stringsAsFactors = FALSE)
c_map <- cv_tbl$c_mult
names(c_map) <- as.character(all_sheets)
grid_df <- read.csv(grid_file, stringsAsFactors = FALSE)

c_df <- tibble(
  target_cancer = names(c_map),
  c_mult_keep = as.numeric(c_map)
)

grid_df <- grid_df %>%
  left_join(c_df, by = "target_cancer") %>%
  mutate(
    c_mult_keep = ifelse(
      is.na(c_mult_keep),
      c_default,
      c_mult_keep
    )
  ) %>%
  filter(c_mult == c_mult_keep) %>%
  dplyr::select(-c_mult_keep)

# optional: remove some lambda_eta values from plotting
grid_df <- grid_df %>%
  filter(!lambda_eta %in% c(0.1, 0.5))

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

plot_df <- grid_df %>%
  mutate(
    cancer_label = ifelse(
      target_cancer %in% names(title_map),
      unname(title_map[target_cancer]),
      target_cancer
    )
    ,
    BIC_plot = ifelse(is.finite(BIC), BIC, NA_real_)
  )

# ------------------------------------------------------------
# helper functions
# ------------------------------------------------------------
parse_sources <- function(s) {
  if (is.null(s) || length(s) == 0 || is.na(s) || trimws(s) == "") return(character(0))
  trimws(strsplit(s, ",")[[1]])
}

get_fix_factors <- function(tname, proc, selected_df, c_map, c_default = 1) {
  if (!tname %in% names(proc)) {
    stop("Target cancer not found in proc: ", tname)
  }
  
  target <- proc[[tname]]
  target_full_df <- target$df
  
  sel_row <- selected_df %>%
    filter(target_cancer == tname)
  
  sel_sources <- character(0)
  if (nrow(sel_row) > 0 &&
      !is.na(sel_row$n_selected[1]) &&
      sel_row$n_selected[1] > 0) {
    sel_sources <- parse_sources(sel_row$selected_sources[1])
    sel_sources <- sel_sources[sel_sources %in% names(proc) & sel_sources != tname]
  }
  
  source_full_df <- NULL
  if (length(sel_sources) > 0) {
    source_full_df <- bind_rows(lapply(sel_sources, function(sn) proc[[sn]]$df))
  }
  
  # same preprocessing logic as tuning step
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
  
  c_mult <- if (tname %in% names(c_map)) c_map[[tname]] else c_default
  
  n_internal <- max(1, floor(c_mult * length(unique(all_time))^(1/3)))
  pn <- n_internal + 4
  L  <- length(target$feature_cols)
  
  fix_lz <- sqrt(log(pn + L) / nrow(target_full_df))
  fix_le <- sqrt(log(pn + L) / nrow(target_full_df))
  
  list(
    fix_lz = fix_lz,
    fix_le = fix_le,
    c_mult = c_mult,
    n_internal = n_internal,
    pn = pn,
    L = L,
    n_target = nrow(target_full_df),
    n_unique_time = length(unique(all_time))
  )
}

# prepare_heatmap_df <- function(df_sub, fix_lz, fix_le) {
#   df_sub <- df_sub %>%
#     mutate(
#       lambda_zeta_raw = lambda_zeta / fix_lz,
#       lambda_eta_raw  = lambda_eta  / fix_le
#     )
#   
#   lambda_zeta_levels <- format(sort(unique(df_sub$lambda_zeta_raw)), scientific = TRUE)
#   lambda_eta_levels  <- format(sort(unique(df_sub$lambda_eta_raw)), scientific = TRUE)
#   
#   df_sub %>%
#     mutate(
#       lambda_zeta_f = factor(
#         format(lambda_zeta_raw, scientific = TRUE),
#         levels = lambda_zeta_levels
#       ),
#       lambda_eta_f = factor(
#         format(lambda_eta_raw, scientific = TRUE),
#         levels = lambda_eta_levels
#       )
#     )
# }
prepare_heatmap_df <- function(df_sub, fix_lz, fix_le) {
  
  df_sub <- df_sub %>%
    mutate(
      lambda_zeta_raw = lambda_zeta / fix_lz,
      lambda_eta_raw  = lambda_eta  / fix_le,
      BIC_plot = ifelse(is.finite(BIC), BIC, Inf)
    )
  
  # --------------------------------------------------
  # ❗ STRICT RULE:
  # remove lambda_zeta if ANY BIC is NA
  # --------------------------------------------------
  valid_zeta <- df_sub %>%
    group_by(lambda_zeta_raw) %>%
    summarise(any_na = any(!is.finite(BIC)), .groups = "drop") %>%
    #filter(!any_na) %>%
    pull(lambda_zeta_raw)
  
  # --------------------------------------------------
  # ❗ remove lambda_eta if ANY BIC is NA
  # --------------------------------------------------
  valid_eta <- df_sub %>%
    group_by(lambda_eta_raw) %>%
    summarise(any_na = any(!is.finite(BIC)), .groups = "drop") %>%
    #filter(!any_na) %>%
    pull(lambda_eta_raw)
  
  df_sub <- df_sub %>%
    filter(
      lambda_zeta_raw %in% valid_zeta,
      lambda_eta_raw  %in% valid_eta
    )
  
  # --------------------------------------------------
  # factor levels AFTER filtering
  # --------------------------------------------------
  lambda_zeta_levels <- format(sort(unique(df_sub$lambda_zeta_raw)), scientific = TRUE)
  lambda_eta_levels  <- format(sort(unique(df_sub$lambda_eta_raw)), scientific = TRUE)
  
  df_sub %>%
    mutate(
      lambda_zeta_f = factor(
        format(lambda_zeta_raw, scientific = TRUE),
        levels = lambda_zeta_levels
      ),
      lambda_eta_f = factor(
        format(lambda_eta_raw, scientific = TRUE),
        levels = lambda_eta_levels
      )
    )
}
make_bic_heatmap <- function(df_sub, title_text) {
  best_row <- df_sub %>%
    filter(is.finite(BIC)) %>%
    arrange(BIC, lambda_zeta_raw, lambda_eta_raw) %>%
    slice(1)
  
  best_key <- NULL
  if (nrow(best_row) == 1) {
    best_key <- paste(best_row$lambda_zeta_f, best_row$lambda_eta_f)
  }
  
  df_sub <- df_sub %>%
    mutate(
      is_best = if (!is.null(best_key)) {
        paste(lambda_zeta_f, lambda_eta_f) == best_key
      } else {
        FALSE
      },
      label_text = case_when(
        is.na(BIC) ~ "NA",
        is.infinite(BIC) ~ "Inf",
        TRUE ~ sprintf("%.3f", BIC)
      )
    )
  
  ggplot(df_sub, aes(x = lambda_eta_f, y = lambda_zeta_f, fill = BIC_plot)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_point(
      data = df_sub %>% filter(is_best),
      aes(x = lambda_eta_f, y = lambda_zeta_f),
      inherit.aes = FALSE,
      shape = 21,
      size = 8,
      stroke = 1.5,
      fill = NA,
      color = "red"
    ) +
    geom_text(
      aes(
        label = label_text,
        fontface = ifelse(is_best, "bold", "plain")
      ),
      size = 2
    ) +
    scale_fill_gradient(
      low = "#deebf7",
      high = "#08519c",
      na.value = "grey92"
    ) +
    labs(
      title = title_text,
      x = expression(c[eta]),
      y = expression(c[zeta]),
      fill = "BIC"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7)
    )
}

# ------------------------------------------------------------
# load cancer data to recover fix_lz / fix_le
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# 1) single cancer: cancer5
# ------------------------------------------------------------
single_cancer <- "cancer5"

single_df <- plot_df %>%
  filter(target_cancer == single_cancer)

if (nrow(single_df) > 0) {
  fix_obj_single <- get_fix_factors(
    tname = single_cancer,
    proc = proc,
    selected_df = selected_df,
    c_map = c_map,
    c_default = c_default
  )
  
  single_df_plot <- prepare_heatmap_df(
    df_sub = single_df,
    fix_lz = fix_obj_single$fix_lz,
    fix_le = fix_obj_single$fix_le
  )
  
  p_single <- make_bic_heatmap(
    single_df_plot,
    title_text = paste0("Cox-SieveTL BIC heatmap: ", unique(single_df_plot$cancer_label))
  )
  
  ggsave(
    filename = file.path(results_dir, paste0("bic_heatmap_", single_cancer, ".pdf")),
    plot = p_single,
    width = 7,
    height = 5.5
  )
}

# ------------------------------------------------------------
# 2) all cancers combined
# ------------------------------------------------------------
all_cancers <- unique(plot_df$target_cancer)
all_cancers <- all_cancers[order(match(all_cancers, names(title_map)))]

plot_list <- vector("list", length(all_cancers))

for (i in seq_along(all_cancers)) {
  cc <- all_cancers[i]
  
  df_sub <- plot_df %>%
    filter(target_cancer == cc)
  
  if (nrow(df_sub) == 0) next
  if (!cc %in% names(proc)) next
  
  fix_obj <- get_fix_factors(
    tname = cc,
    proc = proc,
    selected_df = selected_df,
    c_map = c_map,
    c_default = c_default
  )
  
  df_sub_plot <- prepare_heatmap_df(
    df_sub = df_sub,
    fix_lz = fix_obj$fix_lz,
    fix_le = fix_obj$fix_le
  )
  
  plot_list[[i]] <- make_bic_heatmap(
    df_sub_plot,
    title_text = unique(df_sub_plot$cancer_label)
  )
}

plot_list <- plot_list[!vapply(plot_list, is.null, logical(1))]

combined_plot <- wrap_plots(plotlist = plot_list, ncol = 3) +
  plot_annotation(title = "Cox-SieveTL BIC heatmaps by cancer")

ggsave(
  filename = file.path(results_dir, "bic_heatmap_all_cancers_circle_best.pdf"),
  plot = combined_plot,
  width = 18,
  height = 16
)
