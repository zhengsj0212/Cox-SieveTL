setwd('/Users/yuxisong/Library/CloudStorage/Box-Box/Cox-SieveTL/simulation/results/new_icb_response_result')

source("./results_allTMB_new_hessian_lasso_new/codes/cv_source_selection.R")
source("./results_allTMB_new_hessian_lasso_new/codes/cv_eval_selected_source_bic.R")
# source("./codes/forest_plot_try.R")

data_file   <- "./codes/extracted_cancer_data_by_type.xlsx"
results_dir <- "./results_allTMB_new_hessian_lasso_new"

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

covariate_name <- "LORIS"

x_min <- 0
x_max <- 4

get_sievetl_plugin_transfer_results <- function(covariate_name,
                                                display_name = "Cox-SieveTL") {
  if(display_name == 'Cox-SieveTL') {
    infer_file <- './results_allTMB_new_hessian_lasso_new/realdata_transfer_hessian_inference/all_targets_beta_wald_ci.csv'
  } else {
    infer_file <- file.path(results_dir, "inference_beta_summary_selected_methods.csv")
  }
  if (!file.exists(infer_file)) {
    stop("Cannot find the file in results_dir.")
  }
  
  infer_df <- read.csv(infer_file, stringsAsFactors = FALSE)
  
  # If your file uses target_cancer instead of cancer, rename it
  if (!"cancer" %in% names(infer_df) && "target_cancer" %in% names(infer_df)) {
    infer_df <- infer_df %>% rename(cancer = target_cancer)
  }
  
  if (!"cancer" %in% names(infer_df)) {
    stop("infer_df must contain a column named 'cancer' or 'target_cancer'.")
  }
  
  sheets_available <- excel_sheets(data_file)
  sheets <- intersect(names(title_map), sheets_available)
  
  cancers <- lapply(sheets, function(sh) {
    df <- read_excel(data_file, sheet = sh)
    df <- cap_cancer_df(df)
    df %>% drop_na()
  })
  names(cancers) <- sheets
  proc <- lapply(cancers, prep_df)
  
  n_map <- sapply(proc, function(x) nrow(x$df))
  
  sub <- infer_df %>%
    filter(.data$covariate == covariate_name) %>%
    mutate(
      cancer_label = dplyr::if_else(
        .data$cancer %in% names(title_map),
        unname(title_map[.data$cancer]),
        .data$cancer
      ),
      Samples = as.numeric(n_map[.data$cancer]),
      beta = .data$estimate,
      HR = exp(.data$estimate),
      CI_low = exp(.data$ci_lower),
      CI_high = exp(.data$ci_upper)
    ) %>%
    dplyr::select(cancer, cancer_label, Samples, beta, se, z, p_value, HR, ci_lower, ci_upper)
  
  missing_cancers <- setdiff(names(title_map), sub$cancer)
  
  if (length(missing_cancers) > 0) {
    sub <- bind_rows(
      sub,
      tibble(
        cancer = missing_cancers,
        cancer_label = unname(title_map[missing_cancers]),
        Samples = as.numeric(n_map[missing_cancers]),
        beta = NA_real_,
        se = NA_real_,
        z = NA_real_,
        pvalue = NA_real_,
        HR = NA_real_,
        CI_low = NA_real_,
        CI_high = NA_real_
      )
    )
  }
  
  pooled <- fixed_effect_pool(sub, beta_col = "beta", se_col = "se")
  
  sub <- bind_rows(
    sub,
    tibble(
      cancer = "All",
      cancer_label = "All",
      Samples = sum(sub$Samples, na.rm = TRUE),
      beta = pooled$beta,
      se = pooled$se,
      z = pooled$z,
      pvalue = pooled$pvalue,
      HR = pooled$HR,
      CI_low = pooled$CI_low,
      CI_high = pooled$CI_high
    )
  )
  
  sub$method_display <- display_name
  sub
}
get_size_map <- function() {
  sheets_available <- excel_sheets(data_file)
  sheets <- intersect(names(title_map), sheets_available)
  
  cancers <- lapply(sheets, function(sh) {
    df <- read_excel(data_file, sheet = sh)
    df <- cap_cancer_df(df)
    df %>% drop_na()
  })
  names(cancers) <- sheets
  proc <- lapply(cancers, prep_df)
  
  selected_df <- read.csv(file.path(results_dir, "selected_sources_summary.csv"), stringsAsFactors = FALSE)
  
  size_df <- lapply(names(proc), function(tname) {
    target <- proc[[tname]]
    target_df <- target$df
    
    sel_row <- selected_df[selected_df$target_cancer == tname, , drop = FALSE]
    sel_sources <- character(0)
    if (nrow(sel_row) > 0 &&
        !is.na(sel_row$n_selected[1]) &&
        sel_row$n_selected[1] > 0) {
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
    
    tibble(
      cancer = tname,
      n_target = nrow(pre$train_df),
      n_source = if (is.null(source_df_pool)) 0 else nrow(pre$apply_df_list$src)
    )
  }) %>% bind_rows()
  
  size_df
}

get_coxph_target_results <- function(covariate_name) {
  sheets_available <- excel_sheets(data_file)
  sheets <- intersect(names(title_map), sheets_available)
  
  cancers <- lapply(sheets, function(sh) {
    df <- read_excel(data_file, sheet = sh)
    df <- cap_cancer_df(df)
    df %>% drop_na()
  })
  names(cancers) <- sheets
  
  proc <- lapply(cancers, prep_df)
  
  rows <- list()
  
  for (cc in names(proc)) {
    obj <- proc[[cc]]
    dat <- obj$df
    feature_cols <- obj$feature_cols
    
    n_target <- nrow(dat)
    
    if (!(covariate_name %in% feature_cols)) {
      rows[[length(rows) + 1]] <- tibble(
        cancer = cc,
        cancer_label = unname(title_map[cc]),
        Samples = n_target,
        beta = NA_real_,
        se = NA_real_,
        z = NA_real_,
        pvalue = NA_real_,
        HR = NA_real_,
        CI_low = NA_real_,
        CI_high = NA_real_
      )
      next
    }
    
    fml <- as.formula(
      paste0(
        "Surv(", obj$time_col, ",", obj$event_col, ") ~ ",
        paste(feature_cols, collapse = " + ")
      )
    )
    
    fit <- tryCatch(
      survival::coxph(fml, data = dat),
      error = function(e) NULL
    )
    
    if (is.null(fit)) {
      rows[[length(rows) + 1]] <- tibble(
        cancer = cc,
        cancer_label = unname(title_map[cc]),
        Samples = n_target,
        beta = NA_real_,
        se = NA_real_,
        z = NA_real_,
        pvalue = NA_real_,
        HR = NA_real_,
        CI_low = NA_real_,
        CI_high = NA_real_
      )
      next
    }
    
    sm <- summary(fit)
    coef_tab <- as.data.frame(sm$coefficients)
    coef_tab$term <- rownames(coef_tab)
    
    this_row <- coef_tab %>% filter(term == covariate_name)
    
    if (nrow(this_row) == 0) {
      rows[[length(rows) + 1]] <- tibble(
        cancer = cc,
        cancer_label = unname(title_map[cc]),
        Samples = n_target,
        beta = NA_real_,
        se = NA_real_,
        z = NA_real_,
        pvalue = NA_real_,
        HR = NA_real_,
        CI_low = NA_real_,
        CI_high = NA_real_
      )
      next
    }
    
    beta_hat <- this_row$coef[1]
    se_hat   <- this_row$`se(coef)`[1]
    z_hat    <- this_row$`z`[1]
    p_hat    <- this_row$`Pr(>|z|)`[1]
    
    rows[[length(rows) + 1]] <- tibble(
      cancer = cc,
      cancer_label = unname(title_map[cc]),
      Samples = n_target,
      beta = beta_hat,
      se = se_hat,
      z = z_hat,
      pvalue = p_hat,
      HR = exp(beta_hat),
      CI_low = exp(beta_hat - 1.96 * se_hat),
      CI_high = exp(beta_hat + 1.96 * se_hat)
    )
  }
  
  out <- bind_rows(rows)
  
  pooled <- fixed_effect_pool(out, beta_col = "beta", se_col = "se")
  out <- bind_rows(
    out,
    tibble(
      cancer = "All",
      cancer_label = "All",
      Samples = sum(out$Samples, na.rm = TRUE),
      beta = pooled$beta,
      se = pooled$se,
      z = pooled$z,
      pvalue = pooled$pvalue,
      HR = pooled$HR,
      CI_low = pooled$CI_low,
      CI_high = pooled$CI_high
    )
  )
  
  out
}

fixed_effect_pool <- function(df, beta_col = "beta", se_col = "se") {
  dd <- df %>%
    filter(is.finite(.data[[beta_col]]), is.finite(.data[[se_col]]), .data[[se_col]] > 0)
  
  if (nrow(dd) == 0) {
    return(tibble(
      beta = NA_real_,
      se = NA_real_,
      z = NA_real_,
      pvalue = NA_real_,
      HR = NA_real_,
      CI_low = NA_real_,
      CI_high = NA_real_
    ))
  }
  
  w <- 1 / (dd[[se_col]]^2)
  beta_pool <- sum(w * dd[[beta_col]]) / sum(w)
  se_pool <- sqrt(1 / sum(w))
  z_pool <- beta_pool / se_pool
  p_pool <- 2 * pnorm(-abs(z_pool))
  
  tibble(
    beta = beta_pool,
    se = se_pool,
    z = z_pool,
    pvalue = p_pool,
    HR = exp(beta_pool),
    CI_low = exp(beta_pool - 1.96 * se_pool),
    CI_high = exp(beta_pool + 1.96 * se_pool)
  )
}

get_bootstrap_method_results <- function(covariate_name,
                                         method_name,
                                         display_name) {
  infer_file <- file.path(results_dir, "inference_beta_summary_selected_methods.csv")
  if (!file.exists(infer_file)) {
    stop("Cannot find inference_beta_summary_selected_methods.csv in results_dir.")
  }
  
  infer_df <- read.csv(infer_file, stringsAsFactors = FALSE)
  
  # sample sizes from processed target datasets
  sheets_available <- excel_sheets(data_file)
  sheets <- intersect(names(title_map), sheets_available)
  
  cancers <- lapply(sheets, function(sh) {
    df <- read_excel(data_file, sheet = sh)
    df <- cap_cancer_df(df)
    df %>% drop_na()
  })
  names(cancers) <- sheets
  proc <- lapply(cancers, prep_df)
  
  n_map <- sapply(proc, function(x) nrow(x$df))
  
  sub <- infer_df %>%
    filter(method == method_name, covariate == covariate_name) %>%
    mutate(
      cancer_label = ifelse(target_cancer %in% names(title_map), unname(title_map[target_cancer]), target_cancer),
      Samples = as.numeric(n_map[target_cancer]),
      beta = beta_hat,
      HR = exp(beta_hat),
      CI_low = exp(ci_lo),
      CI_high = exp(ci_hi)
    ) %>%
    dplyr::select(target_cancer, cancer_label, Samples, beta, se, z, pvalue, HR, CI_low, CI_high)
  
  # add missing cancers if any
  missing_cancers <- setdiff(names(title_map), sub$target_cancer)
  if (length(missing_cancers) > 0) {
    sub <- bind_rows(
      sub,
      tibble(
        target_cancer = missing_cancers,
        cancer_label = unname(title_map[missing_cancers]),
        Samples = as.numeric(n_map[missing_cancers]),
        beta = NA_real_,
        se = NA_real_,
        z = NA_real_,
        pvalue = NA_real_,
        HR = NA_real_,
        CI_low = NA_real_,
        CI_high = NA_real_
      )
    )
  }
  
  pooled <- fixed_effect_pool(sub, beta_col = "beta", se_col = "se")
  
  sub <- bind_rows(
    sub,
    tibble(
      target_cancer = "All",
      cancer_label = "All",
      Samples = sum(sub$Samples, na.rm = TRUE),
      beta = pooled$beta,
      se = pooled$se,
      z = pooled$z,
      pvalue = pooled$pvalue,
      HR = pooled$HR,
      CI_low = pooled$CI_low,
      CI_high = pooled$CI_high
    )
  )
  
  sub$method_display <- display_name
  sub
}
get_sievetl_hessian_results <- function(covariate_name,
                                        display_name = "Cox-SieveTL Hessian") {
  
  infer_file <- file.path(
    results_dir,
    "realdata_hessian_beta_cumhaz_inference",
    "all_targets_beta_wald_ci_hessian.csv"
  )
  
  if (!file.exists(infer_file)) {
    stop("Cannot find: ", infer_file)
  }
  
  infer_df <- read.csv(infer_file, stringsAsFactors = FALSE)
  
  if (!"cancer" %in% names(infer_df) && "target_cancer" %in% names(infer_df)) {
    infer_df <- infer_df %>% rename(cancer = target_cancer)
  }
  
  sheets_available <- excel_sheets(data_file)
  sheets <- intersect(names(title_map), sheets_available)
  
  cancers <- lapply(sheets, function(sh) {
    df <- read_excel(data_file, sheet = sh)
    df <- cap_cancer_df(df)
    df %>% drop_na()
  })
  names(cancers) <- sheets
  proc <- lapply(cancers, prep_df)
  
  n_map <- sapply(proc, function(x) nrow(x$df))
  
  sub <- infer_df %>%
    filter(covariate == covariate_name) %>%
    mutate(
      cancer_label = ifelse(
        cancer %in% names(title_map),
        unname(title_map[cancer]),
        cancer
      ),
      Samples = as.numeric(n_map[cancer]),
      beta = estimate,
      HR = exp(estimate),
      CI_low = exp(ci_lower),
      CI_high = exp(ci_upper),
      pvalue = p_value
    ) %>%
    dplyr::select(
      cancer, cancer_label, Samples,
      beta, se, z, pvalue, HR, CI_low, CI_high
    )
  
  missing_cancers <- setdiff(names(title_map), sub$cancer)
  
  if (length(missing_cancers) > 0) {
    sub <- bind_rows(
      sub,
      tibble(
        cancer = missing_cancers,
        cancer_label = unname(title_map[missing_cancers]),
        Samples = as.numeric(n_map[missing_cancers]),
        beta = NA_real_,
        se = NA_real_,
        z = NA_real_,
        pvalue = NA_real_,
        HR = NA_real_,
        CI_low = NA_real_,
        CI_high = NA_real_
      )
    )
  }
  
  pooled <- fixed_effect_pool(sub, beta_col = "beta", se_col = "se")
  
  sub <- bind_rows(
    sub,
    tibble(
      cancer = "All",
      cancer_label = "All",
      Samples = sum(sub$Samples, na.rm = TRUE),
      beta = pooled$beta,
      se = pooled$se,
      z = pooled$z,
      pvalue = pooled$pvalue,
      HR = pooled$HR,
      CI_low = pooled$CI_low,
      CI_high = pooled$CI_high
    )
  )
  
  sub$method_display <- display_name
  sub
}


make_forest_panel <- function(df, panel_title, x_min = 0, x_max = 4) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(cowplot)
  })
  
  # --------------------------------
  # unify cancer column for TransCox / DiscreteKL
  # --------------------------------
  if (!"cancer" %in% names(df) && "target_cancer" %in% names(df)) {
    df <- df %>% rename(cancer = target_cancer)
  }
  
  # if cancer_label missing, create it
  if (!"cancer_label" %in% names(df)) {
    df <- df %>%
      mutate(cancer_label = ifelse(cancer %in% names(title_map), unname(title_map[cancer]), cancer))
  }
  
  # if Samples missing, set NA
  if (!"Samples" %in% names(df)) {
    df <- df %>% mutate(Samples = NA_real_)
  }
  
  df <- df %>%
    filter(!cancer %in% c("cancer11", "cancer18", "All")) %>%
    filter(!is.na(cancer_label))
  
  order_vec <- setdiff(names(title_map), "cancer18")
  order_vec <- setdiff(order_vec, "cancer11")
  label_vec <- unname(title_map[order_vec])
  
  df <- df %>%
    mutate(
      cancer = factor(cancer, levels = order_vec),
      cancer_label = factor(cancer_label, levels = label_vec),
      hr_label = ifelse(is.finite(HR) & is.finite(CI_low) & is.finite(CI_high),
                        format_hr(HR, CI_low, CI_high), ""),
      p_label = ifelse(is.finite(pvalue), format_pval(pvalue), "")
    ) %>%
    arrange(cancer)
  
  n <- nrow(df)
  df$y <- rev(seq_len(n))
  
  # ---------- Left table ----------
  p_left <- ggplot(df, aes(y = y)) +
    geom_text(aes(x = 0.02, label = cancer_label), hjust = 0, size = 6) +
    geom_text(aes(x = 0.98, label = Samples), hjust = 1, size = 6) +
    annotate("text", x = 0.02, y = max(df$y) + 1.1, label = "Cancer type",
             hjust = 0, fontface = "bold", size = 6) +
    annotate("text", x = 0.98, y = max(df$y) + 1.1, label = "Samples",
             hjust = 1, fontface = "bold", size = 6) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, max(df$y) + 1.5), expand = c(0, 0)) +
    theme_void()
  
  # ---------- Forest Plot ----------
  df_mid <- df %>%
    mutate(
      valid_estimate = is.finite(HR),
      valid_ci = valid_estimate & is.finite(CI_low) & is.finite(CI_high),
      lower_plot = ifelse(valid_ci, pmax(CI_low, x_min), NA_real_),
      upper_plot = ifelse(valid_ci, pmin(CI_high, x_max), NA_real_),
      trunc_right = valid_ci & CI_high > x_max
    )
  
  p_mid <- ggplot(df_mid, aes(y = y)) +
    geom_vline(xintercept = 1, color = "grey70") +
    geom_segment(
      data = df_mid %>% filter(valid_ci, !trunc_right),
      aes(x = lower_plot, xend = upper_plot, yend = y),
      linewidth = 0.6
    ) +
    geom_segment(
      data = df_mid %>% filter(valid_ci, trunc_right),
      aes(x = lower_plot, xend = x_max, yend = y),
      linewidth = 0.6,
      arrow = arrow(length = unit(0.08, "inches"), type = "closed")
    ) +
    geom_point(
      data = df_mid %>% filter(valid_estimate),
      aes(x = HR),
      shape = 15,
      size = 3.8
    ) +
    scale_x_continuous(
      limits = c(x_min, x_max),
      breaks = seq(x_min, x_max, by = 1),
      name = "Hazard ratio"
    ) +
    scale_y_continuous(limits = c(0.5, max(df$y) + 1.5), expand = c(0, 0)) +
    labs(title = panel_title) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 22),
      axis.title.x = element_text(size = 18),
      axis.text.x = element_text(size = 15)
    )
  
  # ---------- Right table ----------
  p_right <- ggplot(df, aes(y = y)) +
    geom_text(aes(x = 0.02, label = hr_label), hjust = 0, size = 6) +
    geom_text(aes(x = 0.98, label = p_label), hjust = 1, size = 6) +
    annotate("text", x = 0.02, y = max(df$y) + 1.1, label = "HR (95% CI)",
             hjust = 0, fontface = "bold", size = 6) +
    annotate("text", x = 0.98, y = max(df$y) + 1.1, label = "P value",
             hjust = 1, fontface = "bold", size = 6) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, max(df$y) + 1.5), expand = c(0, 0)) +
    theme_void()
  
  # ---------- Combine plot ----------
  combined <- cowplot::plot_grid(
    p_left, p_mid, p_right,
    nrow = 1,
    rel_widths = c(1.8, 4.8, 2.6),
    align = "h",
    axis = "tb"
  )
  
  combined
}

cox_df <- get_coxph_target_results(covariate_name) %>%
  mutate(method_display = "CoxPH Target")

sieve_df <- get_sievetl_plugin_transfer_results(
  covariate_name = covariate_name,
  display_name = "Cox-SieveTL"
)

sieve_hessian_df <- get_sievetl_hessian_results(
  covariate_name = covariate_name,
  display_name = "Cox-SieveTL Hessian"
)

transcox_df <- get_bootstrap_method_results(
  covariate_name = covariate_name,
  method_name = "transcox",
  display_name = "TransCox"
) %>%
  rename(cancer = target_cancer)

discretekl_df <- get_bootstrap_method_results(
  covariate_name = covariate_name,
  method_name = "discretekl",
  display_name = "DiscreteKL"
) %>%
  rename(cancer = target_cancer)

sieve_df$pvalue <- sieve_df$p_value
sieve_df$CI_low <- exp(sieve_df$ci_lower)
sieve_df$CI_high <- exp(sieve_df$ci_upper)

all_methods_df <- bind_rows(
  cox_df %>% dplyr::select(cancer, cancer_label, Samples, beta, se, z, pvalue, HR, CI_low, CI_high, method_display),
  sieve_df %>% dplyr::select(cancer, cancer_label, Samples, beta, se, z, pvalue, HR, CI_low, CI_high, method_display),
  sieve_hessian_df %>% dplyr::select(cancer, cancer_label, Samples, beta, se, z, pvalue, HR, CI_low, CI_high, method_display),
  transcox_df %>% dplyr::select(cancer, cancer_label, Samples, beta, se, z, pvalue, HR, CI_low, CI_high, method_display),
  discretekl_df %>% dplyr::select(cancer, cancer_label, Samples, beta, se, z, pvalue, HR, CI_low, CI_high, method_display)
)

size_df <- get_size_map()

all_methods_df <- all_methods_df %>%
  left_join(size_df, by = "cancer")
make_forest_panel_all_methods <- function(
    df,
    panel_title,
    x_min = 0,
    x_max = 4,
    row_spacing = 1.8,
    offset_scale = 0.16
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(cowplot)
  })
  
  df <- df %>%
    filter(!cancer %in% c("cancer11", "cancer18", "All")) %>%
    filter(!is.na(cancer_label))
  
  order_vec <- setdiff(names(title_map), c("cancer11", "cancer18"))
  label_vec <- unname(title_map[order_vec])
  
  method_levels <- c(
    "CoxPH Target",
    "Cox-SieveTL",
    "Cox-SieveTL Hessian",
    "TransCox",
    "DiscreteKL"
  )
  
  df <- df %>%
    mutate(
      cancer = factor(cancer, levels = order_vec),
      cancer_label = factor(cancer_label, levels = label_vec),
      method_display = factor(method_display, levels = method_levels)
    ) %>%
    arrange(cancer, method_display)
  
  cancer_rows <- tibble(
    cancer = factor(order_vec, levels = order_vec),
    cancer_label = factor(unname(title_map[order_vec]), levels = label_vec),
    y = rev(seq_along(order_vec)) * row_spacing
  )
  
  df <- df %>%
    left_join(cancer_rows %>% dplyr::select(cancer, y), by = "cancer") %>%
    group_by(cancer) %>%
    mutate(
      method_index = match(method_display, method_levels),
      y_plot = y + c(-2, -1, 0, 1, 2)[method_index] * offset_scale * row_spacing
    ) %>%
    ungroup() %>%
    mutate(
      valid_estimate = is.finite(HR),
      valid_ci = valid_estimate & is.finite(CI_low) & is.finite(CI_high),
      
      lower_plot = ifelse(valid_ci, pmax(CI_low, x_min), NA_real_),
      upper_plot = ifelse(valid_ci, pmin(CI_high, x_max), NA_real_),
      
      trunc_right = valid_ci & CI_high > x_max,
      trunc_left  = valid_ci & CI_low < x_min,
      
      # at least some part of CI overlaps plotting window
      ci_in_range = valid_ci & !(CI_high < x_min | CI_low > x_max),
      
      # final flag used for geom_segment
      draw_segment = ci_in_range & is.finite(lower_plot) & is.finite(upper_plot) &
        (upper_plot >= lower_plot)
    )
  
  left_df <- df %>%
    distinct(cancer, cancer_label, y, n_target, n_source) %>%
    arrange(cancer)
  
  bg_df <- left_df %>%
    mutate(
      ymin = y - 0.48 * row_spacing,
      ymax = y + 0.48 * row_spacing,
      fill_id = rep(c("a", "b"), length.out = n())
    )
  
  color_map <- c(
    "CoxPH Target" = "#000000",
    "Cox-SieveTL" = "#D55E00",
    "Cox-SieveTL Hessian" = "#CC79A7",
    "TransCox" = "#0072B2",
    "DiscreteKL" = "#009E73"
  )
  
  shape_map <- c(
    "CoxPH Target" = 15,
    "Cox-SieveTL" = 16,
    "Cox-SieveTL Hessian" = 8,
    "TransCox" = 17,
    "DiscreteKL" = 18
  )
  
  # left table
  p_left <- ggplot() +
    geom_rect(
      data = bg_df,
      aes(xmin = 0, xmax = 1, ymin = ymin, ymax = ymax, fill = fill_id),
      inherit.aes = FALSE,
      alpha = 0.6
    ) +
    scale_fill_manual(values = c(a = "#f5f5f5", b = "white"), guide = "none") +
    geom_text(
      data = left_df,
      aes(x = 0.02, y = y, label = cancer_label),
      hjust = 0, size = 6
    ) +
    geom_text(
      data = left_df,
      aes(x = 0.72, y = y, label = n_target),
      hjust = 1, size = 6
    ) +
    geom_text(
      data = left_df,
      aes(x = 0.98, y = y, label = n_source),
      hjust = 1, size = 6
    ) +
    annotate(
      "text", x = 0.02, y = max(left_df$y) + 1.1,
      label = "Cancer type", hjust = 0, fontface = "bold", size = 6
    ) +
    annotate(
      "text", x = 0.72, y = max(left_df$y) + 1.1,
      label = "n_target", hjust = 1, fontface = "bold", size = 6
    ) +
    annotate(
      "text", x = 0.98, y = max(left_df$y) + 1.1,
      label = "n_source", hjust = 1, fontface = "bold", size = 6
    ) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(
      limits = c(min(bg_df$ymin), max(left_df$y) + 0.9 * row_spacing),
      expand = c(0, 0)
    ) +
    theme_void()
  
  # forest panel
  p_mid <- ggplot() +
    geom_rect(
      data = bg_df,
      aes(xmin = x_min, xmax = x_max, ymin = ymin, ymax = ymax, fill = fill_id),
      inherit.aes = FALSE,
      alpha = 0.6
    ) +
    scale_fill_manual(values = c(a = "#f5f5f5", b = "white"), guide = "none") +
    geom_vline(xintercept = 1, color = "grey60", linewidth = 0.6) +
    geom_segment(
      data = df %>% filter(draw_segment, !trunc_right),
      aes(
        x = lower_plot, xend = upper_plot,
        y = y_plot, yend = y_plot,
        color = method_display
      ),
      linewidth = 0.7,
      alpha = 0.95
    ) +
    geom_segment(
      data = df %>% filter(draw_segment, trunc_right),
      aes(
        x = lower_plot, xend = x_max,
        y = y_plot, yend = y_plot,
        color = method_display
      ),
      linewidth = 0.7,
      alpha = 0.95,
      arrow = arrow(length = unit(0.08, "inches"), type = "closed")
    ) +
    geom_point(
      data = df %>% filter(valid_estimate),
      aes(
        x = HR, y = y_plot,
        color = method_display,
        shape = method_display
      ),
      size = 2.8,
      alpha = 0.95
    ) +
    scale_color_manual(values = color_map, breaks = method_levels, drop = FALSE) +
    scale_shape_manual(values = shape_map, breaks = method_levels, drop = FALSE) +
    scale_x_continuous(
      limits = c(x_min, x_max),
      breaks = seq(x_min, x_max, by = 1),
      expand = c(0, 0),
      name = "Hazard ratio"
    ) +
    scale_y_continuous(
      limits = c(min(bg_df$ymin), max(left_df$y) + 0.9 * row_spacing),
      breaks = left_df$y,
      labels = rep("", nrow(left_df)),
      expand = c(0, 0)
    ) +
    labs(title = panel_title, color = NULL, shape = NULL) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "grey60", linewidth = 0.35),
      panel.border = element_rect(linewidth = 0.6, color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
      axis.title.x = element_text(size = 16),
      axis.text.x = element_text(size = 13),
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      legend.key.width = unit(1.2, "cm"),
      legend.box = "horizontal"
    ) +
    guides(
      color = guide_legend(nrow = 1, byrow = TRUE),
      shape = guide_legend(nrow = 1, byrow = TRUE)
    )
  
  cowplot::plot_grid(
    p_left, p_mid,
    nrow = 1,
    rel_widths = c(2.6, 5.8),
    align = "h",
    axis = "tb"
  )
}

forest_all_methods_hessian_plot <- make_forest_panel_all_methods(
  df = all_methods_df,
  panel_title = paste0("Forest Plot of Hazard Ratio for All Methods: ", covariate_name),
  x_min = x_min,
  x_max = x_max,
  row_spacing = 1.8,
  offset_scale = 0.14
)

forest_all_methods_hessian_plot

ggsave(
  file.path(results_dir, paste0("forest_all_methods_with_hessian_", covariate_name, ".pdf")),
  forest_all_methods_hessian_plot,
  width = 15,
  height = 15
)

# ggsave(
#   file.path(results_dir, paste0("forest_all_methods_", covariate_name, "_new.pdf")),
#   forest_all_methods_plot,
#   width = 15,
#   height = 15
# )
