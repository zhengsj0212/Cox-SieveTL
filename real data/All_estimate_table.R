#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(knitr)
  library(kableExtra)
  library(stringr)
})

results_dir <- "./results_allTMB_new_hessian_lasso"

# =========================================================
# file paths
# =========================================================
infer_file      <- file.path(results_dir, "inference_beta_summary_selected_methods.csv")
all_beta_file   <- file.path(results_dir, "all_methods_beta_estimates.csv")
sievetl_file    <- '/Users/yuxisong/Library/CloudStorage/Box-Box/Cox-SieveTL/simulation/results/new_icb_response_result/results_allTMB_new_hessian_lasso/realdata_transfer_hessian_inference/all_targets_beta_wald_ci.csv'

boot_file       <- file.path(results_dir, "bootstrap_beta_all_transcox_discretekl.csv")

stopifnot(file.exists(infer_file))
stopifnot(file.exists(all_beta_file))
stopifnot(file.exists(boot_file))

infer_df <- read.csv(infer_file, stringsAsFactors = FALSE)
all_beta_df <- read.csv(all_beta_file, stringsAsFactors = FALSE)
boot_df <- read.csv(boot_file, stringsAsFactors = FALSE)

sievetl_df <- if (file.exists(sievetl_file)) {
  read.csv(sievetl_file, stringsAsFactors = FALSE)
} else {
  NULL
}

# =========================================================
# cancer display names
# =========================================================
cancer_map <- c(
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

cancer_order <- c(
  "cancer1","cancer2","cancer3","cancer4","cancer5","cancer6","cancer7",
  "cancer8","cancer9","cancer10","cancer11","cancer12","cancer13",
  "cancer14","cancer16","cancer15"
)

method_order <- c(
  "Std CoxPH (Target)",
  "Std CoxPH (Combined)",
  "Cox-SieveTL",
  "TransCox",
  "DiscreteKL"
)

# =========================================================
# helper formatting
# =========================================================
fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x) | !is.finite(x), NA_character_, sprintf(paste0("%.", digits, "f"), x))
}

latex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([%&_#${}])", "\\\\\\1", x)
  x <- gsub("~", "\\\\textasciitilde{}", x)
  x <- gsub("\\^", "\\\\textasciicircum{}", x)
  x
}

fmt_p <- function(x) {
  ifelse(
    is.na(x) | !is.finite(x),
    "--",
    ifelse(x < 0.001, "<0.001", formatC(x, digits = 3, format = "f"))
  )
}

p_star <- function(p) {
  case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "^{***}",
    p < 0.01  ~ "^{**}",
    p < 0.05  ~ "^{*}",
    p < 0.1   ~ "^{.}",
    TRUE      ~ ""
  )
}

# ---------------------------------------------------------
# p-value stars for CoxPH / SieveTL
# ---------------------------------------------------------
sig_star_from_p <- function(p) {
  out <- rep("", length(p))
  out[!is.na(p) & p < 0.1]   <- "."
  out[!is.na(p) & p < 0.05]  <- "*"
  out[!is.na(p) & p < 0.01]  <- "**"
  out[!is.na(p) & p < 0.001] <- "***"
  out
}

# ---------------------------------------------------------
# bootstrap CI stars for TransCox / DiscreteKL
# 95% CI excludes 0   -> *
# 99% CI excludes 0   -> **
# 99.9% CI excludes 0 -> ***
# ---------------------------------------------------------
sig_star_from_boot_ci <- function(ci95_lo, ci95_hi, ci99_lo, ci99_hi, ci999_lo, ci999_hi) {
  out <- rep("", length(ci95_lo))
  
  idx95 <- is.finite(ci95_lo) & is.finite(ci95_hi) & (ci95_lo > 0 | ci95_hi < 0)
  idx99 <- is.finite(ci99_lo) & is.finite(ci99_hi) & (ci99_lo > 0 | ci99_hi < 0)
  idx999 <- is.finite(ci999_lo) & is.finite(ci999_hi) & (ci999_lo > 0 | ci999_hi < 0)
  
  out[idx95] <- "*"
  out[idx99] <- "**"
  out[idx999] <- "***"
  out
}

# ---------------------------------------------------------
# method-specific formatter
# ---------------------------------------------------------
fmt_beta_ci_with_star <- function(beta, ci_lo, ci_hi, star, digits = 3) {
  beta_chr <- fmt_num(beta, digits)
  
  ifelse(
    is.na(beta) | !is.finite(beta),
    NA_character_,
    ifelse(
      is.na(ci_lo) | is.na(ci_hi) | !is.finite(ci_lo) | !is.finite(ci_hi),
      paste0(beta_chr, star),
      paste0(beta_chr, " [", fmt_num(ci_lo, digits), ", ", fmt_num(ci_hi, digits), "]", star)
    )
  )
}

# =========================================================
# bootstrap summaries for TransCox / DiscreteKL
# =========================================================
covariate_name <- "LORIS"

boot_sig_df <- boot_df %>%
  filter(
    covariate == covariate_name,
    method %in% c("transcox", "discretekl")
  ) %>%
  group_by(target_cancer, method) %>%
  summarise(
    beta_boot_mean = mean(beta, na.rm = TRUE),
    ci95_lo  = quantile(beta, 0.025,  na.rm = TRUE, names = FALSE, type = 7),
    ci95_hi  = quantile(beta, 0.975,  na.rm = TRUE, names = FALSE, type = 7),
    ci99_lo  = quantile(beta, 0.005,  na.rm = TRUE, names = FALSE, type = 7),
    ci99_hi  = quantile(beta, 0.995,  na.rm = TRUE, names = FALSE, type = 7),
    ci999_lo = quantile(beta, 0.0005, na.rm = TRUE, names = FALSE, type = 7),
    ci999_hi = quantile(beta, 0.9995, na.rm = TRUE, names = FALSE, type = 7),
    .groups = "drop"
  ) %>%
  mutate(
    star_boot = sig_star_from_boot_ci(ci95_lo, ci95_hi, ci99_lo, ci99_hi, ci999_lo, ci999_hi)
  )

# =========================================================
# 1) methods with inference table
# =========================================================
loris_infer <- infer_df %>%
  filter(covariate == "LORIS") %>%
  transmute(
    target_cancer,
    method,
    beta = beta_hat,
    ci_lo = ci_lo,
    ci_hi = ci_hi,
    pvalue = pvalue
  )

loris_infer$method <- recode(
  loris_infer$method,
  coxph_target   = "Std CoxPH (Target)",
  coxph_combined = "Std CoxPH (Combined)",
  transcox       = "TransCox",
  discretekl     = "DiscreteKL"
)

# add stars:
# - CoxPH methods from p-values
# - TransCox / DiscreteKL from bootstrap CIs
loris_infer <- loris_infer %>%
  left_join(
    boot_sig_df %>%
      transmute(
        target_cancer,
        method = recode(method, transcox = "TransCox", discretekl = "DiscreteKL"),
        star_boot
      ),
    by = c("target_cancer", "method")
  ) %>%
  mutate(
    star = case_when(
      method %in% c("Std CoxPH (Target)", "Std CoxPH (Combined)") ~ sig_star_from_p(pvalue),
      method %in% c("TransCox", "DiscreteKL") ~ ifelse(is.na(star_boot), "", star_boot),
      TRUE ~ ""
    ),
    label = fmt_beta_ci_with_star(beta, ci_lo, ci_hi, star)
  )

# =========================================================
# 2) SieveTL
# =========================================================
loris_tl <- NULL
if (!is.null(sievetl_df) && all(c("covariate", "target_cancer") %in% names(sievetl_df))) {
  loris_tl <- sievetl_df %>%
    filter(covariate == "LORIS") %>%
    transmute(
      target_cancer = target_cancer,
      method = "SieveTL",
      beta = as.numeric(estimate),
      ci_lo = ci_lower,
      ci_hi = ci_upper,
      pvalue = p_value
    ) %>%
    mutate(
      star = sig_star_from_p(pvalue),
      method = 'Cox-SieveTL',
      label = fmt_beta_ci_with_star(beta, ci_lo, ci_hi, star)
    )
}

# =========================================================
# combine
# =========================================================
table_long <- bind_rows(
  loris_infer,
  loris_tl
) %>%
  filter(!is.na(target_cancer), !is.na(method)) %>%
  distinct(target_cancer, method, .keep_all = TRUE)

# =========================================================
# make wide table
# =========================================================
table_wide <- table_long %>%
  mutate(
    Cancer_Type = recode(target_cancer, !!!cancer_map),
    Cancer_Type = factor(Cancer_Type, levels = cancer_map[cancer_order]),
    method = factor(method, levels = method_order)
  ) %>%
  filter(method %in% method_order) %>%
  filter(target_cancer %in% cancer_order) %>%
  arrange(Cancer_Type, method) %>%
  dplyr::select(Cancer_Type, target_cancer, method, label) %>%
  pivot_wider(names_from = method, values_from = label) %>%
  mutate(target_cancer = factor(target_cancer, levels = cancer_order)) %>%
  arrange(target_cancer) %>%
  dplyr::select(-target_cancer)

# =========================================================
# save csv
# =========================================================
write.csv(
  table_wide,
  file.path(results_dir, "table_loris_all_methods_mixed_sig.csv"),
  row.names = FALSE
)

# =========================================================
# latex-style table
# =========================================================
caption_txt <- paste(
  "LORIS across cancer types. Entries are beta [95\\% CI].",
  "For Std CoxPH and SieveTL, significance codes are based on p-values:",
  "$***\\ p<0.001$, $**\\ p<0.01$, $*\\ p<0.05$, $.\\ p<0.1$.",
  "For TransCox and DiscreteKL, significance is based on bootstrap percentile intervals:",
  "$*:$ 95\\% CI excludes 0; $**:$ 99\\% CI excludes 0; $***:$ 99.9\\% CI excludes 0."
)

kbl_obj <- kbl(
  table_wide,
  format = "latex",
  booktabs = TRUE,
  align = c("l", rep("c", ncol(table_wide) - 1)),
  caption = caption_txt,
  escape = FALSE
) %>%
  kable_styling(
    latex_options = c("hold_position", "scale_down"),
    font_size = 9
  )

writeLines(
  as.character(kbl_obj),
  con = file.path(results_dir, "table_loris_all_methods_mixed_sig.tex")
)

# =========================================================
# 3) SieveTL-only detailed table
# =========================================================
if (!is.null(sievetl_df) && all(c("covariate", "cancer") %in% names(sievetl_df))) {
  
  loris_sievetl <- sievetl_df %>%
    filter(covariate == "LORIS") %>%
    transmute(
      target_cancer = cancer,
      beta = as.numeric(estimate),
      se = se,
      ci_lo = ci_lower,
      ci_hi = ci_upper,
      pvalue = pvalue
    )
  
  tab <- loris_sievetl %>%
    filter(target_cancer %in% cancer_order) %>%
    mutate(
      `Cancer Type` = recode(target_cancer, !!!cancer_map),
      `Cancer Type` = factor(`Cancer Type`, levels = cancer_map[cancer_order]),
      star = p_star(pvalue),
      
      `Estimate (SE)` = case_when(
        is.na(beta) ~ "--",
        TRUE ~ paste0(
          "$",
          fmt_num(beta),
          ifelse(is.na(se), "", paste0("(", fmt_num(se), ")")),
          "$"
        )
      ),
      
      `95% CI` = case_when(
        is.na(ci_lo) | is.na(ci_hi) ~ "--",
        TRUE ~ paste0("[", fmt_num(ci_lo), ", ", fmt_num(ci_hi), "]")
      ),
      
      `p-value` = paste0(
        "$",
        fmt_p(pvalue),
        star,
        "$"
      )
    ) %>%
    arrange(`Cancer Type`) %>%
    select(`Cancer Type`, `Estimate (SE)`, `95% CI`, `p-value`)
  
  tab[is.na(tab)] <- "--"
  
  write_safe_tex_table <- function(df, file, caption = NULL, label = NULL) {
    df2 <- as.data.frame(df, stringsAsFactors = FALSE)
    
    for (j in seq_along(df2)) {
      df2[[j]] <- ifelse(
        str_detect(df2[[j]], "^\\$.*\\$$"),
        df2[[j]],
        latex_escape(df2[[j]])
      )
    }
    
    hdr <- latex_escape(names(df2))
    colspec <- paste0("l", paste(rep("c", ncol(df2) - 1), collapse = ""))
    
    lines <- c(
      "\\begin{table}[htbp]",
      "\\centering"
    )
    
    if (!is.null(caption)) {
      lines <- c(lines, paste0("\\caption{", caption, "}"))
    }
    if (!is.null(label)) {
      lines <- c(lines, paste0("\\label{", label, "}"))
    }
    
    lines <- c(
      lines,
      "\\resizebox{\\textwidth}{!}{%",
      paste0("\\begin{tabular}{", colspec, "}"),
      "\\toprule",
      paste0(paste(hdr, collapse = " & "), " \\\\"),
      "\\midrule"
    )
    
    row_lines <- apply(df2, 1, function(r) paste0(paste(r, collapse = " & "), " \\\\"))
    lines <- c(
      lines,
      row_lines,
      "\\bottomrule",
      "\\end{tabular}",
      "}%",
      "\\end{table}"
    )
    
    writeLines(lines, file)
  }
  
  out_tex <- file.path(results_dir, "sievetl_loris_table.tex")
  
  write_safe_tex_table(
    df = tab,
    file = out_tex,
    caption = "SieveTL results for LORIS across cancer types. Entries are Estimate(SE), 95\\% CI, and p-value. Significance levels: $^{*}p<0.05$, $^{**}p<0.01$, $^{***}p<0.001$.",
    label = "tab:sievetl_loris"
  )
  
  write.csv(
    tab,
    file.path(results_dir, "sievetl_loris_table.csv"),
    row.names = FALSE
  )
}

cat("Saved:\n")
cat(" - ", file.path(results_dir, "table_loris_all_methods_mixed_sig.csv"), "\n", sep = "")
cat(" - ", file.path(results_dir, "table_loris_all_methods_mixed_sig.tex"), "\n", sep = "")
if (!is.null(sievetl_df)) {
  cat(" - ", file.path(results_dir, "sievetl_loris_table.csv"), "\n", sep = "")
  cat(" - ", file.path(results_dir, "sievetl_loris_table.tex"), "\n", sep = "")
}

