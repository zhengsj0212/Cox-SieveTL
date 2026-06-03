# Beta
#################### Tables
library(dplyr)
library(tidyr)
library(readr)

results_dir <- "./results_allTMB_new_hessian_lasso_new"
beta_file <- file.path(results_dir, "fullsample_beta_estimates.csv")

beta_df <- read.csv(beta_file, stringsAsFactors = FALSE)

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

loris_table <- beta_df %>%
  filter(method %in% c("target_only", "combined_transfer", "full_method")) %>%
  transmute(
    target_cancer,
    Cancer = ifelse(target_cancer %in% names(title_map), title_map[target_cancer], target_cancer),
    method = recode(
      method,
      target_only = "Target Only",
      combined_transfer = "Combined",
      full_method = "Cox-SieveTL"
    ),
    LORIS = as.numeric(LORIS)
  ) %>%
  pivot_wider(
    names_from = method,
    values_from = LORIS
  )

# write.csv(
#   loris_table,
#   file.path(results_dir, "loris_beta_comparison_table.csv"),
#   row.names = FALSE
# )

################# Connected plots
library(dplyr)
library(tidyr)
library(ggplot2)

beta_path <- file.path(results_dir, "fullsample_beta_estimates.csv")
beta_df <- read.csv(beta_path, stringsAsFactors = FALSE)

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

plot_beta <- beta_df %>%
  filter(method %in% c("combined_transfer", "full_method")) %>%
  mutate(
    Cancer = ifelse(target_cancer %in% names(title_map), title_map[target_cancer], target_cancer),
    Method = recode(
      method,
      combined_transfer = "Combined",
      full_method = "Cox-SieveTL"
    )
  ) %>%
  pivot_longer(
    cols = -c(target_cancer, method, Cancer, Method),
    names_to = "covariate",
    values_to = "beta_value"
  )

# optional: control covariate order
covariate_levels <- unique(plot_beta$covariate)
plot_beta <- plot_beta %>%
  mutate(
    covariate = factor(covariate, levels = covariate_levels)
  )

p_beta <- ggplot(
  plot_beta,
  aes(x = covariate, y = beta_value, color = Method, shape = Method, group = Method)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  facet_wrap(~ Cancer, scales = "free_y", ncol = 5) +
  labs(
    x = "Covariate",
    y = expression(hat(beta)),
    color = "Method",
    shape = "Method"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    strip.text = element_text(size = 9)
  )

print(p_beta)

ggsave(
  file.path(results_dir, "beta_compare_combined_vs_sievetl.pdf"),
  p_beta,
  width = 12,
  height = 9
)

# Gamma
library(dplyr)
library(ggplot2)

gamma_path <- file.path(results_dir, "fullsample_gamma_estimates.csv")
gamma_df <- read.csv(gamma_path, stringsAsFactors = FALSE)

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

plot_gamma <- gamma_df %>%
  filter(method %in% c("combined_transfer", "full_method")) %>%
  mutate(
    Cancer = ifelse(target_cancer %in% names(title_map), title_map[target_cancer], target_cancer),
    Method = recode(
      method,
      combined_transfer = "Combined",
      full_method = "Cox-SieveTL"
    )
  )
plot_gamma <- plot_gamma %>%
  mutate(
    basis_term = factor(
      basis_term,
      levels = paste0("g", 1:max(as.numeric(gsub("g", "", basis_term))))
    )
  )
p_gamma <- ggplot(
  plot_gamma,
  aes(x = basis_term, y = gamma_value, color = Method, shape = Method, group = Method)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  facet_wrap(~ Cancer, scales = "free_y", ncol = 5) +
  labs(
    x = "Basis term",
    y = expression(hat(gamma)),
    color = "Method",
    shape = "Method"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    strip.text = element_text(size = 9)
  )


print(p_gamma)

ggsave(
  file.path(results_dir, "gamma_compare_combined_vs_sievetl.pdf"),
  p_gamma,
  width = 12,
  height = 9
)

##### All covariates connecte plots
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# 
# beta_path <- file.path(results_dir, "fullsample_beta_estimates.csv")
# beta_df <- read.csv(beta_path, stringsAsFactors = FALSE)
# 
# title_map <- c(
#   cancer1  = "Bladder",
#   cancer2  = "Breast",
#   cancer3  = "Colorectal",
#   cancer4  = "Endometrial",
#   cancer5  = "Esophageal",
#   cancer6  = "Gastric",
#   cancer7  = "Head & Neck",
#   cancer8  = "Hepatobiliary",
#   cancer9  = "Melanoma",
#   cancer10 = "Mesothelioma",
#   cancer11 = "NSCLC",
#   cancer12 = "Ovarian",
#   cancer13 = "Pancreatic",
#   cancer14 = "Renal",
#   cancer15 = "Sarcoma",
#   cancer16 = "SCLC",
#   cancer18 = "CNS"
# )
# 
# # ------------------------------------------------------------
# # 1. Identify covariate columns (IMPORTANT)
# # ------------------------------------------------------------
# non_cov_cols <- c("target_cancer", "method")
# 
# covariate_cols <- setdiff(names(beta_df), non_cov_cols)
# 
# # OPTIONAL: remove unwanted columns if needed
# covariate_cols <- covariate_cols[!grepl("Intercept", covariate_cols, ignore.case = TRUE)]
# 
# # ------------------------------------------------------------
# # 2. Long format
# # ------------------------------------------------------------
# plot_beta <- beta_df %>%
#   filter(method %in% c("combined_transfer", "full_method")) %>%
#   mutate(
#     Cancer = ifelse(target_cancer %in% names(title_map),
#                     title_map[target_cancer],
#                     target_cancer),
#     Method = recode(
#       method,
#       combined_transfer = "Combined",
#       full_method = "Cox-SieveTL"
#     )
#   ) %>%
#   pivot_longer(
#     cols = all_of(covariate_cols),
#     names_to = "covariate",
#     values_to = "beta_value"
#   )
# 
# # ------------------------------------------------------------
# # 3. Clean / reorder covariates
# # ------------------------------------------------------------
# 
# # Example: put key variables first if they exist
# priority_vars <- c("loris", "tmb")
# 
# covariate_levels <- c(
#   intersect(priority_vars, unique(plot_beta$covariate)),
#   setdiff(unique(plot_beta$covariate), priority_vars)
# )
# 
# plot_beta <- plot_beta %>%
#   mutate(
#     covariate = factor(covariate, levels = covariate_levels)
#   )
# 
# # ------------------------------------------------------------
# # 4. Plot
# # ------------------------------------------------------------
# p_beta <- ggplot(
#   plot_beta,
#   aes(x = covariate, y = beta_value,
#       color = Method, shape = Method, group = Method)
# ) +
#   geom_line(linewidth = 0.8) +
#   geom_point(size = 2.2) +
#   facet_wrap(~ Cancer, scales = "free_y", ncol = 5) +
#   labs(
#     x = "Covariate",
#     y = expression(hat(beta)),
#     color = "Method",
#     shape = "Method"
#   ) +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     legend.position = "bottom",
#     strip.text = element_text(size = 9)
#   )
# 
# print(p_beta)
# 
# ggsave(
#   file.path(results_dir, "beta_compare_combined_vs_sievetl.pdf"),
#   p_beta,
#   width = 12,
#   height = 9
# )
# 
# 
# 
# 
# ##
# 
# 
# # assume your dataframe is called gamma_df
# # delta = Combined - SieveTL for each basis term within each cancer
# 
# gamma_delta_df <- plot_gamma %>%
#   mutate(
#     Method = case_when(
#       method == "combined_transfer" ~ "Combined",
#       method == "full_method" ~ "SieveTL",
#       TRUE ~ Method
#     )
#   ) %>%
#   filter(Method %in% c("Combined", "SieveTL")) %>%
#   dplyr::select(
#     target_cancer,
#     Cancer,
#     basis_term,
#     Method,
#     gamma_value
#   ) %>%
#   pivot_wider(
#     id_cols = c(target_cancer, Cancer, basis_term),
#     names_from = Method,
#     values_from = gamma_value
#   ) %>%
#   mutate(
#     delta = Combined - SieveTL
#   )
# 
# gamma_delta_df
# 
