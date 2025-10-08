here::i_am("plot_performance_generic.R")

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork) # for combining plots

cargs <- commandArgs(trailingOnly = TRUE)
setting <- cargs[[1]]

cfg <- yaml::read_yaml("config_generic.yml")
config <- cfg[[setting]]

df <- readRDS(paste0("results/generic_ER/", setting, "_summary.Rds")) %>%
 filter(method != "ipw_ER")

# Pivot additive columns only
df_add <- df %>%
  select(estimand, method, n, bias_additive, var_additive, coverage_additive) %>%
  rename(
    bias = bias_additive,
    var = var_additive,
    coverage = coverage_additive
  ) %>%
  mutate(
    estimand = case_when(
      estimand == "doomed" ~ "Doomed",
      estimand == "nat_inf" ~ "Naturally infected",
      estimand == "pop" ~ "Population",
      TRUE ~ estimand  
    ),
    method = case_when(
      method == "aipw" ~ "AIPW",
      method == "tmle" ~ "TMLE",
      method == "ipw" ~ "IPW",
      method == "gcomp" ~ "G-computation",
      method == "aipw_ER" ~ "AIPW (ER)",
      method == "ipw_ER" ~ "IPW (ER)",
      method == "gcomp_ER" ~ "G-computation (ER)",
      TRUE ~ method  
    )
  )

# Function to create a performance plot
plot_performance <- function(df, y_var, y_label) {
  ggplot(df, aes(x = n, y = .data[[y_var]], color = method, linetype = method, shape = method)) +
    geom_line(size = 1, alpha = 0.9) +
    geom_point(size = 2.5, alpha = 0.9) +
    facet_wrap(~estimand, scales = "free_y") +
    scale_x_log10(breaks = unique(df$n)) + 
    labs(x = "Sample size (n)", y = y_label, color = "Method", linetype = "Method", shape = "Method") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
}

bias_plot <- plot_performance(df_add, "bias", "Bias")
var_plot <- plot_performance(df_add, "var", "Variance")
coverage_plot <- plot_performance(df_add, "coverage", "Coverage")

# Combine plots vertically with shared legend at the bottom
combined_plot <- bias_plot / var_plot / coverage_plot + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Show plot
combined_plot

ggsave(
  filename = here::here(paste0("results/generic_ER/figures/", setting, "_results_additive.png")),
  plot = combined_plot,
  width = 13.5,  
  height = 9, 
  dpi = 300
)
