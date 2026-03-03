# ---------------------------------------------------
# Make latex table for simulation 1 results
# ---------------------------------------------------

here::i_am("make_sim_1_table.R")

library(dplyr)
library(knitr)
library(kableExtra)

# read results 
default_summary <- readRDS(here::here("results/sim_1/default_summary.Rds"))
violate_er_summary <- readRDS(here::here("results/sim_1/violate_er_summary.Rds"))
violate_cw_summary <- readRDS(here::here("results/sim_1/violate_cw_summary.Rds"))
violate_cw_er_summary <- readRDS(here::here("results/sim_1/violate_cw_er_summary.Rds"))
  
summary_list <- list(
  "Default"      = default_summary,
  "Violate ER"   = violate_er_summary,
  "Violate CW"   = violate_cw_summary,
  "Violate Both" = violate_cw_er_summary
)

format_section <- function(df) {
  df %>%
    # filter(n == n_value) %>%
    mutate(
      method = factor(
        method,
        levels = c("aipw_CW", "aipw_ER", "aipw_ER_CW"),
        labels = c("CW AIPW", "ER AIPW", "CW & ER AIPW")
      ),
      `$sqrt(n) times$ Bias` = sqrt(n) * bias_additive,
      `$n times$ Variance` = n * var_additive, 
      `$n times$ MSE` = n * mse_additive, 
      `$sqrt(n) times$ Bias ` = sqrt(n) * bias_mult,
      `$n times$ Variance `  = n * var_mult,
      `$n times$ MSE ` = n * mse_mult
    ) %>%
    arrange(method) %>%
    select(
      Method = method,
      n,
      `$sqrt(n) times$ Bias`,
      `$n times$ Variance`,
      `$n times$ MSE`,
      `Coverage` = coverage_additive,
      `$sqrt(n) times$ Bias `,
      `$n times$ Variance `,
      `$n times$ MSE `,
      `Coverage ` = coverage_mult
    )
}


make_latex_table <- function(summary_list, caption, label) {
  
  section_tables <- lapply(summary_list, format_section)
  
  combined <- bind_rows(section_tables, .id = "Scenario")
  
  kable(
    combined %>% select(-Scenario),
    format = "latex",
    booktabs = TRUE,
    digits = 3,
    caption = caption,
    label = label,
    align = "lcccccc"
  ) %>%
    kable_styling(latex_options = c("hold_position")) %>%
    add_header_above(c(
      " " = 1,
      "Additive Scale" = 3,
      "Multiplicative Scale" = 3
    )) %>%
    group_rows(
      index = sapply(section_tables, nrow),
      group_label = names(section_tables)
    )
}

make_latex_table(
  summary_list,
  caption = "Bias, Variance, and Coverage",
  label = "tab:bias_var_cov"
)
