# ---------------------------------------------------
# Make latex table for simulation 1 results
# ---------------------------------------------------

here::i_am("make_sim_1_table_bounds.R")

library(dplyr)
library(knitr)
library(kableExtra)

# make list of summary files
prefix_list <- c("cw_er_", "violate_cw_", "violate_er_", "violate_cw_er_")
suffix_list <- c("X1", "X2", "X3", "X1X2", "X1X3", "X2X3", "X1X2X3")
combos <- expand.grid(prefix = prefix_list,
                      suffix = suffix_list)

combos$setting <- paste0(combos$prefix, combos$suffix)
combos$file_name <- paste0(combos$prefix, combos$suffix, "_summary.Rds")
combos$file_name[combos$prefix == "cw_er_" & combos$suffix == "X1"] <- "default_summary.Rds"

results_df <- data.frame()

for(i in 1:length(combos$file_name)){
  fname <- combos$file_name[i]
  setting <- combos$setting[i]
  assumption <- combos$prefix[i]
  covariates <- combos$suffix[i]
  
  res <- readRDS(here::here("results/sim_1/sim_1_bounds", fname))  
  res$setting <- setting
  res$assumption <- assumption
  res$X <- covariates
  
  results_df <- rbind(results_df, res)
  
}

bounds <- results_df[results_df$method == "bound",]
cov_adj_bound <- results_df[results_df$method == "cov_adj_bound",]

format_section <- function(df, bound_type, n_value) {
  
  df_out <- df %>%
    filter(n == n_value, method == bound_type) %>%
    mutate(
      assumption = factor(
        assumption,
        levels = c("cw_er_", "violate_cw_", "violate_er_", "violate_cw_er_"),
        labels = c("CW + ER", "Violate CW", "Violate ER", "Violate CW + ER")
      ),
      `Bias: Lower` = sqrt(n) * bias_lower_additive,
      `Bias: Upper` = sqrt(n) * bias_upper_additive,
      `Bound Width` = paste0(
        round(median_width_additive, 2),
        " (",
        round(q25_width_additive, 2),
        ", ",
        round(q75_width_additive, 2),
        ")"
      )
    )
  
  if (bound_type == "bound") {
    df_out %>%
      arrange(assumption) %>%
      mutate(Covariates = "") %>%
      select(
        Covariates,
        Assumption = assumption,
        `Bias: Lower`,
        `Bias: Upper`,
        `Lower Coverage` = coverage_lower_bound_additive,
        `Upper Coverage` = coverage_upper_bound_additive,
        `Point Estimate Coverage` = coverage_pt_est_additive,
        `Bound Width`
      )
  } else {
    df_out %>%
      arrange(X, assumption) %>%
      select(
        Covariates = X,
        Assumption = assumption,
        `Bias: Lower`,
        `Bias: Upper`,
        `Lower Coverage` = coverage_lower_bound_additive,
        `Upper Coverage` = coverage_upper_bound_additive,
        `Point Estimate Coverage` = coverage_pt_est_additive,
        `Bound Width`
      )
  }
}

bound_500 <- format_section(results_df, "bound", 500)
cov_adj_bound_500 <- format_section(results_df, "cov_adj_bound", 500)

bound_4000 <- format_section(results_df, "bound", 4000)
cov_adj_bound_4000 <- format_section(results_df, "cov_adj_bound", 4000)

tables_500 <- list(
  "Unadjusted Bound" = bound_500,
  "Covariate-Adjusted Bound" = cov_adj_bound_500
)

tables_4000 <- list(
  "Unadjusted Bound" = bound_4000,
  "Covariate-Adjusted Bound" = cov_adj_bound_4000
)

make_latex_table <- function(tables, caption, label) {
  
  combined <- bind_rows(tables, .id = "Method")
  
  kable(
    combined %>% select(-Method),
    format = "latex",
    booktabs = TRUE,
    caption = caption,
    label = label,
    align = "llcccccc",
    digits = 3,
    escape = FALSE
  ) %>%
    kable_styling(
      latex_options = c("hold_position"),
      font_size = 9
    ) %>%
    add_header_above(c(
      " " = 2,
      "$\\sqrt{n} \\times \\text{Bias}$" = 2,
      "Coverage" = 3,
      " " = 1
    )) %>%
    group_rows(
      index = sapply(tables, nrow),
      group_label = names(tables)
    )
}

make_latex_table(
  tables_500,
  caption = "Bias, Coverage, and Bound Width (n = 500; Additive)",
  label = "tab:sim1_bounds_n500"
)

make_latex_table(
  tables_4000,
  caption = "Bias, Coverage, and Bound Width (n = 4000; Additive)",
  label = "tab:sim1_bounds_n4000"
)
