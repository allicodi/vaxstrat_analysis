# ---------------------------------------------------------------------------
# Script for real data analysis of per-protocol PROVIDE data
# ---------------------------------------------------------------------------

library(glmnet)
library(SuperLearner)
library(dplyr)
library(knitr)
library(kableExtra)

here::i_am("real_data_analysis/real_data_analysis.R")

devtools::load_all("../shigella_projects/packages/vaxstrat/")
#library(vaxstrat)

data <- readRDS("real_data_analysis/provide_data/per_protocol_data.Rds")

# get rid of one obs with missing wk10 haz
data <- data[!is.na(data$wk10_haz),]

covariate_list <- c("gender",
                    "wk10_haz",
                    "num_hh_sleep") 

# One-hot encode factors
one_hot_data <- fastDummies::dummy_cols(data,
                                        select_columns = colnames(data)[
                                          !(colnames(data) %in% c("sid", "rotaarm", "enr_date")) &
                                            (sapply(data, is.character) | sapply(data, is.factor))
                                        ],
                                        remove_first_dummy = FALSE,
                                        remove_selected_columns = TRUE,
                                        ignore_na = TRUE)

# Remove spaces from column names in the one-hot encoded data
colnames(one_hot_data) <- gsub(" ", "_", colnames(one_hot_data))

### Get new column names for covariates (covariate_list will have changed if any are factors)
covariate_data <- data[,covariate_list, drop = FALSE]
if(any(sapply(covariate_data, is.factor) == TRUE)){
  one_hot_covariate_data <- fastDummies::dummy_cols(covariate_data,
                                                    remove_first_dummy = FALSE,
                                                    remove_selected_columns = TRUE,
                                                    ignore_na = TRUE)
  colnames(one_hot_covariate_data) <- gsub(" ", "_", colnames(one_hot_covariate_data)) # remove spaces here as well
  one_hot_covariate_colnames <- colnames(one_hot_covariate_data)
} else {
  one_hot_covariate_colnames <- colnames(covariate_data)
}

# main results
results <- vaxstrat(data = one_hot_data, 
                   Y_name = "any_abx_wk52",
                   Z_name = "rotaarm",
                   S_name = "rotaepi",
                   X_name = one_hot_covariate_colnames,
                   estimand = c("nat_inf", "doomed", "pop"),
                   method = c("aipw", "bound", "sens"),
                   exclusion_restriction = c(TRUE, FALSE), 
                   cross_world = c(TRUE, FALSE),
                   n_boot = 1000,
                   seed = 54321,
                   return_se = TRUE,
                   ml = TRUE, 
                   permutation = TRUE,
                   Y_Z_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                   Y_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                   S_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                   S_Z_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                   Y_X_S1_model = "any_abx_wk52 ~ 1",
                   family = "binomial",
                   return_models = FALSE,
                   effect_dir = "negative",
                   epsilon = exp(seq(log(0.55), log(2.2), by = -log(0.55) / 10)))

saveRDS(results, here::here("real_data_analysis/results/main_results_ML_vaxstrat.Rds"))

### Covariate adjusted bounds (do not use models)

data$gender <- as.numeric(data$gender == "Female")

# covariate adjusted bounds: gender
results_gender_bound <- vaxstrat(data = data, 
                                  Y_name = "any_abx_wk52",
                                  Z_name = "rotaarm",
                                  S_name = "rotaepi",
                                  X_name = "gender",
                                  estimand = c("nat_inf"),
                                  method = c("cov_adj_bound"),
                                  n_boot = 1000,
                                  seed = 54321,
                                  return_se = TRUE,
                                  ml = FALSE,
                                  permutation = TRUE,
                                  Y_Z_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                                  Y_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                                  S_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                                  S_Z_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                                  Y_X_S1_model = "any_abx_wk52 ~ 1",
                                  family = "binomial",
                                  return_models = FALSE,
                                  effect_dir = "negative")

# covariate adjusted bounds: enr haz (binary)
data$enr_haz_bin <- ifelse(data$enr_haz < -1, 1, 0)

results_enr_haz_bound <- vaxstrat(data = data, 
                                 Y_name = "any_abx_wk52",
                                 Z_name = "rotaarm",
                                 S_name = "rotaepi",
                                 X_name = "enr_haz_bin",
                                 estimand = c("nat_inf"),
                                 method = c("cov_adj_bound"),
                                 n_boot = 1000,
                                 seed = 54321,
                                 return_se = TRUE,
                                 ml = FALSE,
                                 permutation = TRUE,
                                 Y_Z_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                                 Y_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                                 S_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                                 S_Z_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                                 Y_X_S1_model = "any_abx_wk52 ~ 1",
                                 family = "binomial",
                                 return_models = FALSE,
                                 effect_dir = "negative")


data$num_hh_sleep_bin <- ifelse(data$num_hh_sleep < 5, 1, 0)
results_num_hh_sleep_bound <- vaxstrat(data = data, 
                                  Y_name = "any_abx_wk52",
                                  Z_name = "rotaarm",
                                  S_name = "rotaepi",
                                  X_name = "num_hh_sleep_bin",
                                  estimand = c("nat_inf"),
                                  method = c("cov_adj_bound"),
                                  n_boot = 1000,
                                  seed = 54321,
                                  return_se = TRUE,
                                  ml = FALSE,
                                  permutation = TRUE,
                                  Y_Z_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                                  Y_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                                  S_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                                  S_Z_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                                  Y_X_S1_model = "any_abx_wk52 ~ 1",
                                  family = "binomial",
                                  return_models = FALSE,
                                  effect_dir = "negative")


# gender x enr_haz
data$I_gender_x_enr_haz <- data$gender * data$enr_haz_bin 
results_I_gender_x_enr_haz_bound <- vaxstrat(data = data, 
                                       Y_name = "any_abx_wk52",
                                       Z_name = "rotaarm",
                                       S_name = "rotaepi",
                                       X_name = "I_gender_x_enr_haz",
                                       estimand = c("nat_inf"),
                                       method = c("cov_adj_bound"),
                                       n_boot = 1000,
                                       seed = 54321,
                                       return_se = TRUE,
                                       ml = FALSE,
                                       permutation = TRUE,
                                       Y_Z_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                                       Y_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                                       S_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                                       S_Z_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                                       Y_X_S1_model = "any_abx_wk52 ~ 1",
                                       family = "binomial",
                                       return_models = FALSE,
                                       effect_dir = "negative")

# gender x num sleep
data$I_gender_x_num_hh_sleep <- data$gender * data$num_hh_sleep_bin
results_I_gender_x_num_hh_sleep_bound <- vaxstrat(data = data, 
                                             Y_name = "any_abx_wk52",
                                             Z_name = "rotaarm",
                                             S_name = "rotaepi",
                                             X_name = "I_gender_x_num_hh_sleep",
                                             estimand = c("nat_inf"),
                                             method = c("cov_adj_bound"),
                                             n_boot = 1000,
                                             seed = 54321,
                                             return_se = TRUE,
                                             ml = FALSE,
                                             permutation = TRUE,
                                             Y_Z_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                                             Y_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                                             S_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                                             S_Z_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                                             Y_X_S1_model = "any_abx_wk52 ~ 1",
                                             family = "binomial",
                                             return_models = FALSE,
                                             effect_dir = "negative")

# enr_haz x num_sleep
data$I_enr_haz_x_num_hh_sleep <- data$enr_haz_bin * data$num_hh_sleep_bin
results_I_enr_haz_x_num_hh_sleep_bound <- vaxstrat(data = data, 
                                                  Y_name = "any_abx_wk52",
                                                  Z_name = "rotaarm",
                                                  S_name = "rotaepi",
                                                  X_name = "I_enr_haz_x_num_hh_sleep",
                                                  estimand = c("nat_inf"),
                                                  method = c("cov_adj_bound"),
                                                  n_boot = 1000,
                                                  seed = 54321,
                                                  return_se = TRUE,
                                                  ml = FALSE,
                                                  permutation = TRUE,
                                                  Y_Z_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                                                  Y_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                                                  S_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                                                  S_Z_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                                                  Y_X_S1_model = "any_abx_wk52 ~ 1",
                                                  family = "binomial",
                                                  return_models = FALSE,
                                                  effect_dir = "negative")

data$I_gender_x_enr_haz_x_num_hh_sleep <- data$gender * data$enr_haz_bin * data$num_hh_sleep_bin
results_I_all_bound <- vaxstrat(data = data, 
                                 Y_name = "any_abx_wk52",
                                 Z_name = "rotaarm",
                                 S_name = "rotaepi",
                                 X_name = "I_gender_x_enr_haz_x_num_hh_sleep",
                                 estimand = c("nat_inf"),
                                 method = c("cov_adj_bound"),
                                 n_boot = 1000,
                                 seed = 54321,
                                 return_se = TRUE,
                                 ml = FALSE,
                                 permutation = TRUE,
                                 Y_Z_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                                 Y_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                                 S_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                                 S_Z_X_library = c("SL.glm", "SL.gam", "SL.earth", "SL.step.forward"),
                                 Y_X_S1_model = "any_abx_wk52 ~ 1",
                                 family = "binomial",
                                 return_models = FALSE,
                                 effect_dir = "negative")


# ------------------------------------------------------------------------

# Create Latex tables with results

# Primary table: AIPW results for each estimand

table_1_additive <- data.frame(estimand = c(rep("Naturally Infected", 3),
                                             "Doomed",
                                             "Marginal"),
                                estimator = c("AIPW CW", "AIPW ER", "AIPW CW + ER",
                                              "AIPW", 
                                              "AIPW"),
                                estimate = c(results$nat_inf$aipw_CW$pt_est['additive_effect'],
                                             results$nat_inf$aipw_ER$pt_est['additive_effect'],
                                             results$nat_inf$aipw_ER_CW$pt_est['additive_effect'],
                                             results$doomed$aipw$pt_est['additive_effect'],
                                             results$pop$aipw$pt_est['additive_effect']),
                                lci = c(results$nat_inf$aipw_CW$pt_est['additive_effect'] - 1.96*results$nat_inf$aipw_CW$pt_est['additive_se'],
                                        results$nat_inf$aipw_ER$pt_est['additive_effect'] - 1.96*results$nat_inf$aipw_ER$pt_est['additive_se'],
                                        results$nat_inf$aipw_ER_CW$pt_est['additive_effect'] - 1.96*results$nat_inf$aipw_ER_CW$pt_est['additive_se'],
                                        results$doomed$aipw$pt_est['additive_effect'] - 1.96*results$doomed$aipw$pt_est['additive_se'],
                                        results$pop$aipw$pt_est['additive_effect'] - 1.96*results$pop$aipw$pt_est['additive_se']),
                                uci = c(results$nat_inf$aipw_CW$pt_est['additive_effect'] + 1.96*results$nat_inf$aipw_CW$pt_est['additive_se'],
                                        results$nat_inf$aipw_ER$pt_est['additive_effect'] + 1.96*results$nat_inf$aipw_ER$pt_est['additive_se'],
                                        results$nat_inf$aipw_ER_CW$pt_est['additive_effect'] + 1.96*results$nat_inf$aipw_ER_CW$pt_est['additive_se'],
                                        results$doomed$aipw$pt_est['additive_effect'] + 1.96*results$doomed$aipw$pt_est['additive_se'],
                                        results$pop$aipw$pt_est['additive_effect'] + 1.96*results$pop$aipw$pt_est['additive_se']),
                                p_val = c(results$nat_inf$aipw_CW$p_val$additive,
                                          results$nat_inf$aipw_ER$p_val$additive,
                                          results$nat_inf$aipw_ER_CW$p_val$additive,
                                          results$doomed$aipw$p_val$additive,
                                          results$pop$aipw$p_val$additive))

table_1_mult <- data.frame(estimand = c(rep("Naturally Infected", 3),
                                            "Doomed",
                                            "Marginal"),
                               estimator = c("AIPW CW", "AIPW ER", "AIPW CW + ER",
                                             "AIPW", 
                                             "AIPW"),
                               estimate = exp(c(results$nat_inf$aipw_CW$pt_est['log_multiplicative_effect'],
                                            results$nat_inf$aipw_ER$pt_est['log_multiplicative_effect'],
                                            results$nat_inf$aipw_ER_CW$pt_est['log_multiplicative_effect'],
                                            results$doomed$aipw$pt_est['log_multiplicative_effect'],
                                            results$pop$aipw$pt_est['log_multiplicative_effect'])),
                               lci = exp(c(results$nat_inf$aipw_CW$pt_est['log_multiplicative_effect'] - 1.96*results$nat_inf$aipw_CW$pt_est['log_multiplicative_se'],
                                       results$nat_inf$aipw_ER$pt_est['log_multiplicative_effect'] - 1.96*results$nat_inf$aipw_ER$pt_est['log_multiplicative_se'],
                                       results$nat_inf$aipw_ER_CW$pt_est['log_multiplicative_effect'] - 1.96*results$nat_inf$aipw_ER_CW$pt_est['log_multiplicative_se'],
                                       results$doomed$aipw$pt_est['log_multiplicative_effect'] - 1.96*results$doomed$aipw$pt_est['log_multiplicative_se'],
                                       results$pop$aipw$pt_est['log_multiplicative_effect'] - 1.96*results$pop$aipw$pt_est['log_multiplicative_se'])),
                               uci = exp(c(results$nat_inf$aipw_CW$pt_est['log_multiplicative_effect'] + 1.96*results$nat_inf$aipw_CW$pt_est['log_multiplicative_se'],
                                       results$nat_inf$aipw_ER$pt_est['log_multiplicative_effect'] + 1.96*results$nat_inf$aipw_ER$pt_est['log_multiplicative_se'],
                                       results$nat_inf$aipw_ER_CW$pt_est['log_multiplicative_effect'] + 1.96*results$nat_inf$aipw_ER_CW$pt_est['log_multiplicative_se'],
                                       results$doomed$aipw$pt_est['log_multiplicative_effect'] + 1.96*results$doomed$aipw$pt_est['log_multiplicative_se'],
                                       results$pop$aipw$pt_est['log_multiplicative_effect'] + 1.96*results$pop$aipw$pt_est['log_multiplicative_se'])),
                               p_val = c(results$nat_inf$aipw_CW$p_val$mult,
                                         results$nat_inf$aipw_ER$p_val$mult,
                                         results$nat_inf$aipw_ER_CW$p_val$mult,
                                         results$doomed$aipw$p_val$mult,
                                         results$pop$aipw$p_val$mult))
                      
table_additive <- table_1_additive %>%
  mutate(
    Scale = "Additive"
  )

table_multiplicative <- table_1_mult %>%
  mutate(
    Scale = "Multiplicative"
  )

combined_table <- bind_rows(table_additive, table_multiplicative) %>%
  mutate(
    Scale = factor(Scale, levels = c("Additive", "Multiplicative"))
  ) %>%
  arrange(estimand, estimator, Scale)

table_wide <- combined_table %>%
  select(estimand, estimator, Scale, estimate, lci, uci, p_val) %>%
  tidyr::pivot_wider(
    names_from = Scale,
    values_from = c(estimate, lci, uci, p_val)
  ) %>%
  arrange(
    factor(estimand, levels = c("Marginal", "Doomed", "Naturally Infected")),
    factor(estimator, levels = c("AIPW CW", "AIPW ER", "AIPW CW + ER", "AIPW"))
  ) %>%
  select(estimand, estimator, estimate_Additive, lci_Additive, uci_Additive, p_val_Additive, estimate_Multiplicative, lci_Multiplicative, uci_Multiplicative, p_val_Multiplicative) %>%
  rename(
    `Estimate` = estimate_Additive,
    `LCI` = lci_Additive,
    `UCI` = uci_Additive,
    `p-value` = p_val_Additive,
    `Estimate ` = estimate_Multiplicative,
    `LCI ` = lci_Multiplicative,
    `UCI ` = uci_Multiplicative,
    `p-value ` = p_val_Multiplicative
  )

table_wide_fmt <- table_wide %>%
  mutate(
    `Additive` = sprintf(
      "%.3f (%.3f, %.3f)",
      Estimate, LCI, UCI
    ),
    `Multiplicative` = sprintf(
      "%.3f (%.3f, %.3f)",
      `Estimate `, `LCI `, `UCI `
    ),
    p_val_Additive = round(`p-value`, 3),
    p_val_Multiplicative = round(`p-value `, 3)
  ) %>%
  select(
    estimand,
    estimator,
    Additive,
    p_val_Additive = `p-value`,
    Multiplicative,
    p_val_Multiplicative = `p-value `
  )

section_sizes <- table_wide %>%
  group_by(estimand) %>%
  summarise(n = n()) %>%
  pull(n)

kable(
  table_wide_fmt,
  format = "latex",
  booktabs = TRUE,
  caption = "AIPW Estimates on Additive and Multiplicative Scales",
  label = "tab:aipw_results",
  align = "llcccc"
) %>%
  kable_styling(latex_options = "hold_position") %>%
  add_header_above(c(
    " " = 2,
    "Additive Scale" = 2,
    "Multiplicative Scale" = 2
  )) %>%
  add_header_above(c(
    " " = 2,
    "Estimate (95\\% CI)" = 1,
    "p-value" = 1,
    "Estimate (95\\% CI)" = 1,
    "p-value" = 1
  )) %>%
  collapse_rows(
    columns = 1,
    latex_hline = "major",
    valign = "top"
  )
                      
# Secondary table: Bound results (regular, covariate adjusted)

# -----------------------------------------------------------------------------
# Helper: format bootstrap percentile CI for a given stratum/result object
# -----------------------------------------------------------------------------
fmt_bound_boot <- function(res_obj, scale = c("additive", "mult"), bound = c("lower", "upper"), digits = 3, cov_adj = FALSE) {
  scale <- match.arg(scale)
  bound <- match.arg(bound)
  
  # boot table and names (based on your example)
  if(!cov_adj){
    boot <- res_obj$nat_inf$bound$boot_se
  } else{
    boot <- res_obj$nat_inf$cov_adj_bound$boot_se
  }
  
  
  # estimate names in result object
  est_name <- switch(
    paste0(scale, "_", bound),
    "additive_lower" = "additive_effect_lower",
    "additive_upper" = "additive_effect_upper",
    "mult_lower"     = "mult_effect_lower",
    "mult_upper"     = "mult_effect_upper"
  )
  
  if(!cov_adj){
    est <- res_obj$nat_inf$bound$pt_est[[est_name]]
  } else{
    est <- res_obj$nat_inf$cov_adj_bound$pt_est[[est_name]]
  }
  
  
  # bootstrap CI column names in boot_se (from your example)
  ci_lower_col <- switch(
    paste0(scale, "_", bound),
    "additive_lower" = "lower_ci_additive_lower",
    "additive_upper" = "lower_ci_additive_upper",
    "mult_lower"     = "lower_ci_mult_lower",
    "mult_upper"     = "lower_ci_mult_upper"
  )
  ci_upper_col <- switch(
    paste0(scale, "_", bound),
    "additive_lower" = "upper_ci_additive_lower",
    "additive_upper" = "upper_ci_additive_upper",
    "mult_lower"     = "upper_ci_mult_lower",
    "mult_upper"     = "upper_ci_mult_upper"
  )
  
  ci_low <- as.numeric(boot[[ci_lower_col]])
  ci_high <- as.numeric(boot[[ci_upper_col]])
  
  sprintf(
    "%.*f (%.*f, %.*f)",
    digits, est,
    digits, ci_low,
    digits, ci_high
  )
}

# -----------------------------------------------------------------------------
# Wrapper to build a table-ready data.frame for bounds
# - unadjusted_res: the main `results` object (unadjusted bound)
# - cov_adj_list: named list of cov-adj result objects (names used as Specification labels)
# -----------------------------------------------------------------------------
build_bounds_table_for_latex <- function(unadjusted_res,
                                         cov_adj_list = list(),
                                         digits = 2) {
  
  # Unadjusted row (assumes bound stored at results$nat_inf$bound)
  unadj_row <- data.frame(
    Specification = "Unadjusted",
    `Additive Lower` = fmt_bound_boot(unadjusted_res, scale = "additive", bound = "lower", digits = digits),
    `Additive Upper` = fmt_bound_boot(unadjusted_res, scale = "additive", bound = "upper", digits = digits),
    `Multiplicative Lower` = fmt_bound_boot(unadjusted_res, scale = "mult", bound = "lower", digits = digits),
    `Multiplicative Upper` = fmt_bound_boot(unadjusted_res, scale = "mult", bound = "upper", digits = digits),
    stringsAsFactors = FALSE
  )
  
  # Covariate-adjusted rows: expect cov_adj_list to be a named list, e.g.
  # list("Gender" = results_gender_bound, "Enrollment HAZ" = results_enr_haz_bound, ...)
  cov_rows <- bind_rows(lapply(seq_along(cov_adj_list), function(i) {
    lab <- names(cov_adj_list)[i]
    resi <- cov_adj_list[[i]]
    data.frame(
      Specification = lab,
      `Additive Lower` = fmt_bound_boot(resi, scale = "additive", bound = "lower", digits = digits, cov_adj = TRUE),
      `Additive Upper` = fmt_bound_boot(resi, scale = "additive", bound = "upper", digits = digits, cov_adj = TRUE),
      `Multiplicative Lower` = fmt_bound_boot(resi, scale = "mult", bound = "lower", digits = digits, cov_adj = TRUE),
      `Multiplicative Upper` = fmt_bound_boot(resi, scale = "mult", bound = "upper", digits = digits, cov_adj = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  
  # Combine into named groups for printing
  df_unadj <- unadj_row
  df_cov   <- cov_rows
  
  # Return a list of data frames, so we can use group_rows() like in your sim table
  list(
    "Unadjusted Bound" = df_unadj,
    "Covariate-Adjusted Bound" = df_cov
  )
}


cov_adj_list <- list(
  "Gender"                  = results_gender_bound,
  "Enrollment HAZ (bin)"    = results_enr_haz_bound,
  "Household size (bin)"    = results_num_hh_sleep_bound,
  "Gender × HAZ"            = results_I_gender_x_enr_haz_bound,
  "Gender × Household"      = results_I_gender_x_num_hh_sleep_bound,
  "HAZ × Household"         = results_I_enr_haz_x_num_hh_sleep_bound,
  "All interactions"        = results_I_all_bound
)

tables_list <- build_bounds_table_for_latex(results, cov_adj_list, digits = 3)

# Bind them (preserve group labels)
combined_bounds <- bind_rows(tables_list, .id = "Method")

# -----------------------------------------------------------------------------
# Print LaTeX table for Overleaf
# -----------------------------------------------------------------------------
section_sizes <- sapply(tables_list, nrow)

kable(
  combined_bounds %>% select(-Method),
  format = "latex",
  booktabs = TRUE,
  caption = "Bootstrap percentile bounds for the naturally infected stratum",
  label = "tab:realdata_bounds",
  align = "lcccc",
  escape = FALSE
) %>%
  kable_styling(latex_options = c("hold_position"), font_size = 9) %>%
  add_header_above(c(
    " " = 1,
    "Additive" = 2,
    "Multiplicative" = 2
  )) %>%
  add_header_above(c(
    " " = 1,
    "Lower bound" = 1,
    "Upper bound" = 1,
    "Lower bound" = 1,
    "Upper bound" = 1
  )) %>%
  group_rows(
    index = section_sizes,
    group_label = names(tables_list)
  )
