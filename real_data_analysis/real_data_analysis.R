# ---------------------------------------------------------------------------
# Script for real data analysis of per-protocol PROVIDE data
# ---------------------------------------------------------------------------

library(glmnet)
library(SuperLearner)

here::i_am("real_data_analysis/real_data_analysis.R")

devtools::load_all("../shigella_projects/packages/vegrowth/")

data <- readRDS("provide_data/per_protocol_data.Rds")

# get rid of one obs with missing wk10 haz
data <- data[!is.na(data$wk10_haz),]

# List of covariates to adjust for
covariate_list <- c("gender",
                    "wk10_haz",
                    "num_hh_lt_5",
                    "num_hh_sleep",
                    "fedu_bin",
                    "medu_bin",
                    "inco",
                    "elec",
                    "gas",
                    "tv",
                    "toil",
                    "food_avail")

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

# results <- vegrowth::vegrowth(data = one_hot_data, 
#                               Y_name = "any_abx_wk52",
#                               Z_name = "rotaarm",
#                               S_name = "rotaepi",
#                               X_name = one_hot_covariate_colnames,
#                               estimand = c("nat_inf", "doomed", "pop"),
#                               method = c("gcomp","aipw", "bound", "sens", "ipw"),
#                               n_boot = 1000,
#                               seed = 54321,
#                               return_se = TRUE,
#                               ml = TRUE,
#                               permutation = TRUE,
#                               Y_Z_X_library = c("SL.glm", "SL.ranger", "SL.earth", "SL.step.forward"),
#                               Y_X_library = c("SL.glm", "SL.ranger", "SL.earth", "SL.step.forward"), 
#                               S_X_library = c("SL.glm", "SL.ranger", "SL.earth", "SL.step.forward"),
#                               S_Z_X_library = c("SL.glm", "SL.ranger", "SL.earth", "SL.step.forward"),
#                               Y_X_S1_model = "any_abx_wk52 ~ 1",
#                               Y_Z_X_model = "any_abx_wk52 ~ rotaarm + wk10_haz + num_hh_lt_5 + num_hh_sleep +
#                                             inco + elec + gas + tv + toil_Septic_tank_or_toilet +
#                                             toil_Water_sealed_or_slap_latrine + toil_Pit_latrine + toil_Open_latrine +
#                                             toil_Hanging_latrine + food_avail_Deficit_whole_year + food_avail_Sometimes_deficit +
#                                             food_avail_Neither_deficit_nor_surplus + food_avail_Surplus",
#                               family = "binomial",
#                               return_models = FALSE,
#                               effect_dir = "negative")
# 
# results_fewer_cov <- vegrowth::vegrowth(data = one_hot_data, 
#                               Y_name = "any_abx_wk52",
#                               Z_name = "rotaarm",
#                               S_name = "rotaepi",
#                               X_name = c("wk10_haz", "num_hh_lt_5", "gender_Female", "gender_Male"),
#                               estimand = c("nat_inf", "doomed", "pop"),
#                               method = c("gcomp","aipw", "bound", "sens", "ipw"),
#                               n_boot = 1000,
#                               seed = 54321,
#                               return_se = TRUE,
#                               ml = TRUE,
#                               permutation = TRUE,
#                               Y_Z_X_library = c("SL.glm", "SL.ranger", "SL.earth", "SL.step.forward"),
#                               Y_X_library = c("SL.glm", "SL.ranger", "SL.earth", "SL.step.forward"), 
#                               S_X_library = c("SL.glm", "SL.ranger", "SL.earth", "SL.step.forward"),
#                               S_Z_X_library = c("SL.glm", "SL.ranger", "SL.earth", "SL.step.forward"),
#                               Y_X_S1_model = "any_abx_wk52 ~ 1",
#                               # Y_Z_X_model = "any_abx_wk52 ~ rotaarm + wk10_haz + num_hh_lt_5 + num_hh_sleep +
#                               #               inco + elec + gas + tv + toil_Septic_tank_or_toilet +
#                               #               toil_Water_sealed_or_slap_latrine + toil_Pit_latrine + toil_Open_latrine +
#                               #               toil_Hanging_latrine + food_avail_Deficit_whole_year + food_avail_Sometimes_deficit +
#                               #               food_avail_Neither_deficit_nor_surplus + food_avail_Surplus",
#                               family = "binomial",
#                               return_models = FALSE,
#                               effect_dir = "negative")

# added to principal strata paper
results_fewer_cov_same_lib <- vegrowth::vegrowth(data = one_hot_data, 
                                        Y_name = "any_abx_wk52",
                                        Z_name = "rotaarm",
                                        S_name = "rotaepi",
                                        X_name = c("wk10_haz", "num_hh_lt_5", "gender_Female", "gender_Male"),
                                        estimand = c("nat_inf", "doomed", "pop"),
                                        method = c("gcomp","aipw", "bound", "sens", "ipw"),
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
                                        # Y_Z_X_model = "any_abx_wk52 ~ rotaarm + wk10_haz + num_hh_lt_5 + num_hh_sleep +
                                        #               inco + elec + gas + tv + toil_Septic_tank_or_toilet +
                                        #               toil_Water_sealed_or_slap_latrine + toil_Pit_latrine + toil_Open_latrine +
                                        #               toil_Hanging_latrine + food_avail_Deficit_whole_year + food_avail_Sometimes_deficit +
                                        #               food_avail_Neither_deficit_nor_surplus + food_avail_Surplus",
                                        family = "binomial",
                                        return_models = FALSE,
                                        effect_dir = "negative")


saveRDS(results_fewer_cov_same_lib, here::here("real_data_analysis/results/results_fewer_cov_same_lib.Rds"))

# # these results in ppt for ACIC 
# results_abx_wk52 <- vegrowth(data = one_hot_data_rm_miss_52, 
#                              G_name = "any_abx_wk52",
#                              V_name = "rotaarm",
#                              X_name = one_hot_covariate_colnames,
#                              Y_name = "any_rotads_elisa_wk52",
#                              ml = TRUE,
#                              return_se = TRUE, 
#                              Y_X_library = c("SL.glm", "SL.ranger", "SL.earth", "SL.step.forward"),
#                              G_X_library = c("SL.glm", "SL.ranger", "SL.earth",  "SL.step.forward"),
#                              G_V_X_library = c("SL.glm", "SL.ranger", "SL.earth", "SL.step.forward"),
#                              G_V_X_model = "any_abx_wk52 ~ rotaarm + wk10_haz + num_hh_lt_5 + num_hh_sleep +
#                                             inco + elec + gas + tv + toil_Septic_tank_or_toilet +
#                                             toil_Water_sealed_or_slap_latrine + toil_Pit_latrine + toil_Open_latrine +
#                                             toil_Hanging_latrine + food_avail_Deficit_whole_year + food_avail_Sometimes_deficit +
#                                             food_avail_Neither_deficit_nor_surplus + food_avail_Surplus",
#                              G_X_Y1_model = "any_abx_wk52 ~ 1",
#                              est = c("gcomp", "gcomp_pop_estimand",
#                                      "efficient_aipw", "efficient_tmle",
#                                      "choplump", "hudgens_lower", "hudgens_upper",
#                                      "hudgens_upper_doomed", "hudgens_lower_doomed"),
#                              family = "binomial",
#                              effect_dir = "negative")
# 
# 
# # look at "_new" fuctions for pt estimands of bounds updated (need to finish integrating into the main vegrowth fn)
# 
# nat_inf_bounds <- get_hudgens_stat_new(data = one_hot_data_rm_miss_52, 
#                                        G_name = "any_abx_wk52", 
#                                        V_name = "rotaarm", 
#                                        Y_name = "any_rotads_elisa_wk52", 
#                                        family = "binomial")
# 
# # additive: -0.42 lower, 0.08 upper
# 
# doomed_bounds <- get_hudgens_stat_doomed_new(data = one_hot_data_rm_miss_52, 
#                                              G_name = "any_abx_wk52", 
#                                              V_name = "rotaarm", 
#                                              Y_name = "any_rotads_elisa_wk52", 
#                                              family = "binomial")
# 
# # additive: -0.04 lower, 0.15 upper
# 
# # ----------------------------------------------------------------------------
# # NEW: Sensitivity analysis function
# # ----------------------------------------------------------------------------
# 
# models <- fit_models(
#   data = one_hot_data_rm_miss_52,
#   G_name = "any_abx_wk52",
#   V_name = "rotaarm",
#   X_name = one_hot_covariate_colnames,
#   Y_name = "any_rotads_elisa_wk52",
#   est = "efficient_aipw"
# )
# 
# out <- do_sens_aipw(
#   one_hot_data_rm_miss_52,
#   models,
#   G_name = "any_abx_wk52",
#   V_name = "rotaarm",
#   Y_name = "rotaepi",
#   return_se = TRUE
# )
# 
# plot.sens(out, se = TRUE, effect_type = "additive")
# plot.sens(out, se = TRUE, effect_type = "multiplicative")
# 
