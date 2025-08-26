# ------------------------------------------------------
# Script to run simulation on cluster
# ------------------------------------------------------

here::i_am("run_analysis.R")

source(here::here("simulate_data.R"))

#devtools::load_all("../shigella_projects/packages/vegrowth/")
library(vegrowth)
library(future)
library(future.apply)
library(SuperLearner)

# For initial debugging scratch file
options(echo = TRUE)

# this was for protect, probably can make smaller
options(future.globals.maxSize = 5 * 1024^3) # 5 GB
options(future.globals.onReference = "ignore")

# Path to installed packages on cluster
.libPaths(c("/apps/R/4.4.0/lib64/R/site/library","/apps/R/4.4.0/lib64/R/library", "~/Rlibs_ve_trial"))

ncores <- parallelly::availableCores()
future::plan("multicore", workers = ncores)

# Path to projects folder where results will be saved
project_dir <- "/projects/dbenkes/allison/ve_growth_sims/results/"
seed <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
setting <- Sys.getenv("SETTING")

# read in config file
config <- config::get(file = here::here("config.yml"), config = setting)

# Create dir to save results if does not exist
if (!file.exists(paste0(project_dir, setting))) {
  dir.create(here::here(paste0(project_dir, setting)), recursive = TRUE)
}

# get all parameter grid combinations
grid <- expand.grid(seed = seed, 
                    n_sample_size = as.numeric(config$n_sample_size),
                    effect_protect = config$effect_protect,
                    inflation = as.numeric(config$inflation),
                    nat_inf_epsilon = as.numeric(config$nat_inf_epsilon),
                    doomed_epsilon = as.numeric(config$doomed_epsilon))

results <- future.apply::future_lapply(1:nrow(grid), function(i, grid, data, seed){
  
  data <- simulate_data(seed = seed,
                        effect_protect = grid$effect_protect[i],
                        inflation = grid$inflation[i],
                        nat_inf_epsilon = grid$nat_inf_epsilon[i], 
                        doomed_epsilon = grid$doomed_epsilon[i],
                        n = grid$n_sample_size[i])

  results <- vegrowth::vegrowth(data = data, 
                                Y_name = "any_abx_wk52",
                                Z_name = "rotaarm",
                                S_name = "rotaepi",
                                X_name = c("wk10_haz", "gender", "num_hh_sleep"),
                                estimand = c("nat_inf", "doomed", "pop"),
                                method = c("aipw", "bound"),
                                n_boot = 1000,
                                seed = seed,
                                return_se = TRUE,
                                ml = TRUE,
                                Y_Z_X_library = config$Y_Z_X_library,
                                Y_X_library = config$Y_X_library, 
                                S_X_library = config$S_X_library,
                                S_Z_X_library = config$S_Z_X_library,
                                family = "binomial",
                                return_models = FALSE,
                                effect_dir = "negative")
  
  results_df <- grid[i,]
  results_df <- cbind(results_df, data.frame(# Naturally infected AIPW
    nat_inf_pt_est = results$nat_inf$aipw$pt_est["additive_effect"],
    nat_inf_se = results$nat_inf$aipw$pt_est["additive_se"],
    nat_inf_reject = results$nat_inf$aipw$reject$additive,
    # Naturally infected upper bound
    nat_inf_upper_bound_pt_est = results$nat_inf$bound$pt_est["additive_effect_upper"],
    nat_inf_upper_bound_se = results$nat_inf$bound$boot_se$se_additive_upper,
    nat_inf_upper_bound_lower_ci = results$nat_inf$bound$boot_se$lower_ci_additive_upper,
    nat_inf_upper_bound_upper_ci = results$nat_inf$bound$boot_se$upper_ci_additive_upper,
    # Doomed AIPW
    doomed_pt_est = results$doomed$aipw$pt_est["additive_effect"],
    doomed_se = results$doomed$aipw$pt_est["additive_se"],
    doomed_reject = results$doomed$aipw$reject$additive,
    # Doomed upper bound
    doomed_upper_bound_pt_est = results$doomed$bound$pt_est["additive_effect_upper"],
    doomed_upper_bound_se = results$doomed$bound$boot_se$se_additive_upper,
    doomed_upper_bound_lower_ci = results$doomed$bound$boot_se$lower_ci_additive_upper,
    doomed_upper_bound_upper_ci = results$doomed$bound$boot_se$upper_ci_additive_upper,
    # Pop AIPW
    pop_pt_est = results$pop$aipw$pt_est["additive_effect"],
    pop_se = results$pop$aipw$pt_est["additive_se"],
    pop_reject = results$pop$aipw$reject$additive))

  # save individual row too in case job fails, can splice together later
  saveRDS(results, paste0(project_dir, setting, "/",
                          "n_", grid$n_sample_size[i],
                          "_inflation_", grid$inflation[i],
                          "_doomedepsilon_", grid$doomed_epsilon[i],
                          "_natinfepsilon_", grid$nat_inf_epsilon[i],
                          "_effectprotect_", grid$effect_protect[i], ".Rds"))
  
  return(results_df)

}, grid = grid, data = data, seed = seed, future.seed = TRUE)

full_results <- do.call(rbind, results)

saveRDS(full_results, paste0(project_dir, setting, "/", setting, "_overall.Rds"))
                        