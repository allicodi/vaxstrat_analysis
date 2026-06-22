# ----------------------------------------------------------------------------
# Script to run cross-fitting simulation
# 
# cross-fit more complex data generation
# ----------------------------------------------------------------------------

# Path to installed packages on cluster
.libPaths(c("/apps/R/4.4.0/lib64/R/site/library","/apps/R/4.4.0/lib64/R/library", "~/Rlibs_ve_trial"))

here::i_am("R/run_cross_fit_simulation.R")

source(here::here("R/simulate_data.R"))

#devtools::load_all("~/Documents/shigella_projects/packages/vaxstrat")

# dev version of vaxstrat on crossfit branch on the cluster
devtools::load_all("~/vaxstrat")

# library(vaxstrat)
# library(future.apply)
library(SuperLearner)

# For initial debugging scratch file
options(echo = TRUE)

project_dir <- "/projects/dbenkes/allison/vegrowth_analysis/results/cross_fit/"
seed <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
setting <- Sys.getenv("SETTING")

cfg <- yaml::read_yaml("config_cross_fit.yml")
config <- cfg[[setting]]

# Create dir to save results if does not exist
if (!file.exists(paste0(project_dir, setting))) {
  dir.create(here::here(paste0(project_dir, setting)), recursive = TRUE)
}

# get all parameter grid combinations
grid <- expand.grid(seed = seed, 
                    n_sample_size = as.numeric(config$n_sample_size),
                    effect_protect = config$effect_protect,
                    doomed_inflation = as.numeric(config$doomed_inflation),
                    protected_epsilon = as.numeric(config$protected_epsilon),
                    doomed_epsilon = as.numeric(config$doomed_epsilon),
                    immune_epsilon = as.numeric(config$immune_epsilon))

results <- lapply(1:nrow(grid), function(i, grid){

  data <- simulate_data_cross_fit(seed = grid$seed[i],
                                  effect_protect = grid$effect_protect[i],
                                  doomed_inflation = grid$doomed_inflation[i],
                                  protected_epsilon = grid$protected_epsilon[i], 
                                  doomed_epsilon = grid$doomed_epsilon[i],
                                  immune_epsilon = grid$immune_epsilon[i],
                                  n = grid$n_sample_size[i])
  
  results <- vaxstrat(data = data, 
                      Y_name = "Y",
                      Z_name = "Z",
                      S_name = "S",
                      X_name = c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10"),
                      estimand = config$estimand,
                      method = config$method,
                      exclusion_restriction = c(TRUE, FALSE), # do for both exclusion restriction scenarios
                      cross_world = c(TRUE, FALSE),           # do for cross-world both true and false
                      seed = seed,
                      return_se = TRUE,
                      ml = TRUE,
                      Y_Z_X_library = config$Y_Z_X_library,
                      Y_X_library =  config$Y_X_library, 
                      S_X_library = config$S_X_library,
                      S_Z_X_library = config$S_Z_X_library, 
                      Z_X_library = config$Z_X_library,
                      family = "binomial",
                      return_models = FALSE,
                      effect_dir = "positive",
                      epsilon = grid$protected_epsilon[i], 
                      cross_fit = TRUE, 
                      cf_folds = config$cf_folds, 
                      id_name = "id")
  return(results)
  
}, grid = grid)

saveRDS(results, paste0(project_dir, setting, "/", setting, "_cross_fit_seed_", seed, ".Rds"))

# Repeat for non-cross fit version
results_og <- lapply(1:nrow(grid), function(i, grid){
  
  data <- simulate_data_cross_fit(seed = grid$seed[i],
                                  effect_protect = grid$effect_protect[i],
                                  doomed_inflation = grid$doomed_inflation[i],
                                  protected_epsilon = grid$protected_epsilon[i], 
                                  doomed_epsilon = grid$doomed_epsilon[i],
                                  immune_epsilon = grid$immune_epsilon[i],
                                  n = grid$n_sample_size[i])
  
  results <- vaxstrat(data = data, 
                      Y_name = "Y",
                      Z_name = "Z",
                      S_name = "S",
                      X_name = c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10"),
                      estimand = config$estimand,
                      method = config$method,
                      exclusion_restriction = c(TRUE, FALSE), # do for both exclusion restriction scenarios
                      cross_world = c(TRUE, FALSE),           # do for cross-world both true and false
                      seed = seed,
                      return_se = TRUE,
                      ml = TRUE,
                      Y_Z_X_library = config$Y_Z_X_library,
                      Y_X_library =  config$Y_X_library, 
                      S_X_library = config$S_X_library,
                      S_Z_X_library = config$S_Z_X_library, 
                      Z_X_library = config$Z_X_library,
                      family = "binomial",
                      return_models = FALSE,
                      effect_dir = "positive",
                      epsilon = grid$protected_epsilon[i], 
                      cross_fit = FALSE)
  return(results)
  
}, grid = grid)

saveRDS(results_og, paste0(project_dir, setting, "/", setting, "_single_fit_seed_", seed, ".Rds"))

