# ----------------------------------------------------------------------------
# Script to run 1st simulation in final manuscript
# 
# Generic data generation, AIPW estimators w different sets of assumptions,
# & violations to assumptions in DGP
# ----------------------------------------------------------------------------

# Path to installed packages on cluster
.libPaths(c("/apps/R/4.4.0/lib64/R/site/library","/apps/R/4.4.0/lib64/R/library", "~/Rlibs_ve_trial"))

here::i_am("run_simulation_1.R")

source(here::here("simulate_data.R"))

devtools::load_all("~/Documents/shigella_projects/packages/vegrowth")

# library(vegrowth)
library(future.apply)
library(SuperLearner)

# For initial debugging scratch file
options(echo = TRUE)

# this was for protect, probably can make smaller
options(future.globals.maxSize = 5 * 1024^3) # 5 GB
options(future.globals.onReference = "ignore")

ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
print(ncores)
plan(multisession, workers = ncores)

project_dir <- "/projects/dbenkes/allison/vegrowth_analysis/results/sim_1/"
seed <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
setting <- Sys.getenv("SETTING")

cfg <- yaml::read_yaml("config_sim_1.yml")
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

results <- future.apply::future_lapply(1:nrow(grid), function(i, grid){
  
  library(SuperLearner)
  
  data <- simulate_data_generic(seed = grid$seed[i],
                                effect_protect = grid$effect_protect[i],
                                doomed_inflation = grid$doomed_inflation[i],
                                protected_epsilon = grid$protected_epsilon[i], 
                                doomed_epsilon = grid$doomed_epsilon[i],
                                immune_epsilon = grid$immune_epsilon[i],
                                n = grid$n_sample_size[i])
  
  results <- vegrowth(data = data, 
                      Y_name = "Y",
                      Z_name = "Z",
                      S_name = "S",
                      X_name = c("X1", "X2", "X3"),
                      estimand = config$estimand,
                      method = config$method,
                      exclusion_restriction = c(TRUE, FALSE), # do for both exclusion restriction scenarios
                      cross_world = c(TRUE, FALSE),           # do for cross-world both true and false
                      n_boot = 1000,
                      seed = seed,
                      return_se = TRUE,
                      ml = FALSE,
                      Y_Z_X_model = config$Y_Z_X,
                      Y_X_S1_model = config$Y_X_S1, 
                      Y_X_S0_model = config$Y_X_S0,
                      S_X_model = config$S_X,
                      S_Z_X_model = config$S_Z_X,
                      Z_X_model = config$Z_X, 
                      family = "binomial",
                      return_models = FALSE,
                      effect_dir = "positive",
                      epsilon = grid$protected_epsilon[i])
  
  # save results object, can splice together pieces to get bias, coverage, etc. later
  # except this is unnecessary? object not that big
  # saveRDS(results, paste0(project_dir, setting, "/",
  #                         "seed_", grid$seed[i],
  #                         "_n_", grid$n_sample_size[i],
  #                         "_inflation_", grid$inflation[i],
  #                         "_doomedepsilon_", grid$doomed_epsilon[i],
  #                         "_natinfepsilon_", grid$nat_inf_epsilon[i],
  #                         "_effectprotect_", grid$effect_protect[i], ".Rds"))
  
  return(results)
  
}, grid = grid, future.seed = seed)

# full_results <- do.call(rbind, results)

saveRDS(results, paste0(project_dir, setting, "/", setting, "_overall_seed_", seed, ".Rds"))
