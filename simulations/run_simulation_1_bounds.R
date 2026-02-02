# ----------------------------------------------------------------------------
# Script to run 1st simulation in final manuscript - adding bounds
# 
# Generic data generation, bounds with/without covariate adjustment w different sets of assumptions,
# & violations to assumptions in DGP
# ----------------------------------------------------------------------------

# Path to installed packages on cluster
.libPaths(c("/apps/R/4.4.0/lib64/R/site/library","/apps/R/4.4.0/lib64/R/library", "~/Rlibs_ve_trial"))

here::i_am("run_simulation_1_bounds.R")

source(here::here("simulate_data.R"))

#devtools::load_all("~/Documents/shigella_projects/packages/vegrowth")

library(vegrowth)
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

project_dir <- "/projects/dbenkes/allison/vegrowth_analysis/results/sim_1_bounds/"
seed <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
setting <- Sys.getenv("SETTING")

set.seed(seed)

cfg <- yaml::read_yaml("config_sim_1_bounds.yml")
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
# results <- future.apply::future_lapply(1:nrow(grid), function(i, grid){
  
  data <- simulate_data_generic(seed = grid$seed[i],
                                effect_protect = grid$effect_protect[i],
                                doomed_inflation = grid$doomed_inflation[i],
                                protected_epsilon = grid$protected_epsilon[i], 
                                doomed_epsilon = grid$doomed_epsilon[i],
                                immune_epsilon = grid$immune_epsilon[i],
                                n = grid$n_sample_size[i])
  
  data$X1X2 <- as.numeric(interaction(data$X1, data$X2))
  data$X2X3 <- as.numeric(interaction(data$X2, data$X3))
  data$X1X3 <- as.numeric(interaction(data$X1, data$X3))
  data$X1X2X3 <- as.numeric(interaction(data$X1, data$X2, data$X3))
  
  results <- vegrowth(data = data, 
                      Y_name = "Y",
                      Z_name = "Z",
                      S_name = "S",
                      X_name = config$X,
                      estimand = config$estimand,
                      method = config$method,
                      exclusion_restriction = c(TRUE, FALSE), # do for both exclusion restriction scenarios
                      cross_world = c(TRUE, FALSE),           # do for cross-world both true and false
                      n_boot = 1000,
                      seed = seed,
                      return_se = TRUE,
                      ml = FALSE,
                      family = "binomial",
                      return_models = FALSE,
                      effect_dir = "positive",
                      epsilon = grid$protected_epsilon[i])
  
  return(results)
  
}, grid = grid)
#}, grid = grid, future.seed = seed)

# full_results <- do.call(rbind, results)

saveRDS(results, paste0(project_dir, setting, "/", setting, "_overall_seed_", seed, ".Rds"))
