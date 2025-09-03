# ------------------------------------------------------
# Script to run simulation on cluster
# ------------------------------------------------------

# Path to installed packages on cluster
.libPaths(c("/apps/R/4.4.0/lib64/R/site/library","/apps/R/4.4.0/lib64/R/library", "~/Rlibs_ve_trial"))

here::i_am("run_analysis_generic.R")

source(here::here("simulate_data.R"))

#devtools::load_all("../../shigella_projects/packages/vegrowth/")

#library(future)
library(future.apply)
library(SuperLearner)
library(vegrowth)

# For initial debugging scratch file
options(echo = TRUE)

# this was for protect, probably can make smaller
options(future.globals.maxSize = 5 * 1024^3) # 5 GB
options(future.globals.onReference = "ignore")

ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
print(ncores)
plan(multisession, workers = ncores)

# ncores <- parallelly::availableCores()
# future::plan("multicore", workers = ncores)

# Path to projects folder where results will be saved
project_dir <- "/projects/dbenkes/allison/vegrowth_analysis/results/generic/"
seed <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
setting <- Sys.getenv("SETTING")

# read in config file
#config <- config::get(file = here::here("config_generic.yml"), config = setting)

cfg <- yaml::read_yaml("config_generic.yml")
config <- cfg[[setting]]

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

results <- future.apply::future_lapply(1:nrow(grid), function(i, grid){
  
  library(SuperLearner)
  
  data <- simulate_data_generic(seed = grid$seed[i],
                                effect_protect = grid$effect_protect[i],
                                inflation = grid$inflation[i],
                                nat_inf_epsilon = grid$nat_inf_epsilon[i], 
                                doomed_epsilon = grid$doomed_epsilon[i],
                                n = grid$n_sample_size[i])
  
  results <- vegrowth::vegrowth(data = data, 
                                Y_name = "Y",
                                Z_name = "Z",
                                S_name = "S",
                                X_name = c("X1", "X2", "X3"),
                                estimand = c("nat_inf", "doomed", "pop"),
                                method = c("gcomp", "ipw", "aipw", "tmle", "bound"),
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
                                effect_dir = "positive")
  
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
  
}, grid = grid, future.seed = TRUE)

# full_results <- do.call(rbind, results)

saveRDS(results, paste0(project_dir, setting, "/", setting, "_overall_seed_", seed, ".Rds"))
