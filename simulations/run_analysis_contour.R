# ------------------------------------------------------
# Script to run simulation for power contour plots on cluster
# ------------------------------------------------------
options(echo = TRUE)

# Path to installed packages on cluster
.libPaths(c("/apps/R/4.4.0/lib64/R/site/library","/apps/R/4.4.0/lib64/R/library", "~/Rlibs_ve_trial"))

here::i_am("run_analysis_contour.R")

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

# Read config setting
setting <- Sys.getenv("SETTING")
cfg <- yaml::read_yaml("config_contour.yml")
config <- cfg[[setting]]

# Get seed
seed <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Generic or provide
sim_type <- config$sim_type

grid <- expand.grid(seed = seed,
                    effect_protect = config$effect_protect,
                    doomed_inflation = as.numeric(config$doomed_inflation),
                    protected_inflation = as.numeric(config$protected_inflation),
                    doomed_epsilon = as.numeric(config$doomed_epsilon),
                    protected_epsilon = as.numeric(config$protected_epsilon),
                    immune_delta = as.numeric(config$immune_delta),
                    protected_delta = as.numeric(config$protected_delta))

# eliminate combos where inflation in doomed > inflation in nat_inf (in generic version, this is always true for provide)
if(sim_type == "generic"){
  grid <- subset(grid, doomed_inflation <= protected_inflation)
}

results <- future.apply::future_lapply(1:nrow(grid), function(i, grid, sim_type, config){
  
  library(SuperLearner)
  
  if(sim_type == "generic"){
    big_data <- simulate_data_contour(seed = grid$seed[i],
                                      effect_protect = grid$effect_protect[i],
                                      doomed_inflation = grid$doomed_inflation[i],
                                      protected_inflation = grid$protected_inflation[i],
                                      protected_epsilon = grid$protected_epsilon[i], 
                                      doomed_epsilon = grid$doomed_epsilon[i],
                                      n = 2000)
    
    results <- vegrowth::vegrowth(data = big_data, 
                                  Y_name = "Y",
                                  Z_name = "Z",
                                  S_name = "S",
                                  X_name = c("X1", "X2", "X3"),
                                  estimand = c("nat_inf", "doomed", "pop"),
                                  method = c("aipw"),
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
    
  } else{
    big_data <- simulate_data_provide(seed = grid$seed[i],
                                      effect_protect = grid$effect_protect[i],
                                      doomed_inflation = grid$doomed_inflation[i],
                                      protected_inflation = grid$protected_inflation[i],
                                      protected_epsilon = grid$protected_epsilon[i], 
                                      doomed_epsilon = grid$doomed_epsilon[i],
                                      immune_delta = grid$immune_delta[i],
                                      protected_delta = grid$protected_delta[i],
                                      n = 700)
    
    results <- vegrowth::vegrowth(data = big_data, 
                                  Y_name = "any_abx_wk52",
                                  Z_name = "rotaarm",
                                  S_name = "rotaepi",
                                  X_name = c("wk10_haz", "gender", "num_hh_sleep"),
                                  estimand = c("nat_inf", "doomed", "pop"),
                                  method = c("aipw"),
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
  }
  
  results_df <- grid[i,]
  results_df <- cbind(results_df, data.frame(# Naturally infected AIPW
    nat_inf_pt_est = results$nat_inf$aipw$pt_est["additive_effect"],
    nat_inf_se = results$nat_inf$aipw$pt_est["additive_se"],
    nat_inf_reject = results$nat_inf$aipw$reject$additive,
    # Doomed AIPW
    doomed_pt_est = results$doomed$aipw$pt_est["additive_effect"],
    doomed_se = results$doomed$aipw$pt_est["additive_se"],
    doomed_reject = results$doomed$aipw$reject$additive,
    # Pop AIPW
    pop_pt_est = results$pop$aipw$pt_est["additive_effect"],
    pop_se = results$pop$aipw$pt_est["additive_se"],
    pop_reject = results$pop$aipw$reject$additive))
  
  return(results_df)

}, grid = grid, sim_type = sim_type, config = config, future.seed = TRUE)

results <- do.call(rbind, results)

saveRDS(results, paste0("results/contour/", setting, "_seed_", seed, ".Rds"))
