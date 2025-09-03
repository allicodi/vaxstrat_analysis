# ------------------------------------------------------
# Script to run simulation on cluster (same as truth but using future lapply)
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

# Path to projects folder where results will be saved
setting <- Sys.getenv("SETTING")
cfg <- yaml::read_yaml("config_contour.yml")
config <- cfg[[setting]]

grid <- expand.grid(seed = 12345,
                    effect_protect = config$effect_protect,
                    doomed_inflation = as.numeric(config$doomed_inflation),
                    protected_inflation = as.numeric(config$protected_inflation),
                    doomed_epsilon = as.numeric(config$doomed_epsilon),
                    protected_epsilon = as.numeric(config$protected_epsilon))

# eliminate combos where inflation in doomed > inflation in nat_inf
grid <- subset(grid, doomed_inflation <= protected_inflation)

results <- future.apply::future_lapply(1:nrow(grid), function(i, grid){
  big_data <- simulate_data_contour(seed = grid$seed[i],
                                    effect_protect = grid$effect_protect[i],
                                    doomed_inflation = grid$doomed_inflation[i],
                                    protected_inflation = grid$protected_inflation[i],
                                    protected_epsilon = grid$protected_epsilon[i], 
                                    doomed_epsilon = grid$doomed_epsilon[i],
                                    n = 5e6)
  
  truth <- cbind(grid[i,], data.frame(E_Y1__protected_or_doomed = rep(NA, 1),
                                  E_Y0__protected_or_doomed = rep(NA, 1),
                                  E_Y1__doomed = rep(NA, 1),
                                  E_Y0__doomed = rep(NA, 1),
                                  E_Y1__protected = rep(NA, 1),
                                  E_Y0__protected = rep(NA, 1),
                                  E_Y1__pop = rep(NA, 1),
                                  E_Y0__pop = rep(NA, 1)))
  
  # Naturally infected estimand
  truth$E_Y1__protected_or_doomed[1] <- mean(big_data$Y[
    big_data$Z == 1 &
      big_data$stratum %in% c("Protected", "Doomed")
  ])
  
  truth$E_Y0__protected_or_doomed[1] <- mean(big_data$Y[
    big_data$Z == 0 &
      big_data$stratum %in% c("Protected", "Doomed")
  ])
  
  # Doomed estimand
  truth$E_Y1__doomed[1] <- mean(big_data$Y[
    big_data$Z == 1 &
      big_data$stratum %in% c("Doomed")
  ])
  
  truth$E_Y0__doomed[1] <- mean(big_data$Y[
    big_data$Z == 0 &
      big_data$stratum %in% c("Doomed")
  ])
  
  # Protected estimand
  truth$E_Y1__protected[1] <- mean(big_data$Y[
    big_data$Z == 1 &
      big_data$stratum %in% c("Protected")
  ])
  
  truth$E_Y0__protected[1] <- mean(big_data$Y[
    big_data$Z == 0 &
      big_data$stratum %in% c("Protected")
  ])
  
  # Population estimand
  truth$E_Y1__pop[1] <- mean(big_data$Y[
    big_data$Z == 1 
  ])
  
  truth$E_Y0__pop[1] <- mean(big_data$Y[
    big_data$Z == 0 
  ])
  
  return(truth)
  
}, grid = grid, future.seed = TRUE)


truth <- do.call(rbind, results)

## Point estimates effects
truth$effect_nat_inf <- truth$E_Y1__protected_or_doomed - truth$E_Y0__protected_or_doomed
truth$effect_doomed <- truth$E_Y1__doomed - truth$E_Y0__doomed
truth$effect_protected <- truth$E_Y1__protected - truth$E_Y0__protected
truth$effect_pop <- truth$E_Y1__pop - truth$E_Y0__pop


## Point estimates effects multiplicative
truth$effect_nat_inf_mult <- truth$E_Y1__protected_or_doomed / truth$E_Y0__protected_or_doomed
truth$effect_doomed_mult <- truth$E_Y1__doomed / truth$E_Y0__doomed
truth$effect_protected_mult <- truth$E_Y1__protected / truth$E_Y0__protected
truth$effect_pop_mult <- truth$E_Y1__pop / truth$E_Y0__pop

# full_results <- do.call(rbind, results)

saveRDS(truth, paste0("results/contour/", setting, "_truth.Rds"))
