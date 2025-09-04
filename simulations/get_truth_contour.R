# ------------------------------------------------------
# Script to get truth for contour plots (run locally)
# ------------------------------------------------------

here::i_am("get_truth_contour.R")

source(here::here("simulate_data.R"))

devtools::load_all("../../shigella_projects/packages/vegrowth/")

library(SuperLearner)

# set to "provide_contour_plot" or "generic_contour_plot"; set n large (but if too large won't run well)
setting <- "provide_contour_plot"
n <- 1e6
seed <- 12345

cfg <- yaml::read_yaml("config_contour.yml")
config <- cfg[[setting]]

# Generic or provide
sim_type <- config$sim_type

grid <- expand.grid(seed = seed,
                    effect_protect = config$effect_protect,
                    doomed_inflation = as.numeric(config$doomed_inflation),
                    protected_inflation = as.numeric(config$protected_inflation),
                    doomed_epsilon = as.numeric(config$doomed_epsilon),
                    protected_epsilon = as.numeric(config$protected_epsilon))

# eliminate combos where inflation in doomed > inflation in nat_inf (in generic version, this is always true for provide)
if(sim_type == "generic"){
  grid <- subset(grid, doomed_inflation <= protected_inflation)
}

results <- lapply(1:nrow(grid), function(i, grid, sim_type, n){
  
  truth <- cbind(grid[i,], data.frame(E_Y1__protected_or_doomed = rep(NA, 1),
                                      E_Y0__protected_or_doomed = rep(NA, 1),
                                      E_Y1__doomed = rep(NA, 1),
                                      E_Y0__doomed = rep(NA, 1),
                                      E_Y1__protected = rep(NA, 1),
                                      E_Y0__protected = rep(NA, 1),
                                      E_Y1__pop = rep(NA, 1),
                                      E_Y0__pop = rep(NA, 1)))
  
  if(sim_type == "generic"){
    big_data <- simulate_data_contour(seed = grid$seed[i],
                                      effect_protect = grid$effect_protect[i],
                                      doomed_inflation = grid$doomed_inflation[i],
                                      protected_inflation = grid$protected_inflation[i],
                                      protected_epsilon = grid$protected_epsilon[i], 
                                      doomed_epsilon = grid$doomed_epsilon[i],
                                      n = n)
    
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
    
  } else{
    big_data <- simulate_data_provide(seed = grid$seed[i],
                                      effect_protect = grid$effect_protect[i],
                                      doomed_inflation = grid$doomed_inflation[i],
                                      protected_inflation = grid$protected_inflation[i],
                                      protected_epsilon = grid$protected_epsilon[i], 
                                      doomed_epsilon = grid$doomed_epsilon[i],
                                      n = n)
    
    # Naturally infected estimand
    truth$E_Y1__protected_or_doomed[1] <- mean(big_data$any_abx_wk52[
      big_data$rotaarm == 1 &
        big_data$stratum %in% c("Protected", "Doomed")
    ])
    
    truth$E_Y0__protected_or_doomed[1] <- mean(big_data$any_abx_wk52[
      big_data$rotaarm == 0 &
        big_data$stratum %in% c("Protected", "Doomed")
    ])
    
    # Doomed estimand
    truth$E_Y1__doomed[1] <- mean(big_data$any_abx_wk52[
      big_data$rotaarm == 1 &
        big_data$stratum %in% c("Doomed")
    ])
    
    truth$E_Y0__doomed[1] <- mean(big_data$any_abx_wk52[
      big_data$rotaarm == 0 &
        big_data$stratum %in% c("Doomed")
    ])
    
    # Protected estimand
    truth$E_Y1__protected[1] <- mean(big_data$any_abx_wk52[
      big_data$rotaarm == 1 &
        big_data$stratum %in% c("Protected")
    ])
    
    truth$E_Y0__protected[1] <- mean(big_data$any_abx_wk52[
      big_data$rotaarm == 0 &
        big_data$stratum %in% c("Protected")
    ])
    
    # Population estimand
    truth$E_Y1__pop[1] <- mean(big_data$any_abx_wk52[
      big_data$rotaarm == 1 
    ])
    
    truth$E_Y0__pop[1] <- mean(big_data$any_abx_wk52[
      big_data$rotaarm == 0 
    ])

  }
  
  return(truth)
  
}, grid = grid, sim_type = sim_type, n = n)

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

saveRDS(truth, paste0("results/contour/", setting, "_truth.Rds"))