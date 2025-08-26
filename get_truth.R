# -----------------------------------------------------------------------------
# Script to get truth for given simulation settings
# -----------------------------------------------------------------------------

here::i_am("get_truth.R")

get_truth <- function(config, n = 1e7, seed = 12345){
  
  grid <- expand.grid(effect_protect = config$effect_protect,
                      inflation = as.numeric(config$inflation),
                      nat_inf_epsilon = as.numeric(config$nat_inf_epsilon),
                      doomed_epsilon = as.numeric(config$doomed_epsilon))
  
  truth <- cbind(grid, data.frame(E_Y1__protected_or_doomed = rep(NA, nrow(grid)),
                                  E_Y0__protected_or_doomed = rep(NA, nrow(grid)),
                                  E_Y1__doomed = rep(NA, nrow(grid)),
                                  E_Y0__doomed = rep(NA, nrow(grid)),
                                  E_Y1__protected = rep(NA, nrow(grid)),
                                  E_Y0__protected = rep(NA, nrow(grid)),
                                  E_Y1__pop = rep(NA, nrow(grid)),
                                  E_Y0__pop = rep(NA, nrow(grid))))
  
  for(i in 1:nrow(grid)){
    big_data <- simulate_data(seed = seed,
                              effect_protect = grid$effect_protect[i],
                              inflation = grid$inflation[i],
                              nat_inf_epsilon = grid$nat_inf_epsilon[i], 
                              doomed_epsilon = grid$doomed_epsilon[i],
                              n = n)
    
    # Naturally infected estimand
    truth$E_Y1__protected_or_doomed[i] <- mean(big_data$any_abx_wk52[
      big_data$rotaarm == 1 &
        big_data$stratum %in% c("Protected", "Doomed")
    ])
    
    truth$E_Y0__protected_or_doomed[i] <- mean(big_data$any_abx_wk52[
      big_data$rotaarm == 0 &
        big_data$stratum %in% c("Protected", "Doomed")
    ])
    
    # Doomed estimand
    truth$E_Y1__doomed[i] <- mean(big_data$any_abx_wk52[
      big_data$rotaarm == 1 &
        big_data$stratum %in% c("Doomed")
    ])
    
    truth$E_Y0__doomed[i] <- mean(big_data$any_abx_wk52[
      big_data$rotaarm == 0 &
        big_data$stratum %in% c("Doomed")
    ])
    
    # Protected estimand
    truth$E_Y1__protected[i] <- mean(big_data$any_abx_wk52[
      big_data$rotaarm == 1 &
        big_data$stratum %in% c("Protected")
    ])
    
    truth$E_Y0__protected[i] <- mean(big_data$any_abx_wk52[
      big_data$rotaarm == 0 &
        big_data$stratum %in% c("Protected")
    ])
    
    # Population estimand
    truth$E_Y1__pop[i] <- mean(big_data$any_abx_wk52[
      big_data$rotaarm == 1 
    ])
    
    truth$E_Y0__pop[i] <- mean(big_data$any_abx_wk52[
      big_data$rotaarm == 0 
    ])
    
  }
  
  truth$effect_nat_inf <- truth$E_Y1__protected_or_doomed - truth$E_Y0__protected_or_doomed
  truth$effect_doomed <- truth$E_Y1__doomed - truth$E_Y0__doomed
  truth$effect_protected <- truth$E_Y1__protected - truth$E_Y0__protected
  truth$effect_pop <- truth$E_Y1__pop - truth$E_Y0__pop
  
  return(truth)
}


# test for series of inflations to find which ones make effect size 0 for each estimand
config <- config::get(file = here::here("config.yml"), config = "vary_inflation")
truth <- get_truth(config, n = 1e7)
