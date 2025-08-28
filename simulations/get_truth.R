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
  
  rhobar_v_truth <- vector("list", length = nrow(grid))
  mubar_vs_truth <- vector("list", length = nrow(grid))
  
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
    
    # Bounds
    rhobar_v <- rep(NA,2)
    names(rhobar_v) <- c("rhobar_0", "rhobar_1")
    
    mubar_vs <- rep(NA, 4)
    names(mubar_vs) <- c("mubar_00", "mubar_01", "mubar_10", "mubar_11")
    
    for(v in c(0,1)){
      rhobar_v[[paste0("rhobar_",v)]] <- mean(big_data$rotaepi[big_data$rotaarm == v])
      
      for(s in c(0,1)){
        mubar_vs[[paste0("mubar_",v,s)]] <- mean(big_data$any_abx_wk52[big_data$rotaarm == v & big_data$rotaepi == s])
      }
    }
    
    rhobar_v_truth[[i]] <- rhobar_v
    mubar_vs_truth[[i]] <- mubar_vs
    
  }
  
  ## Point estimates effects
  truth$effect_nat_inf <- truth$E_Y1__protected_or_doomed - truth$E_Y0__protected_or_doomed
  truth$effect_doomed <- truth$E_Y1__doomed - truth$E_Y0__doomed
  truth$effect_protected <- truth$E_Y1__protected - truth$E_Y0__protected
  truth$effect_pop <- truth$E_Y1__pop - truth$E_Y0__pop
  
  ## Bounds effects
  
  rhobar_v <- do.call(rbind, rhobar_v_truth)
  mubar_vs <- do.call(rbind, mubar_vs_truth)
  
  # Naturally infected
  # upper bound on E[Y(1) | S(0) == 1]
  truth$upper_bound_E_Y1__protected_or_doomed <- 
    mubar_vs[,'mubar_11'] * (rhobar_v[,'rhobar_1'] / rhobar_v[,'rhobar_0']) +
    mubar_vs[,'mubar_10'] * (1 - (rhobar_v[,'rhobar_1'] / rhobar_v[,'rhobar_0']))
  
  truth$nat_inf_upper_bound <- truth$upper_bound_E_Y1__protected_or_doomed - truth$E_Y0__protected_or_doomed
  
  # Doomed
  # upper bound on the effect means we want the lower bound on E[Y(0) | S(0) == 1] 
  # in the package it looks like we're just using mubar_01??? 
  truth$lower_bound_E_Y0__doomed <- mubar_vs[,'mubar_01']
  truth$doomed_upper_bound <- truth$E_Y1__doomed - truth$lower_bound_E_Y0__doomed
  
  return(truth)
}

# test for series of inflations to find which ones make effect size 0 for each estimand
config <- config::get(file = here::here("config.yml"), config = "vary_inflation")
truth <- get_truth(config, n = 1e5)
