# -----------------------------------------------------------------------------
# Script to get truth for given simulation settings
# -----------------------------------------------------------------------------

here::i_am("get_truth.R")

get_truth_provide <- function(config, n = 1e7, seed = 12345){
  
  grid <- expand.grid(effect_protect = config$effect_protect,
                      doomed_inflation = as.numeric(config$doomed_inflation),
                      protected_inflation = as.numeric(config$protected_inflation),
                      doomed_epsilon = as.numeric(config$doomed_epsilon),
                      protected_epsilon = as.numeric(config$protected_epsilon))
  
  truth <- cbind(grid, data.frame(E_Y1__protected_or_doomed = rep(NA, nrow(grid)),
                                  E_Y0__protected_or_doomed = rep(NA, nrow(grid)),
                                  E_Y1__doomed = rep(NA, nrow(grid)),
                                  E_Y0__doomed = rep(NA, nrow(grid)),
                                  E_Y1__protected = rep(NA, nrow(grid)),
                                  E_Y0__protected = rep(NA, nrow(grid)),
                                  E_Y1__pop = rep(NA, nrow(grid)),
                                  E_Y0__pop = rep(NA, nrow(grid)),
                                  nat_inf_upper_bound = rep(NA, nrow(grid)),
                                  doomed_upper_bound = rep(NA, nrow(grid))))
  
  rhobar_v_truth <- vector("list", length = nrow(grid))
  mubar_vs_truth <- vector("list", length = nrow(grid))
  
  for(i in 1:nrow(grid)){
    big_data <- simulate_data_provide(seed = seed,
                              effect_protect = grid$effect_protect[i],
                              protected_inflation = grid$protected_inflation[i],
                              doomed_inflation = grid$doomed_inflation[i],
                              protected_epsilon = grid$protected_epsilon[i], 
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
    
    nat_inf_bound <- get_bound_nat_inf(big_data, Y_name = "any_abx_wk52", Z_name = "rotaarm", S_name = "rotaepi", family = "binomial")
    doomed_bound <- get_bound_doomed(big_data, Y_name = "any_abx_wk52", Z_name = "rotaarm", S_name = "rotaepi", family = "binomial")
    
    truth$nat_inf_upper_bound[i] <- nat_inf_bound['additive_effect_upper']
    truth$doomed_upper_bound[i] <- doomed_bound['additive_effect_upper']
    
  }
  
  ## Point estimates effects
  truth$effect_nat_inf <- truth$E_Y1__protected_or_doomed - truth$E_Y0__protected_or_doomed
  truth$effect_doomed <- truth$E_Y1__doomed - truth$E_Y0__doomed
  truth$effect_protected <- truth$E_Y1__protected - truth$E_Y0__protected
  truth$effect_pop <- truth$E_Y1__pop - truth$E_Y0__pop
  
  return(truth)
}

get_truth_generic <- function(config, n = 1e7, seed = 12345){
  
  grid <- expand.grid(effect_protect = config$effect_protect,
                      doomed_inflation = as.numeric(config$doomed_inflation),
                      protected_epsilon = as.numeric(config$protected_epsilon),
                      doomed_epsilon = as.numeric(config$doomed_epsilon),
                      immune_epsilon = as.numeric(config$immune_epsilon))
  
  truth <- cbind(grid, data.frame(E_Y1__protected_or_doomed = rep(NA, nrow(grid)),
                                  E_Y0__protected_or_doomed = rep(NA, nrow(grid)),
                                  E_Y1__doomed = rep(NA, nrow(grid)),
                                  E_Y0__doomed = rep(NA, nrow(grid)),
                                  E_Y1__protected = rep(NA, nrow(grid)),
                                  E_Y0__protected = rep(NA, nrow(grid)),
                                  E_Y1__pop = rep(NA, nrow(grid)),
                                  E_Y0__pop = rep(NA, nrow(grid)),
                                  nat_inf_upper_bound = rep(NA, nrow(grid)),
                                  doomed_upper_bound = rep(NA, nrow(grid)),
                                  nat_inf_upper_bound_mult = rep(NA, nrow(grid)),
                                  doomed_upper_bound_mult = rep(NA, nrow(grid)),
                                  nat_inf_lower_bound = rep(NA, nrow(grid)),
                                  doomed_lower_bound = rep(NA, nrow(grid)),
                                  nat_inf_lower_bound_mult = rep(NA, nrow(grid)),
                                  doomed_lower_bound_mult = rep(NA, nrow(grid))))
  
  rhobar_v_truth <- vector("list", length = nrow(grid))
  mubar_vs_truth <- vector("list", length = nrow(grid))
  
  for(i in 1:nrow(grid)){
    big_data <- simulate_data_generic(seed = seed,
                                      effect_protect = grid$effect_protect[i],
                                      doomed_inflation = grid$doomed_inflation[i],
                                      protected_epsilon = grid$protected_epsilon[i], 
                                      doomed_epsilon = grid$doomed_epsilon[i],
                                      immune_epsilon = grid$immune_epsilon[i],
                                      n = n)
    
    # Naturally infected estimand
    truth$E_Y1__protected_or_doomed[i] <- mean(big_data$Y1[
        big_data$stratum %in% c("Protected", "Doomed")
    ])
    
    truth$E_Y0__protected_or_doomed[i] <- mean(big_data$Y0[
        big_data$stratum %in% c("Protected", "Doomed")
    ])
    
    # Doomed estimand
    truth$E_Y1__doomed[i] <- mean(big_data$Y1[
        big_data$stratum %in% c("Doomed")
    ])
    
    truth$E_Y0__doomed[i] <- mean(big_data$Y0[
        big_data$stratum %in% c("Doomed")
    ])
    
    # Protected estimand
    truth$E_Y1__protected[i] <- mean(big_data$Y1[
        big_data$stratum %in% c("Protected")
    ])
    
    truth$E_Y0__protected[i] <- mean(big_data$Y0[
        big_data$stratum %in% c("Protected")
    ])
    
    # Population estimand
    truth$E_Y1__pop[i] <- mean(big_data$Y1)
    
    truth$E_Y0__pop[i] <- mean(big_data$Y0)
    
    nat_inf_bound <- get_bound_nat_inf(big_data, Y_name = "Y", Z_name = "Z", S_name = "S", family = "binomial")
    doomed_bound <- get_bound_doomed(big_data, Y_name = "Y", Z_name = "Z", S_name = "S", family = "binomial")
    
    truth$nat_inf_upper_bound[i] <- nat_inf_bound['additive_effect_upper']
    truth$doomed_upper_bound[i] <- doomed_bound['additive_effect_upper']
    
    truth$nat_inf_upper_bound_mult[i] <- nat_inf_bound['mult_effect_upper']
    truth$doomed_upper_bound_mult[i] <- doomed_bound['mult_effect_upper']
    
    truth$nat_inf_lower_bound[i] <- nat_inf_bound['additive_effect_lower']
    truth$doomed_lower_bound[i] <- doomed_bound['additive_effect_lower']
    
    truth$nat_inf_lower_bound_mult[i] <- nat_inf_bound['mult_effect_lower']
    truth$doomed_lower_bound_mult[i] <- doomed_bound['mult_effect_lower']
    
  }
  
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
  
  return(truth)
}

get_truth_contour <- function(config, n = 1e7, seed = 12345){
  
  grid <- expand.grid(effect_protect = config$effect_protect,
                      doomed_inflation = as.numeric(config$doomed_inflation),
                      protected_inflation = as.numeric(config$protected_inflation),
                      doomed_epsilon = as.numeric(config$doomed_epsilon),
                      protected_epsilon = as.numeric(config$protected_epsilon))
  
  truth <- cbind(grid, data.frame(E_Y1__protected_or_doomed = rep(NA, nrow(grid)),
                                  E_Y0__protected_or_doomed = rep(NA, nrow(grid)),
                                  E_Y1__doomed = rep(NA, nrow(grid)),
                                  E_Y0__doomed = rep(NA, nrow(grid)),
                                  E_Y1__protected = rep(NA, nrow(grid)),
                                  E_Y0__protected = rep(NA, nrow(grid)),
                                  E_Y1__pop = rep(NA, nrow(grid)),
                                  E_Y0__pop = rep(NA, nrow(grid)),
                                  nat_inf_upper_bound = rep(NA, nrow(grid)),
                                  doomed_upper_bound = rep(NA, nrow(grid)),
                                  nat_inf_upper_bound_mult = rep(NA, nrow(grid)),
                                  doomed_upper_bound_mult = rep(NA, nrow(grid)),
                                  nat_inf_lower_bound = rep(NA, nrow(grid)),
                                  doomed_lower_bound = rep(NA, nrow(grid)),
                                  nat_inf_lower_bound_mult = rep(NA, nrow(grid)),
                                  doomed_lower_bound_mult = rep(NA, nrow(grid))))
  
  rhobar_v_truth <- vector("list", length = nrow(grid))
  mubar_vs_truth <- vector("list", length = nrow(grid))
  
  for(i in 1:nrow(grid)){
    big_data <- simulate_data_contour(seed = seed,
                                      effect_protect = grid$effect_protect[i],
                                      doomed_inflation = grid$doomed_inflation[i],
                                      protected_inflation = grid$protected_inflation[i],
                                      protected_epsilon = grid$protected_epsilon[i], 
                                      doomed_epsilon = grid$doomed_epsilon[i],
                                      n = n)
    
    # Naturally infected estimand
    truth$E_Y1__protected_or_doomed[i] <- mean(big_data$Y[
      big_data$Z == 1 &
        big_data$stratum %in% c("Protected", "Doomed")
    ])
    
    truth$E_Y0__protected_or_doomed[i] <- mean(big_data$Y[
      big_data$Z == 0 &
        big_data$stratum %in% c("Protected", "Doomed")
    ])
    
    # Doomed estimand
    truth$E_Y1__doomed[i] <- mean(big_data$Y[
      big_data$Z == 1 &
        big_data$stratum %in% c("Doomed")
    ])
    
    truth$E_Y0__doomed[i] <- mean(big_data$Y[
      big_data$Z == 0 &
        big_data$stratum %in% c("Doomed")
    ])
    
    # Protected estimand
    truth$E_Y1__protected[i] <- mean(big_data$Y[
      big_data$Z == 1 &
        big_data$stratum %in% c("Protected")
    ])
    
    truth$E_Y0__protected[i] <- mean(big_data$Y[
      big_data$Z == 0 &
        big_data$stratum %in% c("Protected")
    ])
    
    # Population estimand
    truth$E_Y1__pop[i] <- mean(big_data$Y[
      big_data$Z == 1 
    ])
    
    truth$E_Y0__pop[i] <- mean(big_data$Y[
      big_data$Z == 0 
    ])
    
    nat_inf_bound <- get_bound_nat_inf(big_data, Y_name = "Y", Z_name = "Z", S_name = "S", family = "binomial")
    doomed_bound <- get_bound_doomed(big_data, Y_name = "Y", Z_name = "Z", S_name = "S", family = "binomial")
    
    truth$nat_inf_upper_bound[i] <- nat_inf_bound['additive_effect_upper']
    truth$doomed_upper_bound[i] <- doomed_bound['additive_effect_upper']
    
    truth$nat_inf_upper_bound_mult[i] <- nat_inf_bound['mult_effect_upper']
    truth$doomed_upper_bound_mult[i] <- doomed_bound['mult_effect_upper']
    
    truth$nat_inf_lower_bound[i] <- nat_inf_bound['additive_effect_lower']
    truth$doomed_lower_bound[i] <- doomed_bound['additive_effect_lower']
    
    truth$nat_inf_lower_bound_mult[i] <- nat_inf_bound['mult_effect_lower']
    truth$doomed_lower_bound_mult[i] <- doomed_bound['mult_effect_lower']
    
  }
  
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
  
  return(truth)
}

get_truth_generic_bounds <- function(config, n = 1e7, seed = 12345){
  
  grid <- expand.grid(effect_protect = config$effect_protect,
                      doomed_inflation = as.numeric(config$doomed_inflation),
                      protected_epsilon = as.numeric(config$protected_epsilon),
                      doomed_epsilon = as.numeric(config$doomed_epsilon),
                      immune_epsilon = as.numeric(config$immune_epsilon))
  
  # covariate combinations
  X_combos <- c("X1", "X2", "X3",
                "X1X2", "X2X3", "X1X3",
                "X1X2X3")
  
  # all column names you want
  cov_bound_names <- c(
    paste0("cov_adj_upper_bound_", X_combos),
    paste0("cov_adj_lower_bound_", X_combos),
    paste0("cov_adj_upper_bound_mult_", X_combos),
    paste0("cov_adj_lower_bound_mult_", X_combos)
  )
  
  # baseline truth columns
  base_names <- c(
    "nat_inf_upper_bound",
    "nat_inf_upper_bound_mult",
    "nat_inf_lower_bound",
    "nat_inf_lower_bound_mult"
  )
  
  point_names <- c(
    "E_Y1__protected_or_doomed",
    "E_Y0__protected_or_doomed",
    "effect_nat_inf",
    "effect_nat_inf_mult"
  )
  
  all_truth_names <- c(point_names, base_names, cov_bound_names)
  
  na_cols <- data.frame(
    matrix(
      NA_real_,
      nrow = nrow(grid),
      ncol = length(all_truth_names),
      dimnames = list(NULL, all_truth_names)
    )
  )
  
  truth <- cbind(grid, na_cols)
  
  for(i in 1:nrow(grid)){
    big_data <- simulate_data_generic(seed = seed,
                                      effect_protect = grid$effect_protect[i],
                                      doomed_inflation = grid$doomed_inflation[i],
                                      protected_epsilon = grid$protected_epsilon[i], 
                                      doomed_epsilon = grid$doomed_epsilon[i],
                                      immune_epsilon = grid$immune_epsilon[i],
                                      n = n)
    
    big_data$X1X2 <- as.numeric(interaction(big_data$X1, big_data$X2))
    big_data$X2X3 <- as.numeric(interaction(big_data$X2, big_data$X3))
    big_data$X1X3 <- as.numeric(interaction(big_data$X1, big_data$X3))
    big_data$X1X2X3 <- as.numeric(interaction(big_data$X1, big_data$X2, big_data$X3))
    
    # Naturally infected estimand
    truth$E_Y1__protected_or_doomed[i] <- mean(big_data$Y[
      big_data$Z == 1 &
        big_data$stratum %in% c("Protected", "Doomed")
    ])
    
    truth$E_Y0__protected_or_doomed[i] <- mean(big_data$Y[
      big_data$Z == 0 &
        big_data$stratum %in% c("Protected", "Doomed")
    ])
    
    truth$effect_nat_inf[i] <-
      truth$E_Y1__protected_or_doomed[i] -
      truth$E_Y0__protected_or_doomed[i]
    
    truth$effect_nat_inf_mult[i] <-
      truth$E_Y1__protected_or_doomed[i] /
      truth$E_Y0__protected_or_doomed[i]
    
    ## ----------------------------
    ## Unadjusted bounds
    ## ----------------------------
    nat_inf_bound <- get_bound_nat_inf(
      big_data,
      Y_name = "Y",
      Z_name = "Z",
      S_name = "S",
      family = "binomial"
    )
    
    truth$nat_inf_upper_bound[i]        <- nat_inf_bound["additive_effect_upper"]
    truth$nat_inf_lower_bound[i]        <- nat_inf_bound["additive_effect_lower"]
    truth$nat_inf_upper_bound_mult[i]   <- nat_inf_bound["mult_effect_upper"]
    truth$nat_inf_lower_bound_mult[i]   <- nat_inf_bound["mult_effect_lower"]
    
    ## ----------------------------
    ## Covariate-adjusted bounds
    ## ----------------------------
    for(X in X_combos){
      
      cov_adj_bound <- get_cov_adj_bound_nat_inf(
        big_data,
        X_name = X,
        Y_name = "Y",
        Z_name = "Z",
        S_name = "S",
        family = "binomial"
      )
      
      truth[i, paste0("cov_adj_upper_bound_", X)]       <- cov_adj_bound["additive_effect_upper"]
      truth[i, paste0("cov_adj_lower_bound_", X)]       <- cov_adj_bound["additive_effect_lower"]
      truth[i, paste0("cov_adj_upper_bound_mult_", X)]  <- cov_adj_bound["mult_effect_upper"]
      truth[i, paste0("cov_adj_lower_bound_mult_", X)]  <- cov_adj_bound["mult_effect_lower"]
    }
  }
  
  return(truth)
}

# test for series of inflations to find which ones make effect size 0 for each estimand
# config <- config::get(file = here::here("config_contour.yml"), config = "contour_plot")
#config <- config::get(file = here::here("config_generic.yml"), config = "violate_nat_inf")
#truth <- get_truth_contour(config, n = 1e7)
#truth <- get_truth_provide(config, n = 1e6)
#truth <- get_truth_generic(config, n = 1e6)
