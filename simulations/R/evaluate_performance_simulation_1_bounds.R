# ------------------------------------------------------------------------------
# Script to evaluate performance of generic simulations 
# ------------------------------------------------------------------------------

# Path to installed packages on cluster
options(echo = TRUE)

.libPaths(c("/apps/R/4.4.0/lib64/R/site/library","/apps/R/4.4.0/lib64/R/library", "~/Rlibs_ve_trial"))

here::i_am("R/evaluate_performance_simulation_1_bounds.R")

source(here::here("simulate_data.R"))
source(here::here("get_truth.R"))

library(vaxstrat)
library(dplyr)

# Path to projects folder where results will be saved
project_dir <- "/projects/dbenkes/allison/vegrowth_analysis/results/sim_1_bounds/"

# command args config setting
cargs <- commandArgs(trailingOnly = TRUE)
setting <- cargs[[1]]

cfg <- yaml::read_yaml("config_sim_1_bounds.yml")
config <- cfg[[setting]]

result_dir <- file.path(project_dir, setting)

# Get true effects --- will be same for all Xs, get the prefix

prefix_setting <- sub("_X.*$", "", setting)
if(prefix_setting == "default") prefix_setting <- "cw_er"

truth_file <- here::here(paste0("truth/bound_truth/", prefix_setting, "_truth.Rds"))
if(!file.exists(truth_file)){
  truth <- get_truth_generic_bounds(config)
  saveRDS(truth, truth_file)
} else{
  truth <- readRDS(truth_file)
}

# List all matching .Rds files
files <- list.files(
  path = result_dir,
  pattern = paste0("^", setting, "_overall_seed_.*\\.Rds$"),
  full.names = TRUE
)

all_result_df <- data.frame()

for(fname in files){
  res_obj <- readRDS(fname)
  seed <- as.integer(gsub(paste0(setting, "_overall_seed_|\\.Rds"), "", basename(fname)))
  
  grid <- expand.grid(seed = seed, 
                      n_sample_size = as.numeric(config$n_sample_size),
                      effect_protect = config$effect_protect,
                      doomed_inflation = as.numeric(config$doomed_inflation),
                      protected_epsilon = as.numeric(config$protected_epsilon),
                      doomed_epsilon = as.numeric(config$doomed_epsilon),
                      immune_epsilon = as.numeric(config$immune_epsilon),
                      X = config$X)
  
  for(i in 1:nrow(grid)){
    
    # grab row from truth that matches grid setting
    truth_grid <- truth[which(truth$effect_protect == grid$effect_protect[i] &
                                truth$doomed_inflation == grid$doomed_inflation[i] &
                                truth$protected_epsilon == grid$protected_epsilon[i] &
                                truth$doomed_epsilon == grid$doomed_epsilon[i] &
                                truth$immune_epsilon == grid$immune_epsilon[i]), ]
    
    res <- res_obj[[i]]
    
    result_df <- data.frame(
      seed = seed,
      n = grid$n_sample_size[i],
      effect_protect = grid$effect_protect[i],
      doomed_inflation = grid$doomed_inflation[i],
      protected_epsilon = grid$protected_epsilon[i],
      dooomed_epsilon = grid$doomed_epsilon[i],
      immune_epsilon = grid$immune_epsilon[i],
      X = grid$X[i],
      estimand = "nat_inf",
      method = "cov_adj_bound",
      
      # Point est for lower and upper bounds
      additive_effect_lower = res$nat_inf$cov_adj_bound$pt_est['additive_effect_lower'],
      additive_effect_upper = res$nat_inf$cov_adj_bound$pt_est['additive_effect_upper'],
      mult_effect_lower = res$nat_inf$cov_adj_bound$pt_est['mult_effect_lower'],
      mult_effect_upper = res$nat_inf$cov_adj_bound$pt_est['mult_effect_upper'],
      
      # CIs on lower and upper bounds
      additive_lower_lower_ci = res$nat_inf$cov_adj_bound$boot_se$lower_ci_additive_lower,
      additive_lower_upper_ci = res$nat_inf$cov_adj_bound$boot_se$upper_ci_additive_lower,
      additive_upper_lower_ci = res$nat_inf$cov_adj_bound$boot_se$lower_ci_additive_upper,
      additive_upper_upper_ci = res$nat_inf$cov_adj_bound$boot_se$upper_ci_additive_upper,
      
      mult_lower_lower_ci = res$nat_inf$cov_adj_bound$boot_se$lower_ci_mult_lower,
      mult_lower_upper_ci = res$nat_inf$cov_adj_bound$boot_se$upper_ci_mult_lower,
      mult_upper_lower_ci = res$nat_inf$cov_adj_bound$boot_se$lower_ci_mult_upper,
      mult_upper_upper_ci = res$nat_inf$cov_adj_bound$boot_se$upper_ci_mult_upper,
      
      # Truth Bounds
      additive_effect_lower_truth = truth_grid[[paste0('cov_adj_lower_bound_', grid$X[i])]],
      additive_effect_upper_truth = truth_grid[[paste0('cov_adj_upper_bound_', grid$X[i])]],
      mult_effect_lower_truth = truth_grid[[paste0('cov_adj_lower_bound_mult_', grid$X[i])]],
      mult_effect_upper_truth = truth_grid[[paste0('cov_adj_upper_bound_mult_', grid$X[i])]],
     
      # Truth pt est
      additive_effect_truth = truth_grid$effect_nat_inf,
      mult_effect_truth = truth_grid$effect_nat_inf_mult
    )
    
    ## REPEAT FOR NON-COVARIATE ADJUSTED VERSION IF APPLICABLE (will only be in the X1 configs, removed from others to avoid recalc unnecessarily)
    if(!is.null(res$nat_inf$bound)){
      result_df_2 <- data.frame(
        seed = seed,
        n = grid$n_sample_size[i],
        effect_protect = grid$effect_protect[i],
        doomed_inflation = grid$doomed_inflation[i],
        protected_epsilon = grid$protected_epsilon[i],
        dooomed_epsilon = grid$doomed_epsilon[i],
        immune_epsilon = grid$immune_epsilon[i],
        X = "NA",
        estimand = "nat_inf",
        method = "bound",
        
        # Point est for lower and upper bounds
        additive_effect_lower = res$nat_inf$bound$pt_est['additive_effect_lower'],
        additive_effect_upper = res$nat_inf$bound$pt_est['additive_effect_upper'],
        mult_effect_lower = res$nat_inf$bound$pt_est['mult_effect_lower'],
        mult_effect_upper = res$nat_inf$bound$pt_est['mult_effect_upper'],
        
        # CIs on lower and upper bounds
        additive_lower_lower_ci = res$nat_inf$bound$boot_se$lower_ci_additive_lower,
        additive_lower_upper_ci = res$nat_inf$bound$boot_se$upper_ci_additive_lower,
        additive_upper_lower_ci = res$nat_inf$bound$boot_se$lower_ci_additive_upper,
        additive_upper_upper_ci = res$nat_inf$bound$boot_se$upper_ci_additive_upper,
        
        mult_lower_lower_ci = res$nat_inf$bound$boot_se$lower_ci_mult_lower,
        mult_lower_upper_ci = res$nat_inf$bound$boot_se$upper_ci_mult_lower,
        mult_upper_lower_ci = res$nat_inf$bound$boot_se$lower_ci_mult_upper,
        mult_upper_upper_ci = res$nat_inf$bound$boot_se$upper_ci_mult_upper,
        
        # Truth Bounds
        additive_effect_lower_truth = truth_grid$nat_inf_lower_bound,
        additive_effect_upper_truth = truth_grid$nat_inf_upper_bound,
        mult_effect_lower_truth = truth_grid$nat_inf_lower_bound_mult,
        mult_effect_upper_truth = truth_grid$nat_inf_upper_bound_mult,
        
        # Truth pt est
        additive_effect_truth = truth_grid$effect_nat_inf,
        mult_effect_truth = truth_grid$effect_nat_inf_mult
      )
      
      result_df <- rbind(result_df, result_df_2)
    }

    # Difference between bound est & true bound
    result_df$lower_additive_diff <- result_df$additive_effect_lower - result_df$additive_effect_lower_truth
    result_df$upper_additive_diff <- result_df$additive_effect_upper - result_df$additive_effect_upper_truth
    result_df$lower_mult_diff <- log(result_df$mult_effect_lower) - log(result_df$mult_effect_lower_truth)
    result_df$upper_mult_diff <- log(result_df$mult_effect_upper) - log(result_df$mult_effect_upper_truth)
    
    # Width of bounds
    result_df$bound_width_additive <- result_df$additive_effect_upper - result_df$additive_effect_lower
    result_df$bound_width_mult <- result_df$mult_effect_upper - result_df$mult_effect_lower
    
    # Coverage of bounds
    result_df$lower_additive_coverage <- ifelse(result_df$additive_effect_lower_truth >= result_df$additive_lower_lower_ci &
                                            result_df$additive_effect_lower_truth <= result_df$additive_lower_upper_ci, 1, 0)
    
    result_df$upper_additive_coverage <- ifelse(result_df$additive_effect_upper_truth >= result_df$additive_upper_lower_ci &
                                                  result_df$additive_effect_upper_truth <= result_df$additive_upper_upper_ci, 1, 0)
    
    result_df$lower_mult_coverage <- ifelse(result_df$mult_effect_lower_truth >= result_df$mult_lower_lower_ci &
                                           result_df$mult_effect_lower_truth <= result_df$mult_lower_upper_ci, 1, 0)
    
    result_df$upper_mult_coverage <- ifelse(result_df$mult_effect_upper_truth >= result_df$mult_upper_lower_ci &
                                                  result_df$mult_effect_upper_truth <= result_df$mult_upper_upper_ci, 1, 0)
    
    # Coverage of true effect
    result_df$additive_effect_coverage <- ifelse(result_df$additive_effect_lower <= result_df$additive_effect_truth &
                                                   result_df$additive_effect_upper >= result_df$additive_effect_truth, 1, 0)
    
    result_df$mult_effect_coverage <- ifelse(result_df$mult_effect_lower <= result_df$mult_effect_truth &
                                                   result_df$mult_effect_upper >= result_df$mult_effect_truth, 1, 0)
    
    all_result_df <- rbind(all_result_df, result_df)
  }
  
}

saveRDS(all_result_df, here::here(paste0("results/sim_1_bounds/", setting, "_all.Rds")))

# Summarize bias, variance, MSE, coverage
summary_df <- all_result_df %>%
  group_by(estimand, method, n) %>%
  summarise(
    n = n[1],
    confirm_complete = n(),
    
    # Lower CI
    bias_lower_additive = mean(lower_additive_diff, na.rm = TRUE),
    coverage_lower_bound_additive = mean(lower_additive_coverage, na.rm = TRUE),
    bias_lower_mult = mean(lower_mult_diff, na.rm = TRUE),
    coverage_lower_bound_mult = mean(lower_mult_coverage, na.rm = TRUE),
    
    # Upper CI
    bias_upper_additive = mean(upper_additive_diff, na.rm = TRUE),
    coverage_upper_bound_additive = mean(upper_additive_coverage, na.rm = TRUE),
    bias_upper_mult = mean(upper_mult_diff, na.rm = TRUE),
    coverage_upper_bound_mult = mean(upper_mult_coverage, na.rm = TRUE),
    
    # Median (IQR) Bound Width
    median_width_additive = median(bound_width_additive, na.rm = TRUE),
    median_width_mult = median(bound_width_mult, na.rm = TRUE),
    
    q25_width_additive = quantile(bound_width_additive, 0.25, na.rm = TRUE),
    q25_width_mult = quantile(bound_width_mult, 0.25, na.rm = TRUE),
    q75_width_additive = quantile(bound_width_additive, 0.75, na.rm = TRUE),
    q75_width_mult = quantile(bound_width_mult, 0.75, na.rm = TRUE),
    
    # Point Estimate
    coverage_pt_est_additive = mean(additive_effect_coverage, na.rm = TRUE),
    coverage_pt_est_mult = mean(mult_effect_coverage, na.rm = TRUE),
    .groups = "drop"
  )

saveRDS(summary_df, paste0("results/sim_1_bounds/", setting, "_summary.Rds"))
