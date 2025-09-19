# ------------------------------------------------------------------------------
# Script to evaluate performance of generic simulations 
# ------------------------------------------------------------------------------

# Path to installed packages on cluster
options(echo = TRUE)

.libPaths(c("/apps/R/4.4.0/lib64/R/site/library","/apps/R/4.4.0/lib64/R/library", "~/Rlibs_ve_trial"))

here::i_am("evaluate_performance_generic.R")

source(here::here("simulate_data.R"))
source(here::here("get_truth.R"))

library(vegrowth)
library(dplyr)

# Path to projects folder where results will be saved
project_dir <- "/projects/dbenkes/allison/vegrowth_analysis/results/generic_new/"

# command args config setting
cargs <- commandArgs(trailingOnly = TRUE)
setting <- cargs[[1]]

cfg <- yaml::read_yaml("config_generic.yml")
config <- cfg[[setting]]

result_dir <- file.path(project_dir, setting)

# Get true effects
truth_file <- here::here(paste0("truth/", setting, "_truth.Rds"))
if(!file.exists(truth_file)){
  truth <- get_truth_generic(config)
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
                      doomed_epsilon = as.numeric(config$doomed_epsilon))
  
  for(i in 1:nrow(grid)){
    
    res <- res_obj[[i]]
    
    result_df <- data.frame(
      seed = seed,
      n = grid$n_sample_size[i],
      effect_protect = grid$effect_protect[i],
      doomed_inflation = grid$doomed_inflation[i],
      protected_epsilon = grid$protected_epsilon[i],
      dooomed_epsilon = grid$doomed_epsilon[i],
      estimand = c(rep("nat_inf",4),
                   rep("doomed",3),
                   rep("pop",3)),
      method = c(c("gcomp", "ipw", "aipw", "tmle"),
                 c("gcomp", "ipw", "aipw"),
                 c("gcomp", "ipw", "aipw")),
      # Additive estimate and confidence interval
      additive_estimate = c(
        c(res$nat_inf$gcomp$pt_est['additive_effect'],
          res$nat_inf$ipw$pt_est['additive_effect'],
          res$nat_inf$aipw$pt_est['additive_effect'],
          res$nat_inf$tmle$pt_est['additive_effect']),
        c(res$doomed$gcomp$pt_est['additive_effect'],
          res$doomed$ipw$pt_est['additive_effect'],
          res$doomed$aipw$pt_est['additive_effect']),
        c(res$pop$gcomp$pt_est['additive_effect'],
          res$pop$ipw$pt_est['additive_effect'],
          res$pop$aipw$pt_est['additive_effect'])
      ),
      additive_lower_ci = c(
        c(res$nat_inf$gcomp$boot_se$lower_ci_additive,
          res$nat_inf$ipw$boot_se$lower_ci_additive,
          res$nat_inf$aipw$pt_est['additive_effect'] - 1.96*res$nat_inf$aipw$pt_est['additive_se'],
          res$nat_inf$tmle$pt_est['additive_effect'] - 1.96*res$nat_inf$tmle$pt_est['additive_se']),
        c(res$doomed$gcomp$boot_se$lower_ci_additive,
          res$doomed$ipw$boot_se$lower_ci_additive,
          res$doomed$aipw$pt_est['additive_effect'] - 1.96*res$doomed$aipw$pt_est['additive_se']),
        c(res$pop$gcomp$boot_se$lower_ci_additive,
          res$pop$ipw$boot_se$lower_ci_additive,
          res$pop$aipw$pt_est['additive_effect'] - 1.96*res$pop$aipw$pt_est['additive_se'])
      ),
      additive_upper_ci = c(
        c(res$nat_inf$gcomp$boot_se$upper_ci_additive,
          res$nat_inf$ipw$boot_se$upper_ci_additive,
          res$nat_inf$aipw$pt_est['additive_effect'] + 1.96*res$nat_inf$aipw$pt_est['additive_se'],
          res$nat_inf$tmle$pt_est['additive_effect'] + 1.96*res$nat_inf$tmle$pt_est['additive_se']),
        c(res$doomed$gcomp$boot_se$upper_ci_additive,
          res$doomed$ipw$boot_se$upper_ci_additive,
          res$doomed$aipw$pt_est['additive_effect'] + 1.96*res$doomed$aipw$pt_est['additive_se']),
        c(res$pop$gcomp$boot_se$upper_ci_additive,
          res$pop$ipw$boot_se$upper_ci_additive,
          res$pop$aipw$pt_est['additive_effect'] + 1.96*res$pop$aipw$pt_est['additive_se'])
      ),
      # Multiplicative estimate and confidence interval
      mult_estimate = c(
        c(exp(as.numeric(res$nat_inf$gcomp$pt_est['log_multiplicative_effect'])),
          exp(as.numeric(res$nat_inf$ipw$pt_est['log_multiplicative_effect'])),
          exp(as.numeric(res$nat_inf$aipw$pt_est['log_multiplicative_effect'])),
          exp(as.numeric(res$nat_inf$tmle$pt_est['log_multiplicative_effect']))),
        c(exp(as.numeric(res$doomed$gcomp$pt_est['log_multiplicative_effect'])),
          exp(as.numeric(res$doomed$ipw$pt_est['log_multiplicative_effect'])),
          exp(as.numeric(res$doomed$aipw$pt_est['log_multiplicative_effect']))),
        c(exp(as.numeric(res$pop$gcomp$pt_est['log_multiplicative_effect'])),
          exp(as.numeric(res$pop$ipw$pt_est['log_multiplicative_effect'])),
          exp(as.numeric(res$pop$aipw$pt_est['log_multiplicative_effect'])))
      ),
      mult_lower_ci = c(
        c(res$nat_inf$gcomp$boot_se$lower_ci_mult,
          res$nat_inf$ipw$boot_se$lower_ci_mult,
          exp(res$nat_inf$aipw$pt_est['log_multiplicative_effect'] - 1.96*res$nat_inf$aipw$pt_est['log_multiplicative_se']),
          exp(res$nat_inf$tmle$pt_est['log_multiplicative_effect'] - 1.96*res$nat_inf$tmle$pt_est['log_multiplicative_se'])),
        c(res$doomed$gcomp$boot_se$lower_ci_mult,
          res$doomed$ipw$boot_se$lower_ci_mult,
          exp(res$doomed$aipw$pt_est['log_multiplicative_effect'] - 1.96*res$doomed$aipw$pt_est['log_multiplicative_se'])),
        c(res$pop$gcomp$boot_se$lower_ci_mult,
          res$pop$ipw$boot_se$lower_ci_mult,
          exp(res$pop$aipw$pt_est['log_multiplicative_effect'] - 1.96*res$pop$aipw$pt_est['log_multiplicative_se']))
      ),
      mult_upper_ci = c(
        c(res$nat_inf$gcomp$boot_se$upper_ci_mult,
          res$nat_inf$ipw$boot_se$upper_ci_mult,
          exp(res$nat_inf$aipw$pt_est['log_multiplicative_effect'] + 1.96*res$nat_inf$aipw$pt_est['log_multiplicative_se']),
          exp(res$nat_inf$tmle$pt_est['log_multiplicative_effect'] + 1.96*res$nat_inf$tmle$pt_est['log_multiplicative_se'])),
        c(res$doomed$gcomp$boot_se$upper_ci_mult,
          res$doomed$ipw$boot_se$upper_ci_mult,
          exp(res$doomed$aipw$pt_est['log_multiplicative_effect'] + 1.96*res$doomed$aipw$pt_est['log_multiplicative_se'])),
        c(res$pop$gcomp$boot_se$upper_ci_mult,
          res$pop$ipw$boot_se$upper_ci_mult,
          exp(res$pop$aipw$pt_est['log_multiplicative_effect'] + 1.96*res$pop$aipw$pt_est['log_multiplicative_se']))
      ),
      # NEW psi_1 and psi_0
      psi_1 = c(
        c(res$nat_inf$gcomp$pt_est['psi_1'],
          res$nat_inf$ipw$pt_est['psi_1'],
          res$nat_inf$aipw$pt_est['psi_1'],
          res$nat_inf$tmle$pt_est['psi_1']),
        c(res$doomed$gcomp$pt_est['psi_1'],
          res$doomed$ipw$pt_est['psi_1'],
          res$doomed$aipw$pt_est['psi_1']),
        c(res$pop$gcomp$pt_est['psi_1'],
          res$pop$ipw$pt_est['psi_1'],
          res$pop$aipw$pt_est['psi_1'])
      ),
      psi_0 = c(
        c(res$nat_inf$gcomp$pt_est['psi_0'],
          res$nat_inf$ipw$pt_est['psi_0'],
          res$nat_inf$aipw$pt_est['psi_0'],
          res$nat_inf$tmle$pt_est['psi_0']),
        c(res$doomed$gcomp$pt_est['psi_0'],
          res$doomed$ipw$pt_est['psi_0'],
          res$doomed$aipw$pt_est['psi_0']),
        c(res$pop$gcomp$pt_est['psi_0'],
          res$pop$ipw$pt_est['psi_0'],
          res$pop$aipw$pt_est['psi_0'])
      ),
      # Truth
      additive_truth = c(rep(truth$effect_nat_inf[truth$doomed_inflation == grid$doomed_inflation[i] & truth$protected_epsilon == grid$protected_epsilon[i] & truth$doomed_epsilon == grid$doomed_epsilon[i]], 4), 
                         rep(truth$effect_doomed[truth$doomed_inflation == grid$doomed_inflation[i] & truth$protected_epsilon == grid$protected_epsilon[i] & truth$doomed_epsilon == grid$doomed_epsilon[i]], 3), 
                         rep(truth$effect_pop[truth$doomed_inflation == grid$doomed_inflation[i] & truth$protected_epsilon == grid$protected_epsilon[i] & truth$doomed_epsilon == grid$doomed_epsilon[i]], 3)),
      mult_truth = c(rep(truth$effect_nat_inf_mult[truth$doomed_inflation == grid$doomed_inflation[i] & truth$protected_epsilon == grid$protected_epsilon[i] & truth$doomed_epsilon == grid$doomed_epsilon[i]], 4),
                     rep(truth$effect_doomed_mult[truth$doomed_inflation == grid$doomed_inflation[i] & truth$protected_epsilon == grid$protected_epsilon[i] & truth$doomed_epsilon == grid$doomed_epsilon[i]], 3),
                     rep(truth$effect_pop_mult[truth$doomed_inflation == grid$doomed_inflation[i] & truth$protected_epsilon == grid$protected_epsilon[i] & truth$doomed_epsilon == grid$doomed_epsilon[i]], 3)),
      
      psi_1_truth = c(rep(truth$E_Y1__protected_or_doomed[truth$doomed_inflation == grid$doomed_inflation[i] & truth$protected_epsilon == grid$protected_epsilon[i] & truth$doomed_epsilon == grid$doomed_epsilon[i]], 4), 
                         rep(truth$E_Y1__doomed[truth$doomed_inflation == grid$doomed_inflation[i] & truth$protected_epsilon == grid$protected_epsilon[i] & truth$doomed_epsilon == grid$doomed_epsilon[i]], 3), 
                         rep(truth$E_Y1__pop[truth$doomed_inflation == grid$doomed_inflation[i] & truth$protected_epsilon == grid$protected_epsilon[i] & truth$doomed_epsilon == grid$doomed_epsilon[i]], 3)),
      psi_0_truth = c(rep(truth$E_Y0__protected_or_doomed[truth$doomed_inflation == grid$doomed_inflation[i] & truth$protected_epsilon == grid$protected_epsilon[i] & truth$doomed_epsilon == grid$doomed_epsilon[i]], 4), 
                         rep(truth$E_Y0__doomed[truth$doomed_inflation == grid$doomed_inflation[i] & truth$protected_epsilon == grid$protected_epsilon[i] & truth$doomed_epsilon == grid$doomed_epsilon[i]], 3), 
                         rep(truth$E_Y0__pop[truth$doomed_inflation == grid$doomed_inflation[i] & truth$protected_epsilon == grid$protected_epsilon[i] & truth$doomed_epsilon == grid$doomed_epsilon[i]], 3))
    )
    
    result_df$additive_diff <- result_df$additive_estimate - result_df$additive_truth
    result_df$mult_diff <- result_df$mult_estimate - result_df$mult_truth
    
    result_df$psi_1_diff <- result_df$psi_1 - result_df$psi_1_truth
    result_df$psi_0_diff <- result_df$psi_0 - result_df$psi_0_truth
    
    result_df$additive_coverage <- ifelse(result_df$additive_truth >= result_df$additive_lower_ci &
                                            result_df$additive_truth <= result_df$additive_upper_ci, 1, 0)
    result_df$mult_coverage <- ifelse(result_df$mult_truth >= result_df$mult_lower_ci &
                                        result_df$mult_truth <= result_df$mult_upper_ci, 1, 0)
    
    all_result_df <- rbind(all_result_df, result_df)
  }
  
}

# Summarize bias, variance, MSE, coverage
summary_df <- all_result_df %>%
  group_by(estimand, method, n) %>%
  summarise(
    n = n[1],
    confirm_complete = n(),
    bias_additive = mean(additive_diff, na.rm = TRUE),
    var_additive = var(additive_estimate, na.rm = TRUE),
    mse_additive = mean(additive_diff^2, na.rm = TRUE),
    coverage_additive = mean(additive_coverage, na.rm = TRUE),
    
    bias_mult = mean(mult_diff, na.rm = TRUE),
    var_mult = var(mult_estimate, na.rm = TRUE),
    mse_mult = mean(mult_diff^2, na.rm = TRUE),
    coverage_mult = mean(mult_coverage, na.rm = TRUE),
    
    bias_psi_1 = mean(psi_1_diff, na.rm = TRUE),
    bias_psi_0 = mean(psi_0_diff, na.rm = TRUE),
    
    .groups = "drop"
  )

saveRDS(summary_df, paste0("results/generic_new/", setting, "_summary.Rds"))
