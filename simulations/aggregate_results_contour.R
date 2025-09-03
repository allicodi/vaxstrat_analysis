# Path to installed packages on cluster
.libPaths(c("/apps/R/4.4.0/lib64/R/site/library","/apps/R/4.4.0/lib64/R/library", "~/Rlibs_ve_trial"))

here::i_am("aggregate_results_contour.R")

source(here::here("simulate_data.R"))
source(here::here("get_truth.R"))

library(vegrowth)
library(dplyr)

# Path to projects folder where results will be saved
project_dir <- "/projects/dbenkes/allison/vegrowth_analysis/results/contour/"

# command args config setting
cargs <- commandArgs(trailingOnly = TRUE)
setting <- cargs[[1]]

cfg <- yaml::read_yaml("config_contour.yml")
config <- cfg[[setting]]

result_dir <- file.path(project_dir, setting)

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
                      inflation = as.numeric(config$inflation),
                      nat_inf_epsilon = as.numeric(config$nat_inf_epsilon),
                      doomed_epsilon = as.numeric(config$doomed_epsilon))
  
  for(i in 1:nrow(grid)){
    
    res <- res_obj[[i]]
    
    result_df <- data.frame(
      seed = seed,
      n = grid$n_sample_size[i],
      effect_protect = grid$effect_protect[i],
      inflation = grid$inflation[i],
      nat_inf_epsilon = grid$nat_inf_epsilon[i],
      dooomed_epsilon = grid$doomed_epsilon[i],
      estimand = c(rep("nat_inf",5),
                   rep("doomed",4),
                   rep("pop",3)),
      method = c(c("gcomp", "ipw", "aipw", "tmle", "lower bound"),
                 c("gcomp", "ipw", "aipw", "lower bound"),
                 c("gcomp", "ipw", "aipw")),
      # Additive estimate and confidence interval
      additive_estimate = c(
        c(res$nat_inf$gcomp$pt_est['additive_effect'],
          res$nat_inf$ipw$pt_est['additive_effect'],
          res$nat_inf$aipw$pt_est['additive_effect'],
          res$nat_inf$tmle$pt_est['additive_effect'],
          res$nat_inf$bound$pt_est['additive_effect_lower']),
        c(res$doomed$gcomp$pt_est['additive_effect'],
          res$doomed$ipw$pt_est['additive_effect'],
          res$doomed$aipw$pt_est['additive_effect'],
          res$doomed$bound$pt_est['additive_effect_lower']),
        c(res$pop$gcomp$pt_est['additive_effect'],
          res$pop$ipw$pt_est['additive_effect'],
          res$pop$aipw$pt_est['additive_effect'])
      ),
      additive_lower_ci = c(
        c(res$nat_inf$gcomp$boot_se$lower_ci_additive,
          res$nat_inf$ipw$boot_se$lower_ci_additive,
          res$nat_inf$aipw$pt_est['additive_effect'] - 1.96*res$nat_inf$aipw$pt_est['additive_se'],
          res$nat_inf$tmle$pt_est['additive_effect'] - 1.96*res$nat_inf$tmle$pt_est['additive_se'],
          res$nat_inf$bound$boot_se$lower_ci_additive_lower),
        c(res$doomed$gcomp$boot_se$lower_ci_additive,
          res$doomed$ipw$boot_se$lower_ci_additive,
          res$doomed$aipw$pt_est['additive_effect'] - 1.96*res$doomed$aipw$pt_est['additive_se'],
          res$doomed$bound$boot_se$lower_ci_additive_lower),
        c(res$pop$gcomp$boot_se$lower_ci_additive,
          res$pop$ipw$boot_se$lower_ci_additive,
          res$pop$aipw$pt_est['additive_effect'] - 1.96*res$pop$aipw$pt_est['additive_se'])
      ),
      additive_upper_ci = c(
        c(res$nat_inf$gcomp$boot_se$upper_ci_additive,
          res$nat_inf$ipw$boot_se$upper_ci_additive,
          res$nat_inf$aipw$pt_est['additive_effect'] + 1.96*res$nat_inf$aipw$pt_est['additive_se'],
          res$nat_inf$tmle$pt_est['additive_effect'] + 1.96*res$nat_inf$tmle$pt_est['additive_se'],
          res$nat_inf$bound$boot_se$upper_ci_additive_lower),
        c(res$doomed$gcomp$boot_se$upper_ci_additive,
          res$doomed$ipw$boot_se$upper_ci_additive,
          res$doomed$aipw$pt_est['additive_effect'] + 1.96*res$doomed$aipw$pt_est['additive_se'],
          res$doomed$bound$boot_se$upper_ci_additive_lower),
        c(res$pop$gcomp$boot_se$upper_ci_additive,
          res$pop$ipw$boot_se$upper_ci_additive,
          res$pop$aipw$pt_est['additive_effect'] + 1.96*res$pop$aipw$pt_est['additive_se'])
      ),
      # Multiplicative estimate and confidence interval
      mult_estimate = c(
        c(exp(as.numeric(res$nat_inf$gcomp$pt_est['log_multiplicative_effect'])),
          exp(as.numeric(res$nat_inf$ipw$pt_est['log_multiplicative_effect'])),
          exp(as.numeric(res$nat_inf$aipw$pt_est['log_multiplicative_effect'])),
          exp(as.numeric(res$nat_inf$tmle$pt_est['log_multiplicative_effect'])),
          res$nat_inf$bound$pt_est['mult_effect_lower']),
        c(exp(as.numeric(res$doomed$gcomp$pt_est['log_multiplicative_effect'])),
          exp(as.numeric(res$doomed$ipw$pt_est['log_multiplicative_effect'])),
          exp(as.numeric(res$doomed$aipw$pt_est['log_multiplicative_effect'])),
          res$doomed$bound$pt_est['mult_effect_lower']),
        c(exp(as.numeric(res$pop$gcomp$pt_est['log_multiplicative_effect'])),
          exp(as.numeric(res$pop$ipw$pt_est['log_multiplicative_effect'])),
          exp(as.numeric(res$pop$aipw$pt_est['log_multiplicative_effect'])))
      ),
      mult_lower_ci = c(
        c(res$nat_inf$gcomp$boot_se$lower_ci_mult,
          res$nat_inf$ipw$boot_se$lower_ci_mult,
          exp(res$nat_inf$aipw$pt_est['log_multiplicative_effect'] - 1.96*res$nat_inf$aipw$pt_est['log_multiplicative_se']),
          exp(res$nat_inf$tmle$pt_est['log_multiplicative_effect'] - 1.96*res$nat_inf$tmle$pt_est['log_multiplicative_se']),
          res$nat_inf$bound$boot_se$lower_ci_mult_lower),
        c(res$doomed$gcomp$boot_se$lower_ci_mult,
          res$doomed$ipw$boot_se$lower_ci_mult,
          exp(res$doomed$aipw$pt_est['log_multiplicative_effect'] - 1.96*res$doomed$aipw$pt_est['log_multiplicative_se']),
          res$doomed$bound$boot_se$lower_ci_mult_lower),
        c(res$pop$gcomp$boot_se$lower_ci_mult,
          res$pop$ipw$boot_se$lower_ci_mult,
          exp(res$pop$aipw$pt_est['log_multiplicative_effect'] - 1.96*res$pop$aipw$pt_est['log_multiplicative_se']))
      ),
      mult_upper_ci = c(
        c(res$nat_inf$gcomp$boot_se$upper_ci_mult,
          res$nat_inf$ipw$boot_se$upper_ci_mult,
          exp(res$nat_inf$aipw$pt_est['log_multiplicative_effect'] + 1.96*res$nat_inf$aipw$pt_est['log_multiplicative_se']),
          exp(res$nat_inf$tmle$pt_est['log_multiplicative_effect'] + 1.96*res$nat_inf$tmle$pt_est['log_multiplicative_se']),
          res$nat_inf$bound$boot_se$upper_ci_mult_lower),
        c(res$doomed$gcomp$boot_se$upper_ci_mult,
          res$doomed$ipw$boot_se$upper_ci_mult,
          exp(res$doomed$aipw$pt_est['log_multiplicative_effect'] + 1.96*res$doomed$aipw$pt_est['log_multiplicative_se']),
          res$doomed$bound$boot_se$upper_ci_mult_lower),
        c(res$pop$gcomp$boot_se$upper_ci_mult,
          res$pop$ipw$boot_se$upper_ci_mult,
          exp(res$pop$aipw$pt_est['log_multiplicative_effect'] + 1.96*res$pop$aipw$pt_est['log_multiplicative_se']))
      )
    )
    
    all_result_df <- rbind(all_result_df, result_df)
  }
  
}

saveRDS(all_result_df, paste0("results/contour/", setting, "_overall.Rds"))