# ---------------------------------------------------------------------------
# Script to aggregate results from project space into single csv file
# ---------------------------------------------------------------------------

.libPaths(c("/apps/R/4.4.0/lib64/R/site/library","/apps/R/4.4.0/lib64/R/library", "~/Rlibs_ve_trial"))

here::i_am("aggregate_results.R")

library(dplyr)
library(readr)

# Path to projects folder where results will be saved
project_dir <- "/projects/dbenkes/allison/vegrowth_analysis/results/"

# setting to look for name
setting <- Sys.getenv("SETTING")

# List all matching .rds files
files <- list.files(
  path = project_dir,
  pattern = paste0("^", setting, "_overall_seed_.*\\.rds$"),
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No files found for setting: ", setting)
}

# Read and combine all files
result_df <- files %>%
  lapply(readRDS) %>%
  bind_rows()

# Save aggregated file as CSV
out_file <- file.path("results", paste0(setting, "_results.csv"))
write_csv(result_df, out_file)
