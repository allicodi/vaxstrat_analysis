# ----------------------------------------------------------------------------
# Make contour plots for effects 
# ----------------------------------------------------------------------------

here::i_am("plot_contour.R")

source(here::here("get_truth.R"))

cfg <- yaml::read_yaml("config_contour.yml")
config <- cfg[["contour_plot"]]

results <- readRDS(here::here("results/contour/contour_plot_overall.Rds"))



