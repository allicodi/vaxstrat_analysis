# ---------------------------------------------------------------------------
# Script to plot power for each combination of settings
# ---------------------------------------------------------------------------

# read setting from config
# get each combo of parameters (not n)
# for each combo at each n, calc power, make df
# use df to make power plot

library(dplyr)
library(scales) 
library(ggplot2)

here::i_am("plot_power.R")

setting <- "vary_inflation"

results_df <- read.csv(here::here(paste0("results/", setting, "_results.csv")))

# effect_protect = effect in the protected (FALSE no effect, TRUE effect)
# inflation = effect in doomed (0 no effect, -1 effect)
combo_grid <- expand.grid(
  effect_protect = unique(results_df$effect_protect),
  inflation = unique(results_df$inflation)
)

results_grid <- expand.grid(
  effect_protect = unique(results_df$effect_protect),
  inflation = unique(results_df$inflation),
  n = unique(results_df$n_sample_size)
)

results_grid$nat_inf_power <- NA
results_grid$doomed_power <- NA
results_grid$pop_power <- NA

for(i in 1:nrow(combo_grid)){
  sub_res <- results_df[
    results_df$effect_protect == combo_grid$effect_protect[i] &
      results_df$inflation == combo_grid$inflation[i], ]
  
  for(n in unique(sub_res$n_sample_size)){
    sub_res_n <- sub_res[sub_res$n_sample_size == n, ]
    
    results_grid$nat_inf_power[results_grid$n == n &
                                 results_grid$effect_protect == combo_grid$effect_protect[i] &
                                 results_grid$inflation == combo_grid$inflation[i]] <-
      mean(as.numeric(sub_res_n$nat_inf_reject == TRUE))
    
    results_grid$nat_inf_bound_power[results_grid$n == n &
                                 results_grid$effect_protect == combo_grid$effect_protect[i] &
                                 results_grid$inflation == combo_grid$inflation[i]] <-
      mean(as.numeric(sub_res_n$nat_inf_upper_bound_reject == TRUE), na.rm = TRUE) # some NA with smaller sample size?? 
    
    results_grid$doomed_power[results_grid$n == n &
                                results_grid$effect_protect == combo_grid$effect_protect[i] &
                                results_grid$inflation == combo_grid$inflation[i]] <-
      mean(as.numeric(sub_res_n$doomed_reject == TRUE))
    
    results_grid$doomed_bound_power[results_grid$n == n &
                                results_grid$effect_protect == combo_grid$effect_protect[i] &
                                results_grid$inflation == combo_grid$inflation[i]] <-
      mean(as.numeric(sub_res_n$doomed_upper_bound_reject == TRUE), na.rm = TRUE) # some NA with smaller sample size?? 
    
    results_grid$pop_power[results_grid$n == n &
                             results_grid$effect_protect == combo_grid$effect_protect[i] &
                             results_grid$inflation == combo_grid$inflation[i]] <-
      mean(as.numeric(sub_res_n$pop_reject == TRUE))
  }
}

# Short-hand labels for facets
results_grid <- results_grid %>%
  mutate(
    effect_protect = if_else(effect_protect == TRUE, "Protected effect < 0", "Protected effect = 0"),
    effect_doomed  = if_else(inflation == 0, "Doomed effect = 0", if_else(inflation == -0.5, "Doomed effect < 0\n(Inflation = -0.5)", "Doomed effect < 0\n(Inflation = -1)"))
  )

results_grid$effect_doomed <- factor(results_grid$effect_doomed, levels = c("Doomed effect < 0\n(Inflation = -1)",
                                                                            "Doomed effect < 0\n(Inflation = -0.5)",
                                                                            "Doomed effect = 0"))

# Plot power in percent with Lancet color scheme
ggplot(results_grid, aes(x = n, y = nat_inf_power * 100, color = "Naturally Infected")) +
  geom_line(size = 1.2) +
  geom_line(aes(y = doomed_power * 100, color = "Doomed"), size = 1.2) +
  geom_line(aes(y = pop_power * 100, color = "Population"), size = 1.2) +
  geom_line(aes(y = nat_inf_bound_power * 100, color = "Naturally Infected"), linetype = 3, size = 1.2) + 
  geom_line(aes(y = doomed_bound_power * 100, color = "Doomed"), linetype = 2, size = 1.2) + 
  facet_grid(effect_protect ~ effect_doomed) +
  scale_color_manual(
    name = "Estimand",
    values = c(
      "Naturally Infected" = "#42B540FF", 
      "Doomed" = "#ED0000FF",          
      "Population" = "#00468BFF"        
    )
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    x = "Sample Size (n)",
    y = "Power (%)",
    title = "Provide-like simulations power by sample size and effect settings"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold")
  )
# ^ only show protected > doomed and add one for protected

# contour plot
# fix sample size
# settings w combinations of doomed and protected effect sizes
# throw out combos where doomed > protected
# plot power
# contours are the size of the effect in the naturally effected 
# shade regions/add contours/outline or something where there's 80% power

# repeat w same axes, contours are the size in the population estimand
# repeat for doomed

# fourth panel - where there's 80% power, overlay all on single panel

