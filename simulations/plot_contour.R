# ----------------------------------------------------------------------------
# Make contour plots for effects 
# ----------------------------------------------------------------------------

# todo
# keep ratio protected to doomed same
# another row with more immune
# third row with better VE (this VE is ~50%) -- stick w larger immune and higher VE

here::i_am("plot_contour.R")

library(plotly)
library(RColorBrewer)

cfg <- yaml::read_yaml("config_contour.yml")
config <- cfg[["provide_contour_plot"]]

truth <- readRDS(here::here("results/contour/provide_contour_plot_truth.Rds"))

# Plot contours and effect sizes from truth code ------------------------------

#x <- sort(unique(truth$protected_inflation))
#y <- sort(unique(truth$doomed_inflation))

x <- sort(unique(truth$effect_protected))
y <- sort(unique(truth$effect_doomed))

z_nat_inf <- matrix(truth$effect_nat_inf, nrow = length(x), ncol = length(y))
z_doomed <- matrix(truth$effect_doomed, nrow = length(x), ncol = length(y))
z_pop <- matrix(truth$effect_pop, nrow = length(x), ncol = length(y))

# Define exact contour levels
levels <- pretty(c(truth$effect_nat_inf, truth$effect_doomed, truth$effect_pop), n = 25)
zmin <- min(levels)
zmax <- max(levels)
contour_size <- diff(levels)[1]

# Choose a colorscale and opacity
my_colorscale <- "YlGnBu"  # can replace with "Plasma", "Cividis", or diverging palettes
my_opacity <- 0.8

# Helper function for a single contour plot with trace name
make_contour <- function(z_matrix, trace_name, show_colorbar = TRUE) {
  plot_ly(
    x = rev(x),
    y = rev(y),
    z = z_matrix,
    type = "contour",
    name = trace_name,  # set trace name,
    colorscale = my_colorscale,
    opacity = my_opacity,
    contours = list(
      start = levels[1],
      end = levels[length(levels)],
      size = contour_size,
      showlabels = TRUE,
      labelfont = list(size = 14)
    ),
    line = list(width = 2, color= 'black'),
    showscale = show_colorbar,
  ) %>%
    layout(
      xaxis = list(title = "Protected effect", autorange = "reversed"),
      yaxis = list(title = "Doomed effect", autorange = "reversed")
    ) %>%
    colorbar(title = "Effect\nsize")  %>%
    add_annotations(
      text = trace_name,
      x = 0.5,
      y = 1.025,
      yref = "paper",
      xref = "paper",
      xanchor = "center",
      yanchor = "top",
      showarrow = FALSE,
      font = list(size = 20)
    )
  
}

# Create individual plots with trace names

# fig1 <- make_contour(z_doomed, trace_name = "Doomed", show_colorbar = FALSE)
# fig2 <- make_contour(z_pop, trace_name = "Population", show_colorbar = FALSE)
# fig3 <- make_contour(z_nat_inf, trace_name = "Naturally infected", show_colorbar = TRUE)
# 
# # Combine subplots
# fig <- subplot(fig1, fig2, fig3,
#                nrows = 1,
#                shareX = TRUE, shareY = TRUE,
#                titleX = TRUE, titleY = TRUE) %>%
#   layout(
#     # Shared color scale
#     coloraxis = list(cmin = zmin, cmax = zmax)
#   ) %>%
#   colorbar(title = "Effect size")
# 
# fig


# Overlay power from AIPW at sample size n = 700 ------------------------------

sim_res <- readRDS(here::here("results/contour/combined_provide_contour_plot.Rds"))

combos <- expand.grid(doomed_inflation = unique(sim_res$doomed_inflation),
                      protected_inflation = unique(sim_res$protected_inflation))

combos$effect_protected <- truth$effect_protected[truth$doomed_inflation == combos$doomed_inflation &
                                                    truth$protected_inflation == combos$protected_inflation]
combos$effect_doomed <- truth$effect_doomed[truth$doomed_inflation == combos$doomed_inflation &
                                                    truth$protected_inflation == combos$protected_inflation]

combos$doomed_power <- NA
combos$pop_power <- NA
combos$nat_inf_power <- NA

for(i in 1:nrow(combos)){
  sub <- sim_res[sim_res$doomed_inflation == combos$doomed_inflation[i] &
                   sim_res$protected_inflation == combos$protected_inflation[i],]
  
  combos$doomed_power[i] <- mean(as.numeric(sub$doomed_reject))
  combos$pop_power[i] <- mean(as.numeric(sub$pop_reject))
  combos$nat_inf_power[i] <- mean(as.numeric(sub$nat_inf_reject))
}

# Helper to extract convex hull coordinates
get_hull <- function(df, xcol, ycol) {
  hull_idx <- chull(df[[xcol]], df[[ycol]])
  hull_idx <- c(hull_idx, hull_idx[1])  # close the polygon
  df[hull_idx, ]
}

# threshold
p_thresh <- 0.8

# Get grid points where power > 0.8 for each estimand
pts_doomed <- combos[combos$doomed_power >= p_thresh, ]
pts_pop <- combos[combos$pop_power >= p_thresh, ]
pts_nat_inf <- combos[combos$nat_inf_power >= p_thresh, ]

# Build convex hull polygons
hull_doomed   <- get_hull(pts_doomed, "effect_protected", "effect_doomed")
hull_pop      <- get_hull(pts_pop, "effect_protected", "effect_doomed")
hull_nat_inf  <- get_hull(pts_nat_inf, "effect_protected", "effect_doomed")

fig1 <- make_contour(z_doomed, trace_name = "Doomed", show_colorbar = FALSE) %>%
  add_trace(
    data = hull_doomed,
    x = ~effect_protected,
    y = ~effect_doomed,
    type = "scatter",
    mode = "lines",
    fill = "toself",
    fillcolor = "rgba(255, 0, 0, 0)",   
    line = list(color = "#ED0000FF", width = 5),
    inherit = FALSE,
    name = "Doomed\nPower ≥ 80%\n(n=700)",
    showlegend = TRUE
  )

fig2 <- make_contour(z_pop, trace_name = "Population", show_colorbar = FALSE) %>%
  add_trace(
    data = hull_pop,
    x = ~effect_protected,
    y = ~effect_doomed,
    type = "scatter",
    mode = "lines",
    fill = "toself",
    fillcolor = "rgba(255, 0, 0, 0)",   
    line = list(color = "#00468BFF", width = 5),
    inherit = FALSE,
    name = "Population\nPower ≥ 80%\n(n=700)",
    showlegend = TRUE
  )

fig3 <- make_contour(z_nat_inf, trace_name = "Naturally infected", show_colorbar = TRUE) %>%
  add_trace(
    data = hull_nat_inf,
    x = ~effect_protected,
    y = ~effect_doomed,
    type = "scatter",
    mode = "lines",
    fill = "toself",
    fillcolor = "rgba(255, 0, 0, 0)",   
    line = list(color = "#42B540FF", width = 5),
    inherit = FALSE,
    name = "Naturally Infected\nPower ≥ 80%\n(n=700)",
    showlegend = TRUE
  )

# Combine them again
fig <- subplot(fig1, fig2, fig3,
               nrows = 1,
               shareX = TRUE, shareY = TRUE,
               titleX = TRUE, titleY = TRUE) %>%
  layout(
    coloraxis = list(cmin = zmin, cmax = zmax)
  )

fig
