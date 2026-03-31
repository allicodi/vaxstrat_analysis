# ----------------------------------------------------------------------------
# Make contour plots for effects 
# ----------------------------------------------------------------------------

here::i_am("R/plot_contour.R")

library(plotly)
library(RColorBrewer)
cfg <- yaml::read_yaml("config_contour.yml")

setting_names <- c("provide_immune_40_ve_66",
                   "provide_immune_60_ve_66",
                   "provide_immune_40_ve_50",
                   "provide_immune_40_ve_85")

setting_annotations <- c("P(Protected): 40%\nP(Doomed): 20%\nP(Immune): 40%\nVE: 66%",
                         "P(Protected): 27%\nP(Doomed): 13%\nP(Immune): 60%\nVE: 66%",
                         "P(Protected): 30%\nP(Doomed): 30%\nP(Immune): 40%\nVE: 50%",
                         "P(Protected): 51%\nP(Doomed):  9%\nP(Immune): 40%\nVE: 85%")

all_rows <- list()

# To compute global min/max for shared color scale
all_truth <- lapply(setting_names, function(setting) {
  readRDS(here::here(paste0("results/contour/", setting, "_truth.Rds"))) #%>%
    #filter(doomed_inflation < -0.1)
})
global_range <- range(unlist(lapply(all_truth, function(truth) {
  c(truth$effect_nat_inf, truth$effect_doomed, truth$effect_pop)
})))
levels_global <- levels_global <- pretty(global_range, n = 15)  
contour_size <- diff(levels_global)[1]
zmin <- min(levels_global)
zmax <- max(levels_global)

my_colorscale <- "YlGnBu"

# Define uniform axis settings once, outside make_contour
x_tick_vals <- seq(0, -0.3, by = -0.1)   # adjust to match your actual x range
y_tick_vals <- seq(0, -0.06, by = -0.02)   # adjust to match your actual y range
x_range <- c(max(x_tick_vals) + 0.001, min(x_tick_vals) - 0.01)  # reversed
y_range <- c(max(y_tick_vals) + 0.001, min(y_tick_vals) - 0.001)  # reversed

make_contour <- function(x, y, z_matrix, trace_name) {
  plot_ly(
    x = rev(x), y = rev(y), z = z_matrix,
    type = "contour",
    zmin = zmin, zmax = zmax,
    coloraxis = "coloraxis",
    contours = list(
      start = levels_global[1],
      end = levels_global[length(levels_global)],
      size = contour_size ,   
      showlabels = TRUE,
      labelfont = list(size = 13),
      labeldistance = 50 
    ),
    line = list(width = 2, color = 'black'),
    name = trace_name,
    showscale = FALSE
  ) %>%
    layout(
      xaxis = list(
        title = "Protected RD",
        range = x_range,
        tickvals = x_tick_vals,
        ticktext = as.character(x_tick_vals),
        tickfont = list(size = 16),
        titlefont = list(size = 20)
      ),
      yaxis = list(
        title = "Doomed RD",
        range = y_range,
        tickvals = y_tick_vals,
        ticktext = as.character(y_tick_vals),
        tickfont = list(size = 16),
        titlefont = list(size = 20)
      )
    )
}

show_legend <- TRUE

# Loop over settings
for (row_idx in seq_along(setting_names)) {
  setting <- setting_names[row_idx]
  
  truth <- all_truth[[row_idx]]
  x <- sort(unique(truth$effect_protected))
  y <- sort(unique(truth$effect_doomed))
  
  #z_nat_inf <- matrix(truth$effect_nat_inf, nrow = length(x), ncol = length(y))
  #z_doomed  <- matrix(truth$effect_doomed, nrow = length(x), ncol = length(y))
  #z_pop     <- matrix(truth$effect_pop, nrow = length(x), ncol = length(y))
  
  z_nat_inf <- matrix(truth$effect_nat_inf, nrow = length(y), ncol = length(x))
  z_doomed  <- matrix(truth$effect_doomed, nrow = length(y), ncol = length(x))
  z_pop     <- matrix(truth$effect_pop, nrow = length(y), ncol = length(x))
  
  # Load simulation results
  sim_res <- readRDS(here::here(paste0("results/contour/", setting, "_combined_contour_data.Rds"))) #%>%
    #filter(doomed_inflation < -0.1)
  combos <- expand.grid(doomed_inflation = unique(sim_res$doomed_inflation),
                        protected_inflation = unique(sim_res$protected_inflation))
  
  # make sure order did not get swapped/still using right ones
  combos$effect_protected <- truth$effect_protected[match(
    paste(truth$doomed_inflation, truth$protected_inflation),
    paste(combos$doomed_inflation, combos$protected_inflation)
  )]
  combos$effect_doomed <- truth$effect_doomed[match(
    paste(truth$doomed_inflation, truth$protected_inflation),
    paste(combos$doomed_inflation, combos$protected_inflation)
  )]

  # Compute power
  # combos$doomed_power <- NA
  combos$pop_power <- NA
  combos$nat_inf_er_power <- NA
  combos$nat_inf_cw_power <- NA
  combos$nat_inf_er_cw_power <- NA
  
  for (i in seq_len(nrow(combos))) {
    sub <- sim_res[sim_res$doomed_inflation == combos$doomed_inflation[i] &
                     sim_res$protected_inflation == combos$protected_inflation[i], ]
    # combos$doomed_power[i] <- mean(as.numeric(sub$doomed_reject), na.rm = TRUE) # NOTE some NAs for doomed reject, a handful of seeds where there are some NAs within the object for different combos of deltas, not sure why
    combos$pop_power[i]    <- mean(as.numeric(sub$pop_reject))
    combos$nat_inf_er_power[i]<- mean(as.numeric(sub$nat_inf_er_reject))
    combos$nat_inf_cw_power[i]<- mean(as.numeric(sub$nat_inf_cw_reject))
    combos$nat_inf_er_cw_power[i]<- mean(as.numeric(sub$nat_inf_er_cw_reject))
  }
  
  # making polygon to overlay for power
  get_hull <- function(df, xcol, ycol) {
    if (nrow(df) < 3) return(df[0, ])
    hull_idx <- chull(df[[xcol]], df[[ycol]])
    hull_idx <- c(hull_idx, hull_idx[1])
    df[hull_idx, ]
  }
  
  # gray out where doomed effect > protected effect
  mask_points <- expand.grid(x = seq(min(x) -0.001, 0.005, by = 0.005), 
                             y = seq(min(y) -0.001, 0.001, by = 0.005))
  mask_points <- subset(mask_points, mask_points$y < mask_points$x)
  
  # Convex hull to get the polygon boundary
  if (nrow(mask_points) > 2) {
    hull_idx <- chull(mask_points$x, mask_points$y)
    shade_poly <- mask_points[hull_idx, ]
  } else {
    shade_poly <- data.frame(x = numeric(0), y = numeric(0))
  }
  
  p_thresh <- 0.8
  # hull_doomed   <- get_hull(subset(combos, doomed_power >= p_thresh), "effect_protected", "effect_doomed")
  hull_pop      <- get_hull(subset(combos, pop_power >= p_thresh), "effect_protected", "effect_doomed")
  hull_nat_inf_er  <- get_hull(subset(combos, nat_inf_er_power >= p_thresh), "effect_protected", "effect_doomed")
  hull_nat_inf_cw  <- get_hull(subset(combos, nat_inf_cw_power >= p_thresh), "effect_protected", "effect_doomed")
  hull_nat_inf_er_cw  <- get_hull(subset(combos, nat_inf_er_cw_power >= p_thresh), "effect_protected", "effect_doomed")
  
  # fig1 <- make_contour(x, y, z_doomed, "Doomed") %>%
  #   add_polygons(
  #     data = shade_poly, x = ~x, y = ~y,
  #     fillcolor = "rgba(128,128,128,0.4)",  # semi-transparent gray
  #     line = list(width = 0),
  #     inherit = FALSE,
  #     showlegend = FALSE
  #   ) %>%
  #   add_trace(data = hull_doomed, x = ~effect_protected, y = ~effect_doomed,
  #             type = "scatter", mode = "lines", 
  #             line = list(color = "#ED0000FF", width = 5),
  #             name = "Doomed power ≥80% (n=700)", 
  #             legendgroup = "power", 
  #             showlegend = show_legend)
  
  fig2 <- make_contour(x, y, z_pop, "Population") %>%
    add_polygons(
      data = shade_poly, x = ~x, y = ~y,
      fillcolor = "rgba(128,128,128,0.4)",
      line = list(width = 0),
      inherit = FALSE,
      showlegend = FALSE
    ) %>%
    add_trace(data = hull_pop, x = ~effect_protected, y = ~effect_doomed,
              type = "scatter", mode = "lines", 
              line = list(color = "#ED0000FF", width = 5),
              name = "Population power ≥80%", 
              legendgroup = "power", 
              showlegend = show_legend)
  
  fig3 <- make_contour(x, y, z_nat_inf, "Naturally infected") %>%
    add_polygons(
      data = shade_poly, x = ~x, y = ~y,
      fillcolor = "rgba(128,128,128,0.4)",
      line = list(width = 0),
      inherit = FALSE,
      showlegend = FALSE
    ) %>%
    add_trace(data = hull_nat_inf_er, x = ~effect_protected, y = ~effect_doomed,
              type = "scatter", mode = "lines", 
              line = list(color = "#00468BFF", width = 5),
              name = "Naturally infected (ER) power ≥80%", 
              legendgroup = "power", 
              showlegend = show_legend)
  
  fig4 <- make_contour(x, y, z_nat_inf, "Naturally infected") %>%
    add_polygons(
      data = shade_poly, x = ~x, y = ~y,
      fillcolor = "rgba(128,128,128,0.4)",
      line = list(width = 0),
      inherit = FALSE,
      showlegend = FALSE
    ) %>%
    add_trace(data = hull_nat_inf_cw, x = ~effect_protected, y = ~effect_doomed,
              type = "scatter", mode = "lines", 
              line = list(color = "#925E9FFF", width = 5),
              name = "Naturally infected (PI) power ≥80%", 
              legendgroup = "power", 
              showlegend = show_legend)
  
  fig5 <- make_contour(x, y, z_nat_inf, "Naturally infected") %>%
    add_polygons(
      data = shade_poly, x = ~x, y = ~y,
      fillcolor = "rgba(128,128,128,0.4)",
      line = list(width = 0),
      inherit = FALSE,
      showlegend = FALSE
    ) %>%
    add_trace(data = hull_nat_inf_er_cw, x = ~effect_protected, y = ~effect_doomed,
              type = "scatter", mode = "lines", 
              line = list(color = "#42B540FF", width = 5),
              name = "Naturally infected (ER + PI) power ≥80%", 
              legendgroup = "power", 
              showlegend = show_legend)
  
  # row_fig <- subplot(fig1, fig2, fig3, fig4, fig5, nrows = 1,
  #                    shareX = TRUE, shareY = TRUE,
  #                    titleX = TRUE, titleY = TRUE) %>%
  #   layout(
  #     margin = list(l = 250)  # increase left margin to give space for annotations
  #   )
  # 
  # # store with row label
  # all_rows[[row_idx]] <- row_fig %>%
  #   layout(annotations = list(
  #     list(
  #       text = setting_annotations[row_idx],  # use the full annotation text
  #       x = -0.15, y = 0.5,                   # position to the left of the row
  #       xref = "paper", yref = "paper",
  #       textangle = 0,
  #       font = list(size = 14),
  #       showarrow = FALSE,
  #       xanchor = "center", yanchor = "middle"
  #     )
  #   ))
  
  fig_all <- plot_ly(
    x = x,
    y = y,
    type = "scatter",
    mode = "markers",
    marker = list(opacity = 0),
    showlegend = FALSE
  ) %>%
    layout(
      xaxis = list(
        title = "Protected RD",
        range = x_range,
        tickvals = x_tick_vals,
        ticktext = as.character(x_tick_vals),
        tickfont = list(size = 16),
        titlefont = list(size = 20)
      ),
      yaxis = list(
        title = "Doomed RD",
        range = y_range,
        tickvals = y_tick_vals,
        ticktext = as.character(y_tick_vals),
        tickfont = list(size = 16),
        titlefont = list(size = 20)
      )
    ) %>%
    add_trace(data = hull_nat_inf_er_cw,
              x = ~effect_protected, y = ~effect_doomed,
              type = "scatter", mode = "lines",
              opacity = 0.8,
              line = list(color = "#42B540FF", width = 7.5)) %>%
    add_trace(data = hull_nat_inf_cw,
              x = ~effect_protected, y = ~effect_doomed,
              type = "scatter", mode = "lines",
              opacity = 0.8,
              line = list(color = "#925E9FFF", width = 6)) %>%
    add_trace(data = hull_nat_inf_er,
              x = ~effect_protected, y = ~effect_doomed,
              type = "scatter", mode = "lines",
              opacity = 0.8,
              line = list(color = "#00468BFF", width = 6)) %>%
    add_trace(data = hull_pop,
              x = ~effect_protected, y = ~effect_doomed,
              type = "scatter", mode = "lines",
              opacity = 0.8,
              line = list(color = "#ED0000FF", width = 4.8)) 
  
  # all_rows[[length(all_rows) + 1]] <- fig1
  all_rows[[length(all_rows) + 1]] <- fig2
  all_rows[[length(all_rows) + 1]] <- fig3
  all_rows[[length(all_rows) + 1]] <- fig4
  all_rows[[length(all_rows) + 1]] <- fig5
  all_rows[[length(all_rows) + 1]] <- fig_all
  
  
  show_legend <- FALSE
}


# ------------------------------------------------------------------------------
# Assemble final 4 x 5 plot (NO nested subplots)
# ------------------------------------------------------------------------------

n_rows <- length(setting_names)
n_cols <- 5


final_fig <- do.call(
  subplot,
  c(
    all_rows,
    list(
      nrows = n_rows,
      shareX = TRUE,
      shareY = TRUE,
      titleX = TRUE,
      titleY = TRUE
    )
  )
)


base_panel_height <- 250
extra_vertical_space <- 350  # space for top + bottom annotations

final_fig <- final_fig %>%
  layout(
    height = n_rows * base_panel_height + extra_vertical_space
  )

# final_fig <- final_fig %>%
#   layout(
#     margin = list(l = 0, r = 0, t = 10, b = 10),  # control spacing between plots
#     height = n_rows * 250  # adjust height for multiple rows
#   )
# final_fig

# ------------------------------------------------------------------------------
# Row annotations (left side, centered per row)
# ------------------------------------------------------------------------------

row_ys <- seq(
  from = 1 - 1 / (2 * n_rows),
  to   = 1 / (2 * n_rows),
  length.out = n_rows
)

row_annotations <- lapply(seq_len(n_rows), function(i) {
  list(
    text = setting_annotations[i],
    x = -0.06,
    y = row_ys[i],
    xref = "paper",
    yref = "paper",
    showarrow = FALSE,
    xanchor = "right",
    yanchor = "middle",
    align = "right",
    font = list(size = 18)
  )
})

# ------------------------------------------------------------------------------
# Column titles (top, evenly spaced)
# ------------------------------------------------------------------------------

col_titles <- c(
  "Population\n ",
  "Naturally infected\n(ER)",
  "Naturally infected\n(PI)",
  "Naturally infected\n(ER + PI)",
  "All Power\n "
)

col_xs <- seq(
  from = 1 / (2 * n_cols),
  to   = 1 - 1 / (2 * n_cols),
  length.out = n_cols
)

col_annotations <- lapply(seq_len(n_cols), function(i) {
  list(
    text = col_titles[i],
    x = col_xs[i],
    y = 1,
    xref = "paper",
    yref = "paper",
    showarrow = FALSE,
    xanchor = "center",
    yanchor = "bottom",
    font = list(size = 20)
  )
})

# ------------------------------------------------------------------------------
# Final layout (colorbar, legend, annotations)
# ------------------------------------------------------------------------------

final_fig <- final_fig %>%
  layout(
    coloraxis = list(
      colorscale = my_colorscale,
      cmin = zmin,
      cmax = zmax,
      colorbar = list(
        title = list(text = "Effect size", side = "top"),
        tickfont = list(size = 16)
      )
    ),
    legend = list(
      orientation = "h",
      x = 0.5,
      y = -0.1,
      xanchor = "center",
      yanchor = "top",
      font = list(size = 17),
      bgcolor = "rgba(255,255,255,0)",
      traceorder = "normal",       # <-- add this
      itemwidth = -1,              # <-- controls spacing between items
      tracegroupgap = 0            # <-- collapse group gaps
    ),
    annotations = c(row_annotations, col_annotations),
    margin = list(l = 260, t = 100, b = 240)  # increase bottom margin to fit 2 rows THIS IS WHAT IS BREAKING MY CONTOURS
  )

final_fig

#final_fig <- plotly_build(final_fig)
#final_fig$x$layout$margin <- list(l = 200, t = 120, b = 240)

# Save without margin
htmlwidgets::saveWidget(final_fig, here::here("results/contour/figures/contour_VE.html"), selfcontained = TRUE)

# Then take the screenshot with webshot padding instead
webshot2::webshot(
  here::here("results/contour/figures/contour_VE.html"),
  here::here("results/contour/figures/contour_VE.png"),
  vwidth = 1800, vheight = 1300,
  delay = 2
)

# orca(final_fig, file = here::here("results/contour/figures/contour_VE.png"), 
#      format = tools::file_ext(file),
#      scale = NULL, width = 1200, height = 1000, mathjax = FALSE,
#      parallel_limit = NULL, verbose = FALSE, debug = FALSE,
#      safe = FALSE, more_args = NULL)

#plotly::plotly_IMAGE(final_fig, width = 1200, height = 1000, format = "png", here::here("results/contour/figures/contour_VE.png"))
