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

setting_names <- c("provide_immune_30_ve_66__2",
                   "provide_immune_40_ve_66__2",
                   "provide_immune_50_ve_66__2",
                   "provide_immune_60_ve_66__2",
                   "provide_immune_70_ve_66__2",
                   "provide_immune_80_ve_66__2")#,
                   #"provide_immune_70_ve_50",
                   #"provide_immune_70_ve_85")

setting_annotations <- c("Protected: 47%\nDoomed: 23%\nImmune: 30%\nVE: 66%",
                         "Protected: 33%\nDoomed: 17%\nImmune: 50%\nVE: 66%",
                         "Protected: 20%\nDoomed: 10%\nImmune: 70%\nVE: 66%")#,
                         #"Protected: 15%\nDoomed: 15%\nImmune: 70%\nVE: 50%",
                         #"Protected: 25%\nDoomed: 5%\nImmune: 70%\nVE: 85%")

all_rows <- list()

# To compute global min/max for shared color scale
all_truth <- lapply(setting_names, function(setting) {
  readRDS(here::here(paste0("results/contour/", setting, "_truth.Rds")))
})
global_range <- range(unlist(lapply(all_truth, function(truth) {
  c(truth$effect_nat_inf, truth$effect_doomed, truth$effect_pop)
})))
levels_global <- pretty(global_range, n = 30)
contour_size <- diff(levels_global)[1]
zmin <- min(levels_global)
zmax <- max(levels_global)

my_colorscale <- "YlGnBu"

make_contour <- function(x, y, z_matrix, trace_name) {
  plot_ly(
    x = rev(x), y = rev(y), z = z_matrix,
    type = "contour",
    zmin = zmin, zmax = zmax,
    coloraxis = "coloraxis",  # <- use shared coloraxis
    contours = list(
      start = levels_global[1],
      end = levels_global[length(levels_global)],
      size = contour_size,
      showlabels = TRUE,
      labelfont = list(size = 14)
    ),
    line = list(width = 2, color = 'black'),
    name = trace_name,
    showscale = FALSE  # turn off per-trace colorbars
  ) %>%
    layout(
      xaxis = list(title = "Protected effect", autorange = "reversed"),
      yaxis = list(title = "Doomed effect", autorange = "reversed")
    )
}

show_legend <- TRUE

# Loop over settings
for (row_idx in seq_along(setting_names)) {
  setting <- setting_names[row_idx]
  
  truth <- all_truth[[row_idx]]
  x <- sort(unique(truth$effect_protected))
  y <- sort(unique(truth$effect_doomed))
  
  z_nat_inf <- matrix(truth$effect_nat_inf, nrow = length(x), ncol = length(y))
  z_doomed  <- matrix(truth$effect_doomed, nrow = length(x), ncol = length(y))
  z_pop     <- matrix(truth$effect_pop, nrow = length(x), ncol = length(y))
  
  # Load simulation results
  sim_res <- readRDS(here::here(paste0("results/contour/", setting, "_combined_contour_data.Rds")))
  combos <- expand.grid(doomed_inflation = unique(sim_res$doomed_inflation),
                        protected_inflation = unique(sim_res$protected_inflation))
  combos$effect_protected <- truth$effect_protected[match(
    paste(truth$doomed_inflation, truth$protected_inflation),
    paste(combos$doomed_inflation, combos$protected_inflation)
  )]
  combos$effect_doomed <- truth$effect_doomed[match(
    paste(truth$doomed_inflation, truth$protected_inflation),
    paste(combos$doomed_inflation, combos$protected_inflation)
  )]
  
  # Compute power
  combos$doomed_power <- NA
  combos$pop_power <- NA
  combos$nat_inf_power <- NA
  for (i in seq_len(nrow(combos))) {
    sub <- sim_res[sim_res$doomed_inflation == combos$doomed_inflation[i] &
                     sim_res$protected_inflation == combos$protected_inflation[i], ]
    combos$doomed_power[i] <- mean(as.numeric(sub$doomed_reject)) # NOTE some NAs for doomed reject, a handful of seeds where there are some NAs within the object for different combos of deltas, not sure why
    combos$pop_power[i]    <- mean(as.numeric(sub$pop_reject))
    combos$nat_inf_power[i]<- mean(as.numeric(sub$nat_inf_reject))
  }
  
  get_hull <- function(df, xcol, ycol) {
    if (nrow(df) < 3) return(df[0, ])
    hull_idx <- chull(df[[xcol]], df[[ycol]])
    hull_idx <- c(hull_idx, hull_idx[1])
    df[hull_idx, ]
  }
  
  p_thresh <- 0.8
  hull_doomed   <- get_hull(subset(combos, doomed_power >= p_thresh), "effect_protected", "effect_doomed")
  hull_pop      <- get_hull(subset(combos, pop_power >= p_thresh), "effect_protected", "effect_doomed")
  hull_nat_inf  <- get_hull(subset(combos, nat_inf_power >= p_thresh), "effect_protected", "effect_doomed")
  
  fig1 <- make_contour(x, y, z_doomed, "Doomed") %>%
    add_trace(data = hull_doomed, x = ~effect_protected, y = ~effect_doomed,
              type = "scatter", mode = "lines", 
              line = list(color = "#ED0000FF", width = 5),
              name = "Doomed power ≥80% (n=700)", 
              legendgroup = "power", 
              showlegend = show_legend)
  
  fig2 <- make_contour(x, y, z_pop, "Population") %>%
    add_trace(data = hull_pop, x = ~effect_protected, y = ~effect_doomed,
              type = "scatter", mode = "lines", 
              line = list(color = "#00468BFF", width = 5),
              name = "Population power ≥80% (n=700)", 
              legendgroup = "power", 
              showlegend = show_legend)
  
  fig3 <- make_contour(x, y, z_nat_inf, "Naturally infected") %>%
    add_trace(data = hull_nat_inf, x = ~effect_protected, y = ~effect_doomed,
              type = "scatter", mode = "lines", 
              line = list(color = "#42B540FF", width = 5),
              name = "Naturally infected power ≥80% (n=700)", 
              legendgroup = "power", 
              showlegend = show_legend)
  
  row_fig <- subplot(fig1, fig2, fig3, nrows = 1,
                     shareX = TRUE, shareY = TRUE,
                     titleX = TRUE, titleY = TRUE) %>%
    layout(
      margin = list(l = 250)  # increase left margin to give space for annotations
    )
  
  # store with row label
  all_rows[[row_idx]] <- row_fig %>%
    layout(annotations = list(
      list(
        text = setting_annotations[row_idx],  # use the full annotation text
        x = -0.2, y = 0.5,                   # position to the left of the row
        xref = "paper", yref = "paper",
        textangle = 0,
        font = list(size = 14),
        showarrow = FALSE,
        xanchor = "center", yanchor = "middle"
      )
    ))
  
  show_legend <- FALSE
}

final_fig <- subplot(all_rows, nrows = length(setting_names), shareX = TRUE, shareY = TRUE) %>%
  layout(
    coloraxis = list(
      colorscale = my_colorscale,
      cmin = zmin,
      cmax = zmax,
      colorbar = list(
        title = list(text = "Effect size", side = "top"),
        tickfont = list(size = 14)
      )
    ),
    legend = list(
      orientation = "h",
      x = 1.025,
      y = 0,
      xanchor = "left",
      yanchor = "top",
      font = list(size = 14),
      bgcolor = "rgba(255,255,255,0)"
    )
  )


final_fig

