# ==============================================================================
# plot_tomography_sinogram.R
# Squeezed sinograms pr_delta(x, theta) for all eight states. What the
# resolution-respecting experimentalist would measure across angles.
#
# No existing plot_tools.R helper exists for sinograms (they are new in
# this Letter), but styling matches the conventions there: theme_bw,
# axis text size 8, panel grid blanked, single-direction colormap with
# HEATMAP_COLOR_LOW / HEATMAP_COLOR_HIGH stops.
# ==============================================================================

library(here)
library(data.table)
library(ggplot2)
library(patchwork)

source(here("R", "plot_tools.R"))
source(here("R", "state_builder.R"))

DATA_RDS    <- here("data", "tomography_data.rds")
COMPUTE_R   <- here("R", "tomography_compute.R")
OUTPUT_PDF  <- here("figures", "tomography_sinogram.pdf")

results <- load_tomography_cache(DATA_RDS, COMPUTE_R)

STATE_LABELS <- c(
  squeezed_vacuum = "Squeezed vacuum",
  harmonic_n1     = "Harmonic n=1",
  morse_n8        = "Morse n=8",
  double_well_n5  = "Double-well n=5",
  cat_2           = "2-cat",
  cat_3           = "3-cat",
  cat_4_square    = "4-cat square",
  cat_compass     = "Zurek compass"
)
STATE_ORDER <- names(STATE_LABELS)

latex_font <- "CMU Serif"

build_sinogram_panel <- function(r, label) {
  S      <- r$sinogram
  x      <- r$x_grid
  th_deg <- r$theta_grid * 180 / pi

  dt <- as.data.table(expand.grid(theta=th_deg, x=x))
  # Note expand.grid order: first arg varies fastest, so we transpose
  # the matrix when building the values vector.
  dt[, s := as.vector(t(S))]
  max_s <- max(dt$s, na.rm=TRUE)
  if (max_s <= 0) max_s <- 1
  dt[, s_plot := s / max_s]

  ggplot(dt, aes(x=theta, y=x, fill=s_plot)) +
    geom_raster(interpolate=TRUE) +
    scale_fill_gradient(low=HEATMAP_COLOR_LOW, high=HEATMAP_COLOR_HIGH,
                        limits=c(0, 1), guide="none") +
    coord_cartesian(expand=FALSE) +
    scale_x_continuous(breaks=c(0, 45, 90, 135),
                       labels=function(v) sprintf("%.0f", v)) +
    labs(title=label,
         x=expression(theta * "  (deg)"),
         y=expression(italic(x)[theta])) +
    theme_bw(base_family=latex_font) +
    theme(plot.title=element_text(size=10, hjust=0.5,
                                  margin=margin(b=2)),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          axis.text=element_text(size=8),
          plot.margin=margin(2, 2, 2, 2))
}

panels <- lapply(STATE_ORDER, function(s) build_sinogram_panel(results[[s]],
                                                                STATE_LABELS[[s]]))

fig <- wrap_plots(panels, ncol=4) +
  plot_annotation(theme=theme(plot.margin=margin(10,10,10,10)))

dir.create(dirname(OUTPUT_PDF), showWarnings=FALSE, recursive=TRUE)
ggsave(OUTPUT_PDF, fig,
       width=FIGURE_WIDTH_IN * 4/3,
       height=2 * ROW_HEIGHT_IN * 1.4,
       units="in", device=cairo_pdf)
cat(sprintf("Wrote %s\n", OUTPUT_PDF))
