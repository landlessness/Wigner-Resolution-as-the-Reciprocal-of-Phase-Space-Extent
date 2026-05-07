# ==============================================================================
# plot_tomography_difference.R
# Difference W(q, p) - tilde_W_delta(q, p) as a signed heatmap. Should
# light up exactly where W goes negative (cat states, Fock). Uses
# plot_wigner_heatmap() (diverging colormap) for consistency with the
# left column of plot_wigner.R / plot_cats.R.
# ==============================================================================

library(here)
library(data.table)
library(ggplot2)
library(patchwork)

source(here("R", "plot_tools.R"))
source(here("R", "state_builder.R"))

DATA_RDS    <- here("data", "tomography_data.rds")
COMPUTE_R   <- here("R", "tomography_compute.R")
OUTPUT_PDF  <- here("figures", "tomography_difference.pdf")

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

build_signed_heatmap_dt <- function(W, q_grid, p_grid) {
  dt <- as.data.table(expand.grid(q=q_grid, p=p_grid))
  dt[, w := as.vector(W)]
  max_abs <- max(abs(dt$w), na.rm=TRUE)
  if (max_abs <= 0) max_abs <- 1
  dt[, w_plot := w / max_abs]
  dt
}

build_panel <- function(r, label) {
  D     <- r$W_tomo - r$tilde_W
  hm_dt <- build_signed_heatmap_dt(D, r$q_grid_tomo, r$p_grid_tomo)
  max_D <- max(abs(D), na.rm=TRUE)

  p <- plot_wigner_heatmap(
    hm_dt, r$overlay_layers, df_traj=NULL,
    q_lim=c(r$q_lo, r$q_hi), p_lim=c(r$p_lo, r$p_hi),
    custom_breaks_q=r$custom_breaks_q,
    custom_breaks_p=r$custom_breaks_p,
    label_format=r$label_format, base_font=latex_font)
  p + ggtitle(sprintf("%s (max|D|=%.2g)", label, max_D)) +
    theme(plot.title=element_text(size=10, hjust=0.5,
                                  family=latex_font,
                                  margin=margin(b=2)))
}

panels <- lapply(STATE_ORDER, function(s) build_panel(results[[s]],
                                                       STATE_LABELS[[s]]))

fig <- wrap_plots(panels, ncol=4) +
  plot_annotation(theme=theme(plot.margin=margin(10,10,10,10)))

dir.create(dirname(OUTPUT_PDF), showWarnings=FALSE, recursive=TRUE)
ggsave(OUTPUT_PDF, fig,
       width=FIGURE_WIDTH_IN * 4/3,
       height=2 * ROW_HEIGHT_IN * 1.4,
       units="in", device=cairo_pdf)
cat(sprintf("Wrote %s\n", OUTPUT_PDF))
