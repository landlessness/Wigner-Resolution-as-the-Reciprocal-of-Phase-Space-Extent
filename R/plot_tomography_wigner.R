# ==============================================================================
# plot_tomography_wigner.R
# Bare Wigner W(q, p) for all eight states. Sanity check that the input
# to the tomography pipeline is correct. Uses plot_wigner_heatmap() from
# plot_tools.R so styling matches plot_wigner.R / plot_cats.R exactly.
# ==============================================================================

library(here)
library(data.table)
library(ggplot2)
library(patchwork)

source(here("R", "plot_tools.R"))
source(here("R", "state_builder.R"))   # load_tomography_cache

DATA_RDS    <- here("data", "tomography_data.rds")
COMPUTE_R   <- here("R", "tomography_compute.R")
OUTPUT_PDF  <- here("figures", "tomography_wigner.pdf")

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

# ------------------------------------------------------------------------------
# BUILD ONE PANEL
#
# The state bundle in `r` already carries:
#   r$state$heatmap_dt    -- the data.table plot_wigner_heatmap expects
#   r$overlay_layers      -- symplectic overlay (Fermi blob + conjugate blobs)
#   r$q_lo, r$q_hi, ...    -- display window
#   r$custom_breaks_q, .._p, r$label_format
# So we just call plot_wigner_heatmap directly.
# ------------------------------------------------------------------------------

build_panel <- function(r, label) {
  p <- plot_wigner_heatmap(
    r$state$heatmap_dt, r$overlay_layers, df_traj=NULL,
    q_lim=c(r$q_lo, r$q_hi), p_lim=c(r$p_lo, r$p_hi),
    custom_breaks_q=r$custom_breaks_q,
    custom_breaks_p=r$custom_breaks_p,
    label_format=r$label_format, base_font=latex_font)
  p + ggtitle(label) +
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
