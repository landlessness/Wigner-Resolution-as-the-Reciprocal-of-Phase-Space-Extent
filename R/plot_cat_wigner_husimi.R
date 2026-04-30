# ==============================================================================
# plot_cat_wigner_husimi.R
# Husimi resolution of the n-cat Wigner negativity for four configurations:
#
#   row 1: n=2          (lobes at N, S)
#   row 2: n=3          (apex at N, base lobes at p=-5; triangle pointing up)
#   row 3: n=4 diagonal (lobes at NE, NW, SW, SE — Zurek rotated 45 deg)
#   row 4: n=4 cardinal (lobes at N, E, S, W — Zurek original compass)
#
# Each row has three columns matching plot_compass_wigner_husimi.R:
#   Left:    bare Wigner heatmap with Husimi-circle overlay
#   Middle:  cross-section at p=0 of Husimi-convolved density (non-negative)
#   Right:   2D heatmap of Husimi-convolved density (non-negative)
#
# All rows share an extreme-p anchoring at +/- CAT_P_MAX so that the
# horizontal p=0 cross-section sits halfway between the extreme-p lobes
# in every row.
#
# When run standalone, produces figures/cat_wigner_husimi.pdf at 6.75"x9.6".
# When sourced under ASSEMBLY_MODE, builds named panels but skips ggsave.
#
# Reference: Husimi 1940 Proc. Phys. Math. Soc. Jpn. 22, 264; Lee 1995
#            Phys. Rep. 259, 147.
# ==============================================================================

library(here)
library(patchwork)

source(here("R", "plot_tools.R"))
source(here("R", "cat_system.R"))               # cat_psi, CAT_*
source(here("R", "wigner_density.R"))           # build_wigner_state,
                                                # apply_kernel_density
source(here("R", "husimi_kernel.R"))            # husimi_kernel_matrix,
                                                # husimi_overlay_layers

latex_font      <- "CMU Serif"
dir_figures     <- here("figures")
if (!dir.exists(dir_figures)) dir.create(dir_figures, recursive=TRUE)
file_output_pdf <- file.path(dir_figures, "cat_wigner_husimi.pdf")

# ------------------------------------------------------------------------------
# DISPLAY WINDOW (must match all cat_wigner_*.R scripts)
# ------------------------------------------------------------------------------

Q_DISPLAY <- 8

# ------------------------------------------------------------------------------
# DISPLAY GRIDS AND BREAKS (shared across rows)
# ------------------------------------------------------------------------------

q_lo <- -Q_DISPLAY
q_hi <-  Q_DISPLAY
p_lo <- -Q_DISPLAY
p_hi <-  Q_DISPLAY

custom_breaks_q <- c(-CAT_P_MAX, 0, CAT_P_MAX)
custom_breaks_p <- c(-CAT_P_MAX, 0, CAT_P_MAX)
label_format    <- function(x) sprintf("%.0f", x)
q_display       <- seq(q_lo, q_hi, length.out=500)
custom_breaks_s <- c(-CAT_P_MAX, 0, CAT_P_MAX)

# ------------------------------------------------------------------------------
# CONFIGURATIONS (same as plot_cat_wigner.R)
# ------------------------------------------------------------------------------

cat_configurations <- list(
  list(label="n2",  desc="n=2 (N, S)",                        n_cats=2, orientation="cardinal"),
  list(label="n3",  desc="n=3 (apex at N, base at p=-5)",     n_cats=3, orientation="cardinal"),
  list(label="n4d", desc="n=4 diagonal (NE, NW, SW, SE)",     n_cats=4, orientation="diagonal"),
  list(label="n4",  desc="n=4 cardinal (N, E, S, W = Zurek)", n_cats=4, orientation="cardinal")
)

# ------------------------------------------------------------------------------
# PER-ROW BUILDER
# ------------------------------------------------------------------------------

build_cat_husimi_row <- function(cfg) {
  cat(sprintf("\n=== %s ===\n", cfg$desc))

  cat("Building cat wavefunction...\n")
  q_psi   <- seq(CAT_Q_MIN, CAT_Q_MAX, by=CAT_DQ)
  psi_vec <- cat_psi(q_psi, n_cats=cfg$n_cats, p_max=CAT_P_MAX,
                     xi=CAT_XI, orientation=cfg$orientation)

  cat("\nBuilding Wigner on extended integration grid...\n")
  state <- build_wigner_state(psi_vec, q_psi,
                              q_lo, q_hi, p_lo, p_hi, q_display,
                              n_q_int=1601, n_p_int=1601)

  cat("\nApplying Husimi kernel...\n")
  Q_state <- apply_kernel_density(state, husimi_kernel_matrix,
                                  q_lo, q_hi, p_lo, p_hi, q_display)

  cat("\nExtracting horizontal cross-section at p=0...\n")
  dt_cross <- data.table(q=q_display, Q_husimi=Q_state$P_cross)
  cat(sprintf("  Q cross-section range: [%.4f, %.4f]  min should be >= 0\n",
              min(Q_state$P_cross), max(Q_state$P_cross)))

  q_peak <- max(dt_cross$Q_husimi, na.rm=TRUE)
  y_lim  <- q_peak * 1.1

  overlay_layers <- husimi_overlay_layers(q_center=0)
  empty_overlay  <- list()

  heatmap_overlay <- plot_wigner_heatmap(
    state$heatmap_dt, overlay_layers,
    q_lim=c(q_lo, q_hi), p_lim=c(p_lo, p_hi),
    custom_breaks_q=custom_breaks_q,
    custom_breaks_p=custom_breaks_p,
    label_format=label_format, base_font=latex_font)

  cross_panel <- plot_husimi_cross_section(
    dt_cross,
    q_lim=c(q_lo, q_hi), y_lim=y_lim,
    custom_breaks=custom_breaks_s,
    label_format=label_format,
    base_font=latex_font)

  husimi_heatmap <- plot_wigner_heatmap(
    Q_state$heatmap_dt, empty_overlay,
    q_lim=c(q_lo, q_hi), p_lim=c(p_lo, p_hi),
    custom_breaks_q=custom_breaks_q,
    custom_breaks_p=custom_breaks_p,
    label_format=label_format, base_font=latex_font)

  list(overlay_heatmap=heatmap_overlay,
       cross=cross_panel,
       husimi_heatmap=husimi_heatmap)
}

# ------------------------------------------------------------------------------
# BUILD ALL FOUR ROWS
# ------------------------------------------------------------------------------

rows <- lapply(cat_configurations, build_cat_husimi_row)
names(rows) <- vapply(cat_configurations, function(cfg) cfg$label, character(1))

# Named panels for assemble_cat.R.
cat_wigner_husimi_n2_overlay_heatmap  <- rows$n2$overlay_heatmap
cat_husimi_n2_cross                   <- rows$n2$cross
cat_husimi_n2_heatmap                 <- rows$n2$husimi_heatmap

cat_wigner_husimi_n3_overlay_heatmap  <- rows$n3$overlay_heatmap
cat_husimi_n3_cross                   <- rows$n3$cross
cat_husimi_n3_heatmap                 <- rows$n3$husimi_heatmap

cat_wigner_husimi_n4d_overlay_heatmap <- rows$n4d$overlay_heatmap
cat_husimi_n4d_cross                  <- rows$n4d$cross
cat_husimi_n4d_heatmap                <- rows$n4d$husimi_heatmap

cat_wigner_husimi_n4_overlay_heatmap  <- rows$n4$overlay_heatmap
cat_husimi_n4_cross                   <- rows$n4$cross
cat_husimi_n4_heatmap                 <- rows$n4$husimi_heatmap

# ------------------------------------------------------------------------------
# SAVE STANDALONE PDF (skipped when sourced by assemble_cat.R)
# Four rows of three panels each at 6.75"x2.4" per row -> 6.75"x9.6" total.
# ------------------------------------------------------------------------------

if (!exists("ASSEMBLY_MODE")) {
  p_full <-
    (cat_wigner_husimi_n2_overlay_heatmap |
     cat_husimi_n2_cross |
     cat_husimi_n2_heatmap) /
    (cat_wigner_husimi_n3_overlay_heatmap |
     cat_husimi_n3_cross |
     cat_husimi_n3_heatmap) /
    (cat_wigner_husimi_n4d_overlay_heatmap |
     cat_husimi_n4d_cross |
     cat_husimi_n4d_heatmap) /
    (cat_wigner_husimi_n4_overlay_heatmap |
     cat_husimi_n4_cross |
     cat_husimi_n4_heatmap)
  ggsave(filename=file_output_pdf, plot=p_full, device=cairo_pdf,
         width=6.75, height=9.6, limitsize=FALSE)
  cat("\nDone.", file_output_pdf, "\n")
}
