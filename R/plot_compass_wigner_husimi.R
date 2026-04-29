# ==============================================================================
# plot_compass_wigner_husimi.R
# Compass-state Husimi resolution panel — middle row of paper Fig. 3.
#
# When run standalone, produces a 3-panel PDF showing the full convolution
# story for the Husimi kernel:
#
#   compass_wigner_husimi_overlay_heatmap | compass_husimi_cross | compass_husimi_heatmap
#
#   Left:    Pristine Wigner W(q, p) with the Husimi kernel contour
#            (single unit circle) drawn on top, showing the kernel
#            footprint against the unsmoothed input.
#   Center:  Cross-section Q along the SW-NE diagonal, post-convolution.
#   Right:   Convolved Husimi density Q(q, p) = (W * G_husimi)(q, p), with
#            the same Husimi kernel contour overlay for direct comparison
#            with the left panel.
#
# When sourced by assemble_compass.R (which sets ASSEMBLY_MODE), all three
# panels remain available as named variables but no PDF is written.
#
# Reference: Husimi 1940 Proc. Phys.-Math. Soc. Japan 22, 264;
#            Lee 1995 Phys. Rep. 259, 147.
# ==============================================================================

library(here)
library(patchwork)

source(here("R", "plot_tools.R"))
source(here("R", "compass_system.R"))           # compass_psi, COMPASS_*
source(here("R", "wigner_density.R"))           # build_wigner_state,
                                                # apply_kernel_density
source(here("R", "husimi_kernel.R"))            # husimi_kernel_matrix,
                                                # husimi_overlay_layers

latex_font      <- "CMU Serif"
dir_figures     <- here("figures")
if (!dir.exists(dir_figures)) dir.create(dir_figures, recursive=TRUE)
file_output_pdf <- file.path(dir_figures, "compass_wigner_husimi.pdf")

# ------------------------------------------------------------------------------
# DISPLAY WINDOW (must match all compass_wigner_*.R scripts)
# ------------------------------------------------------------------------------

Q_DISPLAY <- 8

# ------------------------------------------------------------------------------
# COMPASS WAVEFUNCTION
# ------------------------------------------------------------------------------

cat("Building compass state...\n")
q_psi   <- seq(COMPASS_Q_MIN, COMPASS_Q_MAX, by=COMPASS_DQ)
psi_vec <- compass_psi(q_psi, L=COMPASS_L, xi=COMPASS_XI)
cat(sprintf("  L=%.1f  xi=%.1f  nq=%d  dq=%.4f\n",
            COMPASS_L, COMPASS_XI, length(q_psi), COMPASS_DQ))

# ------------------------------------------------------------------------------
# DISPLAY GRIDS AND BREAKS
# ------------------------------------------------------------------------------

q_lo <- -Q_DISPLAY
q_hi <-  Q_DISPLAY
p_lo <- -Q_DISPLAY
p_hi <-  Q_DISPLAY

custom_breaks_q <- c(-COMPASS_L/2, 0, COMPASS_L/2)
custom_breaks_p <- c(-COMPASS_L/2, 0, COMPASS_L/2)
label_format    <- function(x) sprintf("%.0f", x)
q_display       <- seq(q_lo, q_hi, length.out=500)

s_max     <- Q_DISPLAY
s_display <- seq(-s_max, s_max, length.out=1000)
custom_breaks_s <- c(-COMPASS_L/2, 0, COMPASS_L/2)

# ------------------------------------------------------------------------------
# BUILD WIGNER STATE
# ------------------------------------------------------------------------------

cat("\nBuilding Wigner on extended integration grid...\n")
state <- build_wigner_state(psi_vec, q_psi,
                            q_lo, q_hi, p_lo, p_hi, q_display,
                            n_q_int=1601, n_p_int=1601)

# ------------------------------------------------------------------------------
# APPLY THE HUSIMI KERNEL
# ------------------------------------------------------------------------------

cat("\nApplying Husimi kernel...\n")
husimi_for_state <- function(qg, pg) husimi_kernel_matrix(qg, pg)
Q <- apply_kernel_density(state, husimi_for_state,
                          q_lo, q_hi, p_lo, p_hi, q_display)

# ------------------------------------------------------------------------------
# DIAGONAL CROSS-SECTION (post-convolution)
# ------------------------------------------------------------------------------

cat("\nExtracting SW-NE diagonal cross-section...\n")
Q_diag <- extract_diagonal_cross_section(Q$P_matrix,
                                         state$q_int, state$p_int,
                                         s_display)
dt_cross <- data.table(q=s_display, Q_husimi=Q_diag)
cat(sprintf("  Q diagonal cross-section range: [%.4f, %.4f]  min should be >= 0\n",
            min(Q_diag), max(Q_diag)))

q_peak <- max(dt_cross$Q_husimi, na.rm=TRUE)
y_lim  <- q_peak * 1.1

# ------------------------------------------------------------------------------
# OVERLAY: Husimi kernel contour (single solid black unit circle).
# Drawn on the column-1 panel (kernel sitting on the input Wigner) but NOT
# on the column-3 panel (the overlay would be visually redundant on the
# already-convolved output).
# ------------------------------------------------------------------------------

overlay_layers <- husimi_overlay_layers(q_center=0)
empty_overlay  <- list()

# ------------------------------------------------------------------------------
# BUILD NAMED PANELS
# Three panels for the standalone 3-column layout; assemble_compass.R
# uses all three.
# ------------------------------------------------------------------------------

compass_wigner_husimi_overlay_heatmap <- plot_wigner_heatmap(
  state$heatmap_dt, overlay_layers,
  q_lim=c(q_lo, q_hi), p_lim=c(p_lo, p_hi),
  custom_breaks_q=custom_breaks_q,
  custom_breaks_p=custom_breaks_p,
  label_format=label_format, base_font=latex_font)

compass_husimi_cross <- plot_husimi_cross_section(
  dt_cross,
  q_lim=c(-s_max, s_max), y_lim=y_lim,
  custom_breaks=custom_breaks_s,
  label_format=label_format,
  base_font=latex_font)

compass_husimi_heatmap <- plot_wigner_heatmap(
  Q$heatmap_dt, empty_overlay,
  q_lim=c(q_lo, q_hi), p_lim=c(p_lo, p_hi),
  custom_breaks_q=custom_breaks_q,
  custom_breaks_p=custom_breaks_p,
  label_format=label_format, base_font=latex_font)

# ------------------------------------------------------------------------------
# SAVE STANDALONE PDF (skipped when sourced by assemble_compass.R)
# Three 2.25"-wide panels side by side at 2.4" tall.
# ------------------------------------------------------------------------------

if (!exists("ASSEMBLY_MODE")) {
  p_final <- compass_wigner_husimi_overlay_heatmap |
             compass_husimi_cross |
             compass_husimi_heatmap
  ggsave(filename=file_output_pdf, plot=p_final, device=cairo_pdf,
         width=6.75, height=2.4, limitsize=FALSE)
  cat("Done.", file_output_pdf, "\n")
}
