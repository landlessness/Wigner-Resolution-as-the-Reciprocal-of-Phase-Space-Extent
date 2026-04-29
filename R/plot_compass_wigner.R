# ==============================================================================
# plot_compass_wigner.R
# Compass-state Wigner ground-truth panel — top row of paper Fig. 3.
#
# When run standalone, produces a 2-panel PDF:
#   compass_wigner_heatmap | compass_wigner_cross
#
# When sourced by assemble_compass.R (which sets ASSEMBLY_MODE), the
# panels remain available as named variables but no PDF is written.
#
# Reference: Zurek 2001 Nature 412, 712.
# ==============================================================================

library(here)
library(patchwork)

source(here("R", "plot_tools.R"))
source(here("R", "compass_system.R"))           # compass_psi, COMPASS_*
source(here("R", "wigner_density.R"))           # build_wigner_state

latex_font      <- "CMU Serif"
dir_figures     <- here("figures")
if (!dir.exists(dir_figures)) dir.create(dir_figures, recursive=TRUE)
file_output_pdf <- file.path(dir_figures, "compass_wigner.pdf")

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
# DIAGONAL CROSS-SECTION
# ------------------------------------------------------------------------------

cat("\nExtracting SW-NE diagonal cross-section...\n")
W_diag <- extract_diagonal_cross_section(state$W_matrix,
                                         state$q_int, state$p_int,
                                         s_display)
dt_cross <- data.table(q=s_display, W_raw=W_diag)
cat(sprintf("  W diagonal cross-section range: [%.4f, %.4f]\n",
            min(W_diag), max(W_diag)))

w_max <- max(abs(dt_cross$W_raw), na.rm=TRUE)
y_lim <- w_max * 1.1

# ------------------------------------------------------------------------------
# BUILD NAMED PANELS
# Both panels are also referenced by assemble_compass.R.
# ------------------------------------------------------------------------------

empty_overlay <- list()

compass_wigner_heatmap <- plot_wigner_heatmap(
  state$heatmap_dt, empty_overlay,
  q_lim=c(q_lo, q_hi), p_lim=c(p_lo, p_hi),
  custom_breaks_q=custom_breaks_q,
  custom_breaks_p=custom_breaks_p,
  label_format=label_format, base_font=latex_font)

compass_wigner_cross <- plot_wigner_cross_section(
  dt_cross,
  q_lim=c(-s_max, s_max), y_lim=y_lim,
  custom_breaks=custom_breaks_s,
  label_format=label_format,
  base_font=latex_font)

# ------------------------------------------------------------------------------
# SAVE STANDALONE PDF (skipped when sourced by assemble_compass.R)
# ------------------------------------------------------------------------------

if (!exists("ASSEMBLY_MODE")) {
  p_final <- compass_wigner_heatmap | compass_wigner_cross
  ggsave(filename=file_output_pdf, plot=p_final, device=cairo_pdf,
         width=4.5, height=2.4, limitsize=FALSE)
  cat("Done.", file_output_pdf, "\n")
}
