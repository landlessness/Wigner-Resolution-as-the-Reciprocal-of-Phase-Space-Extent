# ==============================================================================
# plot_compass_wigner_symplectic.R
# Compass-state symplectic resolution panel — bottom row of paper Fig. 3.
#
# When run standalone, produces a 3-panel PDF showing the full convolution
# story for the symplectic kernel pair:
#
#   compass_wigner_symplectic_overlay_heatmap | compass_symplectic_cross | compass_symplectic_heatmap
#
#   Left:    Pristine Wigner W(q, p) with the symplectic kernel triple
#            (outer Fermi blob A, inner conjugate quantum blobs a_q and
#            a_p) drawn on top, showing the kinematic structure against
#            the unsmoothed input.
#   Center:  Cross-section P_joint along the SW-NE diagonal, post-
#            convolution.
#   Right:   Joint convolved density
#              P_joint(q, p) = (P_{delta q}(q, p) + P_{delta p}(q, p)) / 2,
#            with the same symplectic overlay. The two conjugate kernels
#            preserve the cardinal-aligned sub-Planck structure (the
#            central chessboard) while smoothing the diagonal-aligned
#            interferences.
#
# When sourced by assemble_compass.R (which sets ASSEMBLY_MODE), all three
# panels remain available as named variables but no PDF is written.
#
# Reference: de Gosson, "Symplectic Geometry and Quantum Mechanics"
#            (Birkhauser 2006); Hudson 1974 Rep. Math. Phys. 6, 249.
# ==============================================================================

library(here)
library(patchwork)

source(here("R", "plot_tools.R"))
source(here("R", "compass_system.R"))           # compass_psi, COMPASS_*
source(here("R", "wigner_density.R"))           # build_wigner_state,
                                                # apply_kernel_density
source(here("R", "symplectic_kernel.R"))        # G_delta_q_kernel_matrix,
                                                # G_delta_p_kernel_matrix,
                                                # symplectic_overlay_layers
source(here("R", "quantum_tools.R"))            # numerical_covariance

latex_font      <- "CMU Serif"
dir_figures     <- here("figures")
if (!dir.exists(dir_figures)) dir.create(dir_figures, recursive=TRUE)
file_output_pdf <- file.path(dir_figures, "compass_wigner_symplectic.pdf")

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

# Compass covariance from Schrödinger pipeline (the figure is in the
# Schrödinger row, so we determine Delta_q, Delta_p from psi rather than
# from a classical orbit).
rs <- numerical_covariance(psi_vec, q_psi, hbar=1.0)
cat(sprintf("  L=%.1f  xi=%.1f\n", COMPASS_L, COMPASS_XI))
cat(sprintf("  A_RS/A0=%.2f | Delta_q=%.3f Delta_p=%.3f | delta_q=%.3f delta_p=%.3f\n",
            rs$A_over_A0, rs$Delta_q, rs$Delta_p, rs$delta_q, rs$delta_p))

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
# APPLY THE TWO CONJUGATE KERNELS
# ------------------------------------------------------------------------------

cat("\nApplying G_delta_q kernel...\n")
gdq_for_state <- function(qg, pg) G_delta_q_kernel_matrix(qg, pg,
                                                          rs$Delta_q, rs$Delta_p)
P_dq <- apply_kernel_density(state, gdq_for_state,
                             q_lo, q_hi, p_lo, p_hi, q_display)

cat("\nApplying G_delta_p kernel...\n")
gdp_for_state <- function(qg, pg) G_delta_p_kernel_matrix(qg, pg,
                                                          rs$Delta_q, rs$Delta_p)
P_dp <- apply_kernel_density(state, gdp_for_state,
                             q_lo, q_hi, p_lo, p_hi, q_display)

# ------------------------------------------------------------------------------
# JOINT DENSITY
# Average the two conjugate-kernel densities.
# ------------------------------------------------------------------------------

cat("\nAveraging to P_joint = (P_delta_q + P_delta_p) / 2 ...\n")
P_joint_matrix <- 0.5 * (P_dq$P_matrix + P_dp$P_matrix)

# Build a heatmap_dt for the joint density. Each input was independently
# normalized for its colormap, so we rebuild the joint heatmap from the
# averaged matrix.
q_mask <- state$q_int >= q_lo & state$q_int <= q_hi
p_mask <- state$p_int >= p_lo & state$p_int <= p_hi
joint_heatmap_dt <- as.data.table(expand.grid(
  q = state$q_int[q_mask],
  p = state$p_int[p_mask]
))
joint_heatmap_dt[, w := as.vector(P_joint_matrix[q_mask, p_mask])]
max_w_joint <- max(joint_heatmap_dt$w, na.rm=TRUE)
joint_heatmap_dt[, w_plot := if (max_w_joint > 0) w/max_w_joint else w]

cat(sprintf("  P_joint min over heatmap region: %.4e (Hudson: should be >= 0)\n",
            min(joint_heatmap_dt$w)))

# ------------------------------------------------------------------------------
# DIAGONAL CROSS-SECTION (post-convolution)
# ------------------------------------------------------------------------------

cat("\nExtracting SW-NE diagonal cross-section...\n")
P_diag <- extract_diagonal_cross_section(P_joint_matrix,
                                         state$q_int, state$p_int,
                                         s_display)
dt_cross <- data.table(q=s_display, P_sympl=P_diag)
cat(sprintf("  P_joint diagonal cross-section range: [%.4e, %.4e]\n",
            min(P_diag), max(P_diag)))

p_peak <- max(dt_cross$P_sympl, na.rm=TRUE)
y_lim  <- p_peak * 1.1

# ------------------------------------------------------------------------------
# OVERLAY: three solid black ellipses (outer A + inner a_q, a_p).
# Drawn on the column-1 panel (kernel sitting on the input Wigner) but
# NOT on the column-3 panel (the overlay would be visually redundant on
# the already-convolved output).
# ------------------------------------------------------------------------------

overlay_layers <- symplectic_overlay_layers(rs$Delta_q, rs$Delta_p, q_center=0)
empty_overlay  <- list()

# ------------------------------------------------------------------------------
# BUILD NAMED PANELS
# Three panels for the standalone 3-column layout; assemble_compass.R
# uses all three.
# ------------------------------------------------------------------------------

compass_wigner_symplectic_overlay_heatmap <- plot_wigner_heatmap(
  state$heatmap_dt, overlay_layers,
  q_lim=c(q_lo, q_hi), p_lim=c(p_lo, p_hi),
  custom_breaks_q=custom_breaks_q,
  custom_breaks_p=custom_breaks_p,
  label_format=label_format, base_font=latex_font)

compass_symplectic_cross <- plot_symplectic_cross_section(
  dt_cross,
  q_lim=c(-s_max, s_max), y_lim=y_lim,
  custom_breaks=custom_breaks_s,
  label_format=label_format,
  base_font=latex_font)

compass_symplectic_heatmap <- plot_wigner_heatmap(
  joint_heatmap_dt, empty_overlay,
  q_lim=c(q_lo, q_hi), p_lim=c(p_lo, p_hi),
  custom_breaks_q=custom_breaks_q,
  custom_breaks_p=custom_breaks_p,
  label_format=label_format, base_font=latex_font)

# ------------------------------------------------------------------------------
# SAVE STANDALONE PDF (skipped when sourced by assemble_compass.R)
# Three 2.25"-wide panels side by side at 2.4" tall.
# ------------------------------------------------------------------------------

if (!exists("ASSEMBLY_MODE")) {
  p_final <- compass_wigner_symplectic_overlay_heatmap |
             compass_symplectic_cross |
             compass_symplectic_heatmap
  ggsave(filename=file_output_pdf, plot=p_final, device=cairo_pdf,
         width=6.75, height=2.4, limitsize=FALSE)
  cat("Done.", file_output_pdf, "\n")
}
