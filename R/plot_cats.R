# ==============================================================================
# plot_cats.R
# Symplectic resolution of Wigner negativity for n-cat states.
# Each row shows one n-cat configuration:
#   row 1: n_cats=2     (N/S lobes; classic two-cat with central interference)
#   row 2: n_cats=3     (triangle pointing up; pair-fringe envelopes)
#   row 3: n_cats=4     (rotated-45 Zurek compass; lobes at NE/NW/SW/SE,
#                        central chessboard plus pair-fringe envelopes
#                        at the diagonal pair midpoints)
#
# Three columns:
#   left:    bare Wigner heatmap with the symplectic quantum of action
#            (Fermi blob + conjugate quantum blobs) overlaid. Widths from
#            the Robertson-Schroedinger covariance of the cat wavefunction
#            itself.
#   center:  raw Wigner cross-section W(q, 0) at p=0 with negativity
#            visible as the curve dipping below zero.
#   right:   symplectic resolution P_{delta q}(q, 0) at p=0 (the Wigner
#            convolved with the symplectic kernel G_{delta q}, sliced at
#            p=0). Husimi cross-section Q(q, 0) overlaid as the prior-art
#            quantum comparator.
#
# QUANTUM UNIVERSE -- WHAT THIS PIPELINE TOUCHES (and does NOT touch):
#
#   This pipeline lives entirely in the quantum universe:
#     * cat wavefunction psi(q) = sum_k exp(-(q-q_k)^2/(2 xi^2)) exp(i p_k q)
#     * Wigner function W(q, p) via Fourier transform of psi
#     * Robertson-Schroedinger covariance of psi: <q^2>_psi, <p^2>_psi
#       compute Delta_q^RS, Delta_p^RS via numerical_covariance(psi)
#     * symplectic kernel G_{delta q} convolved with W
#     * Husimi kernel for the prior-art comparator overlay
#
#   This pipeline does NOT touch:
#     * any classical orbit (cat states have none -- they are not
#       eigenstates of any potential)
#     * orbit covariance (no orbit_covariance)
#     * Bohr-Sommerfeld classical action (cat states have no BS condition)
#     * oscillating WKB densities, Airy/Langer constructions
#
#   The n_cats parameter is a state-construction parameter, NOT a quantum
#   number. There is no semiclassical "action capacity = n_cats * A_0"
#   relation; the cat state's action capacity comes entirely from the
#   wavefunction's RS covariance.
#
# Author: Brian S. Mulloy
# ==============================================================================

library(here)
library(data.table)
library(patchwork)

source(here("R", "plot_tools.R"))

# Quantum-universe modules ONLY. No classical_action_tools, no
# semiclassical_density, no airy_uniform, no bound-state system modules.
source(here("R", "cat_system.R"))                # cat_psi, cat_lobe_positions,
                                                 # CAT_* parameters
source(here("R", "quantum_tools.R"))             # numerical_covariance
source(here("R", "wigner_density.R"))            # build_wigner_state,
                                                 # apply_kernel_cross_section
source(here("R", "symplectic_kernel.R"))         # G_delta_q_kernel_matrix,
                                                 # symplectic_overlay_layers
source(here("R", "husimi_kernel.R"))             # husimi_kernel_matrix

latex_font  <- "CMU Serif"
dir_figures <- here("figures")
if (!dir.exists(dir_figures)) dir.create(dir_figures, recursive=TRUE)
file_output_pdf <- file.path(dir_figures, "cats.pdf")

# ------------------------------------------------------------------------------
# DISPLAY WINDOW (shared across all three rows)
#
# For diag and 2-cat/3-cat rows, lobes anchored at +/- CAT_P_MAX = 5.
# For the 4-cat axis variant, lobes at +/- p_max*sqrt(2) ~ 7.07 (the 45
# deg rotation of diag, preserving adjacent-cat distance). CAT_Q_DISPLAY
# = 10 covers all four rows with margin around the outermost lobes.
# All rows use the same window and the same axis tick marks at +/-5
# so the rows align visually across the figure; the axis-variant row 4
# will have its lobes just past the +/-5 ticks at +/-7.07.
# ------------------------------------------------------------------------------

q_lo <- -CAT_Q_DISPLAY
q_hi <-  CAT_Q_DISPLAY
p_lo <- -CAT_Q_DISPLAY
p_hi <-  CAT_Q_DISPLAY

custom_breaks_q <- c(-CAT_P_MAX, 0, CAT_P_MAX)
custom_breaks_p <- c(-CAT_P_MAX, 0, CAT_P_MAX)
label_format    <- function(x) sprintf("%.0f", x)
q_display       <- seq(q_lo, q_hi, length.out=500)

# ------------------------------------------------------------------------------
# ROW DESCRIPTORS
#
# Cat states don't fit the bound-state descriptor pattern (no E, no
# orbit, no V); they're parameterized only by n_cats and (for n_cats=4)
# a variant naming the lobe orientation. List the four rows explicitly:
#   row 1: n_cats=2          -- N/S
#   row 2: n_cats=3          -- triangle apex up
#   row 3: n_cats=4 diag     -- compass with lobes on the diagonals at
#                               NE/NW/SW/SE. Lobes are off the p=0
#                               slice; cross-section sees only on-slice
#                               cross-fringe envelopes and the central
#                               sub-Planck chessboard.
#   row 4: n_cats=4 axis     -- diag rotated 45 deg, with the same
#                               adjacent-cat distance. Lobes E and W
#                               now sit on p=0; the cross-section
#                               directly resolves them. Same physics
#                               as row 3 up to phase-space rotation
#                               (same A_RS/A_0, same kernel widths,
#                               same QoA size); what changes is what
#                               the p=0 slice projects.
#
# Pairing diag and axis with identical adjacent-cat distance makes the
# rotation invariance of the construction visible: same kernel, same
# action capacity, different cross-section content.
# ------------------------------------------------------------------------------

row_configs <- list(
  list(n_cats=2, variant="diag"),
  list(n_cats=3, variant="diag"),
  list(n_cats=4, variant="diag"),
  list(n_cats=4, variant="axis")
)

# ------------------------------------------------------------------------------
# CAT-WIGNER ROW BUILDER
#
# Pure quantum-universe pipeline. Takes n_cats, returns the three-panel
# row.
# ------------------------------------------------------------------------------

build_cat_row <- function(config, base_font="") {
  n_cats  <- config$n_cats
  variant <- config$variant
  cat(sprintf("\n== n_cats=%d variant=%s ==\n", n_cats, variant))

  # 1. Sample psi on the cat wavefunction grid.
  q_psi   <- seq(CAT_Q_MIN, CAT_Q_MAX, by=CAT_DQ)
  psi_vec <- cat_psi(q_psi, n_cats, variant=variant,
                     p_max=CAT_P_MAX, xi=CAT_XI, hbar=CAT_HBAR)

  # 2. Robertson-Schroedinger covariance.
  rs <- numerical_covariance(psi_vec, q_psi, hbar=CAT_HBAR)
  cat(sprintf("  RS-covariance: A_RS/A0=%.4f Delta_q=%.3f Delta_p=%.3f\n",
              rs$A_over_A0, rs$Delta_q, rs$Delta_p))
  cat(sprintf("                 delta_q=%.3f delta_p=%.3f q_mean=%.3f\n",
              rs$delta_q, rs$delta_p, rs$q_mean))

  # 3. Wigner state. build_wigner_state computes W via FFT of psi.
  cat("  Building Wigner state...\n")
  state <- build_wigner_state(
    psi_vec=psi_vec, psi_q_grid=q_psi,
    q_lo=q_lo, q_hi=q_hi, p_lo=p_lo, p_hi=p_hi,
    q_display=q_display)

  # 4. Symplectic kernel: shape Gaussian, area h/2; widths from RS.
  symplectic_kernel_for_state <- function(qg, pg) {
    G_delta_q_kernel_matrix(qg, pg, rs$Delta_q, rs$Delta_p, hbar=CAT_HBAR)
  }
  cat("  Convolving with symplectic kernel...\n")
  P_sympl_cross <- apply_kernel_cross_section(state, symplectic_kernel_for_state,
                                              q_display)

  # 5. Husimi: same cross-section pipeline, fixed-width Gaussian kernel.
  husimi_kernel_for_state <- function(qg, pg) {
    husimi_kernel_matrix(qg, pg)
  }
  cat("  Computing Husimi cross-section at p=0...\n")
  Q_husimi_cross <- apply_kernel_cross_section(state, husimi_kernel_for_state,
                                               q_display)

  # 6. Y-scaling.
  W_cross_peak <- max(abs(state$W_cross), na.rm=TRUE)
  if (!is.finite(W_cross_peak) || W_cross_peak == 0) W_cross_peak <- 1
  y_lim_W <- W_cross_peak * 1.1

  P_peak_data <- max(P_sympl_cross, na.rm=TRUE)
  if (!is.finite(P_peak_data) || P_peak_data == 0) P_peak_data <- 1
  y_lim_P <- P_peak_data * 1.1

  # 7. Build cross-section data tables.
  dt_W       <- data.table(q=q_display, W_raw=state$W_cross)
  dt_P_sympl <- data.table(q=q_display, rho_sympl=P_sympl_cross)
  dt_Husimi  <- data.frame(q=q_display, rho=Q_husimi_cross)

  # 8. QoA overlay. RS-derived widths, centered at the wavefunction's
  #    <q>. For symmetric cat configurations <q> = 0; non-symmetric
  #    layouts could shift the center.
  overlay_layers <- symplectic_overlay_layers(rs$Delta_q, rs$Delta_p,
                                              q_center=rs$q_mean,
                                              hbar=CAT_HBAR)

  # 9. Husimi overlay on the symplectic resolution panel.
  husimi_overlay <- list(
    list(
      data       = dt_Husimi,
      color      = "gray30",
      linewidth  = 0.35,
      fill       = "gray85",
      fill_alpha = 0.5
    )
  )

  # 10. Assemble three panels.
  list(
    plot_wigner_heatmap(
      state$heatmap_dt, overlay_layers, df_traj=NULL,
      q_lim=c(q_lo, q_hi), p_lim=c(p_lo, p_hi),
      custom_breaks_q=custom_breaks_q,
      custom_breaks_p=custom_breaks_p,
      label_format=label_format, base_font=base_font),
    plot_wigner_cross_section(
      dt_W, q_lim=c(q_lo, q_hi), y_lim=y_lim_W,
      custom_breaks=custom_breaks_q,
      label_format=label_format, base_font=base_font),
    plot_semiclassical_resolution(
      dt_P_sympl, q_lim=c(q_lo, q_hi), y_lim=y_lim_P,
      custom_breaks=custom_breaks_q,
      label_format=label_format, base_font=base_font,
      overlays=husimi_overlay,
      y_label=expression(italic(P)[italic(delta*q)](italic(q)*","*0)))
  )
}

# ------------------------------------------------------------------------------
# DRIVE: build all three rows, assemble, save.
# ------------------------------------------------------------------------------

cat("Computing cats figure (4 cat configurations x 3 panels)...\n")

rows <- lapply(row_configs,
               function(cfg) build_cat_row(cfg, base_font=latex_font))

p_final <- assemble_grid_unlabeled(rows,
                                   COLUMN_TITLE_CENTER_WIGNER,
                                   COLUMN_TITLE_RIGHT_SYMPLECTIC,
                                   base_font=latex_font)

save_figure(p_final, file_output_pdf, length(row_configs))
