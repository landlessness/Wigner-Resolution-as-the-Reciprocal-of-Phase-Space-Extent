# ==============================================================================
# plot_double_well_wigner_symplectic.R
# Symplectic resolution of Wigner negativity — symmetric double well.
#
# Three rows indexed by quantum number:
#   n=0  lowest sub-barrier tunneling-doublet partner (symmetric)
#   n=2  upper sub-barrier doublet partner (still below the barrier)
#   n=7  above-barrier connected orbit
#
# Three columns:
#   Left:    Wigner heatmap W_n(q,p) with QoA overlay (squeezed kernel cells)
#            and the classical orbit at energy E_n drawn on top in gray.
#            Sub-barrier states have two disconnected orbital loops; the
#            above-barrier state has a single peanut-shaped orbit. The
#            classical_trajectory() helper handles both via its NA-row-
#            separator convention.
#   Center:  W_n(q, 0) cross-section. For doublet partners the central
#            interference peak is the phase-space signature of coherent
#            superposition between the two wells.
#   Right:   P_{delta q}(q, 0) cross-section, the W * G_{delta q} convolution
#            evaluated at p = 0. Non-negative everywhere by Hudson's theorem.
#            The Husimi cross-section Q(q, 0) is overlaid with semi-
#            transparent fill so the kernel-capacity contrast is visible
#            in a single panel: where the symplectic curve resolves the
#            three-peak structure (left well, central interference, right
#            well), the Husimi curve flattens it.
#
# The pipeline is identical to the Morse Wigner figure: same kernel-
# agnostic apply_kernel_cross_section() helper, same QoA overlay, same
# orbit_covariance() to derive the kernel widths from the classical
# orbit (which for sub-barrier states integrates over both disconnected
# loops to produce a single set of widths spanning both wells). Only
# the system module differs.
# ==============================================================================

library(here)
library(patchwork)

source(here("R", "plot_tools.R"))
source(here("R", "double_well_system.R"))
source(here("R", "classical_action_tools.R"))   # orbit_covariance,
                                                # classical_trajectory
source(here("R", "wigner_density.R"))           # build_wigner_state,
                                                # apply_kernel_cross_section
source(here("R", "symplectic_kernel.R"))        # G_delta_q_kernel_matrix,
                                                # symplectic_overlay_layers
source(here("R", "husimi_kernel.R"))            # husimi_kernel_matrix
source(here("R", "schroedinger_solver.R"))      # solve_schroedinger

latex_font      <- "CMU Serif"
dir_figures     <- here("figures")
if (!dir.exists(dir_figures)) dir.create(dir_figures, recursive=TRUE)
file_output_pdf <- file.path(dir_figures, "double_well_wigner_symplectic.pdf")

# Quantum numbers for the three rows.
#   n=0, 2: sub-barrier doublet partners. n=0 is the symmetric ground state;
#           n=2 is the upper symmetric (with n=1 being the antisymmetric
#           partner of n=0, n=3 the antisymmetric partner of n=2). We pick
#           the symmetric partners — they have the central peak feature
#           the symplectic kernel is meant to preserve, while the
#           antisymmetric partners would have a node there.
#   n=7:    connected above-barrier orbit. The double-well barrier is at
#           V_0 = 4; for our parameters the spectrum has E_7 well above
#           the barrier so the classical orbit is a single connected loop.
target_quantum_numbers <- c(0, 2, 7)

# Solve Schroedinger once for all rows. Wigner FFT consumes psi_n.
cat("Solving Schroedinger for double well...\n")
dw_soln <- solve_schroedinger(double_well_V,
                              DOUBLE_WELL_Q_MIN, DOUBLE_WELL_Q_MAX,
                              DOUBLE_WELL_DQ,
                              n_states=DOUBLE_WELL_N_STATES)

# ------------------------------------------------------------------------------

build_double_well_row <- function(n, base_font="") {

  # Energy and turning points. Take energy from the Schroedinger solver
  # (matches the wavefunction the heatmap renders); turning points from
  # the polynomial-root solver in double_well_system.R, which returns
  # 4 real roots below the barrier and 2 above.
  E_n <- dw_soln$energies[n + 1]
  tp_obj <- double_well_turning_points(E_n)
  tp     <- tp_obj$roots
  is_sub_barrier <- length(tp) == 4

  if (is_sub_barrier) {
    cat(sprintf("\n== n=%d | E_n=%.4f | sub-barrier | tps=[%.3f, %.3f, %.3f, %.3f] ==\n",
                n, E_n, tp[1], tp[2], tp[3], tp[4]))
  } else {
    cat(sprintf("\n== n=%d | E_n=%.4f | above-barrier | tps=[%.3f, %.3f] ==\n",
                n, E_n, tp[1], tp[2]))
  }

  # Orbit covariance. The helper handles both 2-tp (single orbit) and
  # 4-tp (two disconnected orbits) cases — for sub-barrier states the
  # moments are accumulated over both segments, producing a single
  # (Delta_q, Delta_p) that envelops both wells.
  cov <- orbit_covariance(double_well_V, E_n, tp)
  cat(sprintf("  A_orbit/A0=%.4f | <q>=%.4f | Delta_q=%.3f Delta_p=%.3f\n",
              cov$A_over_A0, cov$q_mean, cov$Delta_q, cov$Delta_p))
  cat(sprintf("  delta_q=%.3f delta_p=%.3f\n",
              cov$delta_q, cov$delta_p))

  # Display window. For sub-barrier doublet states we need a window that
  # spans both wells comfortably; for the above-barrier state we use the
  # outer turning points. Pad with a fraction of the orbit width.
  q_outer_lo <- min(tp)
  q_outer_hi <- max(tp)
  q_span     <- q_outer_hi - q_outer_lo
  q_pad      <- q_span * 0.2
  q_lo       <- q_outer_lo - q_pad
  q_hi       <- q_outer_hi + q_pad

  # Momentum window: kinematic max p set by E_n - V_min (V_min = -V_barrier).
  V_min  <- -double_well_barrier
  p_max  <- sqrt(2 * (E_n - V_min))
  p_lo   <- -1.3 * p_max
  p_hi   <-  1.3 * p_max

  custom_breaks_q <- round(c(q_outer_lo, q_outer_hi), 1)
  custom_breaks_p <- round(c(-p_max, 0, p_max), 1)
  label_format    <- function(x) sprintf("%.1f", x)
  q_display       <- seq(q_lo, q_hi, length.out=500)

  # Build the Wigner state from the n-th Schroedinger eigenstate.
  psi_vec <- dw_soln$psi_matrix[, n + 1]
  state   <- build_wigner_state(psi_vec, dw_soln$q_grid,
                                q_lo, q_hi, p_lo, p_hi, q_display)

  # Kernel cross-sections at p=0. Same kernel-agnostic helper for both;
  # only the kernel matrix differs.
  symplectic_kernel_for_state <- function(q_grid, p_grid) {
    G_delta_q_kernel_matrix(q_grid, p_grid, cov$Delta_q, cov$Delta_p)
  }
  husimi_kernel_for_state <- function(q_grid, p_grid) {
    husimi_kernel_matrix(q_grid, p_grid)
  }
  P_cross <- apply_kernel_cross_section(state, symplectic_kernel_for_state,
                                        q_display)
  Q_cross <- apply_kernel_cross_section(state, husimi_kernel_for_state,
                                        q_display)

  # ----- Build the three panels -----

  # LEFT: Wigner heatmap with QoA cells (centered on q_mean) and the
  # classical orbit drawn over the heatmap. classical_trajectory()
  # handles the multi-segment sub-barrier case via NA-row separators.
  overlay_layers <- symplectic_overlay_layers(cov$Delta_q, cov$Delta_p,
                                              q_center=cov$q_mean)
  df_traj <- classical_trajectory(double_well_V, E_n, tp)

  p_left <- plot_wigner_heatmap(
    state$heatmap_dt, overlay_layers,
    df_traj=df_traj,
    q_lim=c(q_lo,q_hi), p_lim=c(p_lo,p_hi),
    custom_breaks_q=custom_breaks_q,
    custom_breaks_p=custom_breaks_p,
    label_format=label_format, base_font=base_font)

  # CENTER: W_n(q, 0) cross-section. Signed.
  W_peak <- max(abs(state$W_cross), na.rm=TRUE)
  if (!is.finite(W_peak) || W_peak == 0) W_peak <- 1
  y_lim_W <- W_peak * 1.1
  dt_W <- data.table(q=q_display, W_raw=state$W_cross)

  p_center <- plot_wigner_cross_section(
    dt_W, q_lim=c(q_lo,q_hi), y_lim=y_lim_W,
    custom_breaks=custom_breaks_q,
    label_format=label_format, base_font=base_font)

  # RIGHT: P_{delta q}(q, 0) cross-section, non-negative.
  # Husimi cross-section overlaid with semi-transparent fill — same
  # convention as the Morse Wigner figure.
  P_peak <- max(c(P_cross, Q_cross), na.rm=TRUE)
  if (!is.finite(P_peak) || P_peak == 0) P_peak <- 1
  y_lim_P <- P_peak * 1.1
  dt_P <- data.table(q=q_display, rho_sympl=P_cross)

  husimi_overlay <- list(
    data       = data.frame(q = q_display, rho = Q_cross),
    color      = "gray70",
    linewidth  = 0.2,
    fill       = "gray85",
    fill_alpha = 0.5
  )
  overlays <- list(husimi_overlay)

  p_right <- plot_semiclassical_resolution(
    dt_P, q_lim=c(q_lo,q_hi), y_lim=y_lim_P,
    custom_breaks=custom_breaks_q,
    label_format=label_format, base_font=base_font,
    overlays=overlays,
    y_label=expression(italic(P)[italic(delta*q)](italic(q)*","*0)))

  row_label <- sprintf("italic(n)==%d", n)

  list(row_label, p_left, p_center, p_right)
}

cat("Computing double-well Wigner-symplectic grid...\n")
rows    <- lapply(target_quantum_numbers,
                  function(n) build_double_well_row(n, base_font=latex_font))
p_final <- assemble_grid(rows,
                         COLUMN_TITLE_CENTER_WIGNER,
                         COLUMN_TITLE_RIGHT_SYMPLECTIC,
                         base_font=latex_font)

save_figure(p_final, file_output_pdf, length(target_quantum_numbers))
