# ==============================================================================
# plot_morse_wigner_symplectic.R
# Main paper Fig. 3: Symplectic resolution of Wigner negativity — Morse
#
# Three rows indexed by quantum number n = 0, 8, 16 (same selection as the
# semiclassical-symplectic Morse figure). Each row is a single eigenstate of
# the Morse oscillator at its Schroedinger energy.
#
# Three columns:
#   Left:    Wigner heatmap W_n(q,p) with QoA overlay (squeezed kernel cells)
#            and the classical orbit at energy E_n drawn on top in gray.
#            Signed (diverging) colormap.
#   Center:  W_n(q, 0) cross-section. Signed; oscillates negative for n >= 1.
#            The diagnostic for "Wigner negativity to be resolved."
#   Right:   P_{delta q}(q, 0) cross-section, the W * G_{delta q} convolution
#            evaluated at p = 0. Non-negative everywhere by Hudson's theorem.
#            The Husimi cross-section Q(q, 0) is overlaid in gray dotted.
#            Both cross-sections are computed via the same kernel-agnostic
#            apply_kernel_cross_section() helper; only the kernel differs.
#            The contrast makes the kernel-capacity story visible directly
#            on the resolved panel.
#
# The figure pairs naturally with the Wigner cross-section in the middle
# column: cross-section in, cross-section out. The resolution-vs-|psi|^2
# story (which would compare a marginal to a marginal) belongs in the
# semiclassical-symplectic figure, not here.
#
# Pipeline: schroedinger eigenstate -> Wigner via FFT (build_wigner_state) ->
# orbit-derived covariance from V(q), E_n (NOT from psi covariance — see
# manuscript, "Action Capacity and Orbit Action") -> symplectic kernel widths
# delta_q = hbar/Delta_p, delta_p = hbar/Delta_q -> apply_kernel_cross_section
# at p = 0 for the convolved density.
#
# Note on the kernel widths: although the Wigner figure operates on the
# quantum state, the kernel widths come from the same orbit_covariance() the
# semiclassical figure uses, not from the Schroedinger psi covariance. This
# is the manuscript's central design choice: the kernel adapts to the orbit's
# kinematic envelope, regardless of whether the input being convolved is the
# quantum Wigner function or the classical microcanonical shell. The two
# figures share the same kernel and the same QoA overlay; only the input
# differs.
# ==============================================================================

library(here)
library(patchwork)

source(here("R", "plot_tools.R"))
source(here("R", "morse_system.R"))
source(here("R", "classical_action_tools.R"))   # orbit_covariance,
                                                # classical_trajectory
source(here("R", "wigner_density.R"))           # build_wigner_state,
                                                # apply_kernel_cross_section
source(here("R", "symplectic_kernel.R"))        # G_delta_q_kernel_matrix,
                                                # symplectic_overlay_layers,
                                                # symplectic_marginal_density
source(here("R", "husimi_kernel.R"))            # husimi_marginal_density
source(here("R", "schroedinger_solver.R"))      # solve_schroedinger

latex_font      <- "CMU Serif"
dir_figures     <- here("figures")
if (!dir.exists(dir_figures)) dir.create(dir_figures, recursive=TRUE)
file_output_pdf <- file.path(dir_figures, "morse_wigner_symplectic.pdf")

# Quantum numbers for the three rows. Same as the Morse semiclassical figure.
target_quantum_numbers <- c(0, 8, 16)

# Solve Schroedinger once for all rows. The Wigner FFT consumes psi_n.
cat("Solving Schroedinger for Morse...\n")
morse_soln <- solve_schroedinger(morse_V,
                                 MORSE_Q_MIN, MORSE_Q_MAX, MORSE_DQ,
                                 n_states=MORSE_N_STATES)

# ------------------------------------------------------------------------------

build_morse_row <- function(n, base_font="") {

  # Energy and turning points: use the Schroedinger eigenvalue (which for
  # Morse equals the Bohr-Sommerfeld value to numerical precision, but we
  # take it from the same solver that produced psi to keep the Wigner
  # heatmap and the orbit visually exactly co-energetic).
  E_n <- morse_soln$energies[n + 1]
  tp  <- morse_turning_points(E_n)
  q_minus <- tp$q_minus
  q_plus  <- tp$q_plus

  cat(sprintf("\n== n=%d | E_n=%.4f | q-=%.4f | q+=%.4f ==\n",
              n, E_n, q_minus, q_plus))

  # Orbit covariance sets the symplectic kernel widths. Same call as the
  # semiclassical figure — the kernel is orbit-derived, not psi-derived.
  cov <- orbit_covariance(morse_V, E_n, tp)
  cat(sprintf("  A_orbit/A0=%.4f | <q>=%.4f | Delta_q=%.3f Delta_p=%.3f\n",
              cov$A_over_A0, cov$q_mean, cov$Delta_q, cov$Delta_p))
  cat(sprintf("  delta_q=%.3f delta_p=%.3f\n",
              cov$delta_q, cov$delta_p))

  # Display window — same convention as the Morse semiclassical figure.
  q_center <- (q_plus + q_minus) / 2
  q_span   <- q_plus - q_minus
  q_pad    <- q_span * 0.3
  q_lo     <- min(q_minus - q_pad, q_center - 1.3)
  q_hi     <- max(q_plus  + q_pad, q_center + 1.3)
  p_max    <- sqrt(2*E_n)
  p_lo     <- -(p_max * 1.3)
  p_hi     <-   p_max * 1.3

  custom_breaks_q <- round(c(q_minus, q_plus), 1)
  custom_breaks_p <- round(c(-p_max, 0, p_max), 1)
  label_format    <- function(x) sprintf("%.1f", x)
  q_display       <- seq(q_lo, q_hi, length.out=500)

  # Build the Wigner state. Uses the n-th Schroedinger eigenstate as input.
  psi_vec <- morse_soln$psi_matrix[, n + 1]
  state   <- build_wigner_state(psi_vec, morse_soln$q_grid,
                                q_lo, q_hi, p_lo, p_hi, q_display)

  # Apply both kernels and extract their p=0 cross-sections on q_display.
  # Both cross-sections are computed via the same kernel-agnostic helper,
  # apply_kernel_cross_section(); only the kernel closure differs. The
  # right column then plots both side-by-side, making the kernel-capacity
  # contrast visible at a glance.
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

  # LEFT: Wigner heatmap with QoA overlay and classical-orbit overlay.
  # The diverging colormap renders W_n with sign; the orbit (gray) sits
  # over the heatmap as a geometric reference.
  overlay_layers <- symplectic_overlay_layers(cov$Delta_q, cov$Delta_p,
                                              q_center=q_center)
  df_traj <- classical_trajectory(morse_V, E_n, tp)

  p_left <- plot_wigner_heatmap(
    state$heatmap_dt, overlay_layers,
    df_traj=df_traj,
    q_lim=c(q_lo,q_hi), p_lim=c(p_lo,p_hi),
    custom_breaks_q=custom_breaks_q,
    custom_breaks_p=custom_breaks_p,
    label_format=label_format, base_font=base_font)

  # CENTER: W_n(q, 0) cross-section. Signed.
  # y_lim convention required by plot_wigner_cross_section: y_lim = peak * 1.1.
  W_peak <- max(abs(state$W_cross), na.rm=TRUE)
  if (!is.finite(W_peak) || W_peak == 0) W_peak <- 1
  y_lim_W <- W_peak * 1.1
  dt_W <- data.table(q=q_display, W_raw=state$W_cross)

  p_center <- plot_wigner_cross_section(
    dt_W, q_lim=c(q_lo,q_hi), y_lim=y_lim_W,
    custom_breaks=custom_breaks_q,
    label_format=label_format, base_font=base_font)

  # RIGHT: P_{delta q}(q, 0) cross-section, non-negative.
  # Husimi cross-section overlaid with semi-transparent fill so the
  # kernel-capacity contrast is visible: where Husimi extends beyond the
  # symplectic peaks, a lighter-gray fill region is visible above the
  # symplectic ribbon. The dashed line marks the upper boundary of the
  # Husimi fill.
  #
  # We do NOT overlay |psi_n|^2 here. |psi_n|^2 is the q-marginal of W,
  # not a cross-section, so it would be a different kind of object on the
  # same axes. The resolution-vs-|psi|^2 story belongs in the semiclassical
  # figure, whose right column is naturally a marginal.
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

  # Row label: quantum number, matching the semiclassical figure's convention.
  row_label <- sprintf("italic(n)==%d", n)

  list(row_label, p_left, p_center, p_right)
}

cat("Computing Morse Wigner-symplectic grid...\n")
rows    <- lapply(target_quantum_numbers,
                  function(n) build_morse_row(n, base_font=latex_font))
p_final <- assemble_grid(rows,
                         COLUMN_TITLE_CENTER_WIGNER,
                         COLUMN_TITLE_RIGHT_SYMPLECTIC,
                         base_font=latex_font)

save_figure(p_final, file_output_pdf, length(target_quantum_numbers))
