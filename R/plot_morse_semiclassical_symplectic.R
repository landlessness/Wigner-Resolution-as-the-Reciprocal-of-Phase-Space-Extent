# ==============================================================================
# plot_morse_semiclassical_symplectic.R
# Main paper Fig. 2: Symplectic resolution of semiclassical caustics — Morse
#
# Three rows indexed by quantum number n. The harmonic-limit identity
# A/A_0 = 2n+1 gives nominal A/A_0 = 1, 17, 33 for n = 0, 8, 16, but Morse
# anharmonicity makes A_orbit/A_0 deviate appreciably from 2n+1 at high n
# (the orbit covariance is asymmetric, and the A_orbit(E) curve is bounded
# above by a peak near but below dissociation). We therefore use n = 0, 8, 16
# directly — same rows as the Wigner Morse figure — and label each row by
# the actually-achieved A_orbit/A_0 from orbit_covariance().
#
# Three columns:
#   Left:    regularized 2D energy shell heatmap with QoA overlay
#   Center:  analytical WKB caustic 1/sqrt(2*(E-V)) — diverges at turning points
#   Right:   1D position density rho_{delta q}(q), the marginal of
#            (W_cl * G_{delta q})
#
# The pipeline is Schroedinger-free: quantum number -> Bohr-Sommerfeld
# Morse energy E_n -> turning points -> orbit moments give Delta_q, Delta_p,
# hence the symplectic kernel widths. No wavefunction needed.
# ==============================================================================

library(here)
library(patchwork)

source(here("R", "plot_tools.R"))
source(here("R", "morse_system.R"))
source(here("R", "classical_action_tools.R"))   # orbit_covariance, classical_action
source(here("R", "semiclassical_density.R"))    # build_semiclassical_state
source(here("R", "symplectic_kernel.R"))        # G_delta_q_kernel_matrix,
                                                # symplectic_marginal_density,
                                                # symplectic_overlay_layers

GOLDEN_FILL <- 0.6

latex_font      <- "CMU Serif"
dir_figures     <- here("figures")
if (!dir.exists(dir_figures)) dir.create(dir_figures, recursive=TRUE)
file_output_pdf <- file.path(dir_figures, "morse_semiclassical_symplectic.pdf")

# Quantum numbers for the three rows. Same selection as the Wigner Morse
# figure: ground state, mid-anharmonic, strongly anharmonic horseshoe.
target_quantum_numbers <- c(0, 8, 16)

# ------------------------------------------------------------------------------

build_morse_row <- function(n, base_font="") {

  # Bohr-Sommerfeld energy at quantum number n. For Morse the BS spectrum
  # is exact, so E_n is also the analytic Schroedinger eigenvalue.
  E_n <- morse_E_BS(n)
  tp  <- morse_turning_points(E_n)
  q_minus <- tp$q_minus
  q_plus  <- tp$q_plus

  cat(sprintf("\n== n=%d | E_n=%.4f | q-=%.4f | q+=%.4f ==\n",
              n, E_n, q_minus, q_plus))

  # Orbit-derived covariance sizes the symplectic kernel.
  cov <- orbit_covariance(morse_V, E_n, tp)
  cat(sprintf("  A_orbit/A0=%.4f | <q>=%.4f | Delta_q=%.3f Delta_p=%.3f\n",
              cov$A_over_A0, cov$q_mean, cov$Delta_q, cov$Delta_p))
  cat(sprintf("  (harmonic-limit nominal A/A0 = 2n+1 = %d)\n", 2*n+1))

  # Display window — same convention as the Wigner Morse figure.
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

  # Build the semiclassical state: regularized 2D energy shell + 1D WKB.
  state <- build_semiclassical_state(E_n, morse_V,
                                     q_lo, q_hi, p_lo, p_hi, q_display)

  # Apply symplectic kernel and marginalize. Closure binds Delta_q, Delta_p
  # for this state.
  symplectic_kernel_for_state <- function(q_grid, p_grid) {
    G_delta_q_kernel_matrix(q_grid, p_grid, cov$Delta_q, cov$Delta_p)
  }
  rho_sympl <- symplectic_marginal_density(state, symplectic_kernel_for_state,
                                           q_display)

  # Y-scaling, independent per column.
  # Right column: density peak fills (GOLDEN_FILL + 0.2) of vertical.
  rho_peak  <- max(rho_sympl, na.rm=TRUE)
  y_lim_rho <- rho_peak / (GOLDEN_FILL + 0.2)

  # Middle column: caustic bowl floor at orbit center sits at
  # (1 - GOLDEN_FILL) of vertical; infinity arrows handle the divergences
  # at the turning points.
  q_center_idx <- which.min(abs(q_display - q_center))
  caustic_floor <- state$wkb_density[q_center_idx]
  if (!is.finite(caustic_floor) || caustic_floor <= 0) {
    finite_caustic <- state$wkb_density[is.finite(state$wkb_density) &
                                          state$wkb_density > 0]
    caustic_floor <- if (length(finite_caustic) > 0)
      min(finite_caustic, na.rm=TRUE) else 1
  }
  y_lim_caustic <- caustic_floor / (1 - GOLDEN_FILL)

  dt_caustic <- data.table(q=q_display, wkb_density=state$wkb_density)
  dt_rho     <- data.table(q=q_display, rho_sympl=rho_sympl)

  # QoA overlay (centered on orbit, same as Wigner figure).
  overlay_layers <- symplectic_overlay_layers(cov$Delta_q, cov$Delta_p,
                                              q_center=q_center)

  # Row label: actually-achieved A_orbit/A_0 to two decimal places. This
  # tells the reader the kinematic envelope size that the symplectic
  # overlay is expressing — not the (inaccurate-for-Morse) harmonic-limit
  # nominal value.
  row_label <- sprintf("%.2f~italic(A)[0]", cov$A_over_A0)

  list(
    row_label,
    plot_semiclassical_heatmap(
      state$heatmap_dt, overlay_layers,
      q_lim=c(q_lo,q_hi), p_lim=c(p_lo,p_hi),
      custom_breaks_q=custom_breaks_q,
      custom_breaks_p=custom_breaks_p,
      label_format=label_format, base_font=base_font),
    plot_wkb_caustic_cross_section(
      dt_caustic, q_lim=c(q_lo,q_hi), y_lim=y_lim_caustic,
      custom_breaks=custom_breaks_q,
      label_format=label_format,
      q_minus=q_minus, q_plus=q_plus,
      base_font=base_font),
    plot_semiclassical_resolution(
      dt_rho, q_lim=c(q_lo,q_hi), y_lim=y_lim_rho,
      custom_breaks=custom_breaks_q,
      label_format=label_format, base_font=base_font)
  )
}

cat("Computing Morse semiclassical-symplectic grid...\n")
rows    <- lapply(target_quantum_numbers,
                  function(n) build_morse_row(n, base_font=latex_font))
p_final <- assemble_grid(rows,
                         COLUMN_TITLE_CENTER_SEMICLASSICAL,
                         COLUMN_TITLE_RIGHT_SYMPLECTIC,
                         base_font=latex_font)

save_figure(p_final, file_output_pdf, length(target_quantum_numbers))
