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
#   Left:    classical orbit at energy E_n with QoA overlay (squeezed kernel
#            cells). The orbit is a 1D curve (not a heatmap); the QoA cells
#            are the visual representation of delta_q, delta_p, Delta_q,
#            Delta_p — the resolution scales determining the convolution.
#   Center:  oscillating WKB density |psi_WKB|^2 = (A^2/p) cos^2(S/hbar - pi/4),
#            with n+1 interference lobes between the turning points and
#            divergences at the turning points where p -> 0 (rendered as
#            infinity arrows).
#   Right:   1D position density rho_{delta q}(q), the marginal of
#            (W_cl * G_{delta q}). Airy uniform density rho_Airy(q) overlaid
#            in dashed gray as the prior-art semiclassical comparator.
#
# The pipeline is Schroedinger-free: quantum number -> Bohr-Sommerfeld
# Morse energy E_n -> turning points -> orbit moments give Delta_q, Delta_p,
# hence the symplectic kernel widths. No wavefunction needed for either
# the symplectic curve or the Airy overlay (Airy is computed directly from
# V(q), E_n, and the turning points via the Langer/Miller construction;
# see airy_uniform.R).
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
source(here("R", "airy_uniform.R"))             # airy_uniform_density
                                                # (prior-art semiclassical
                                                # comparator for the
                                                # right-column overlay)

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

  # Build the semiclassical state: oscillating WKB phase-space lift +
  # oscillating WKB density. The lift is the input the symplectic kernel
  # operates on; |psi_WKB|^2 feeds the middle column. The lift itself is
  # not drawn (the left column shows the classical orbit with the QoA
  # overlay representing the convolution that turns the lift into the
  # right column's resolved density).
  state <- build_semiclassical_state(E_n, morse_V,
                                     q_lo, q_hi, p_lo, p_hi, q_display,
                                     q_minus=q_minus, q_plus=q_plus)

  # Classical orbit trajectory for the left column.
  df_traj <- classical_trajectory(morse_V, E_n, tp)

  # Apply symplectic kernel and marginalize. Closure binds Delta_q, Delta_p
  # for this state.
  symplectic_kernel_for_state <- function(q_grid, p_grid) {
    G_delta_q_kernel_matrix(q_grid, p_grid, cov$Delta_q, cov$Delta_p)
  }
  rho_sympl <- symplectic_marginal_density(state, symplectic_kernel_for_state,
                                           q_display)

  # Airy uniform density on the display grid. Computed here (rather than
  # just before the overlay is built) so its peak can contribute to the
  # right-column y-axis scaling — without this the overlay can shoot past
  # the panel top in cases where rho_Airy is more peaked than the
  # symplectic-resolved density (e.g. ground state). The Airy density is
  # the standard prior-art semiclassical comparator: Langer's (1937)
  # uniform wavefunction extended to bound states by Miller (1968) via
  # midpoint patching of two single-turning-point Langer functions. See
  # airy_uniform.R for the construction details and citation chain.
  rho_airy <- airy_uniform_density(q_display, E_n, morse_V, tp)

  # Y-scaling for middle column. The right-column scaling is handled
  # internally by plot_semiclassical_resolution_split (independent y-axes
  # for the two sub-panels).
  #
  # Middle column: y_lim is driven by the peak of the oscillating
  # |psi_WKB|^2 within the orbit (excluding the divergent cells right
  # at the turning points, which are rendered separately by the
  # infinity arrows). This is the largest finite amplitude the curve
  # will display; setting y_lim above it gives the oscillation lobes
  # room to breathe and leaves headroom for the infinity arrows. We
  # use the maximum finite value within the strictly-allowed region
  # (q_minus + small_pad, q_plus - small_pad) so a near-turning-point
  # spike doesn't dominate the scale.
  inside_pad   <- 0.02 * (q_plus - q_minus)
  inside_mask  <- (q_display > q_minus + inside_pad) &
                  (q_display < q_plus  - inside_pad)
  finite_inside <- state$wkb_density[inside_mask &
                                     is.finite(state$wkb_density) &
                                     state$wkb_density > 0]
  if (length(finite_inside) > 0) {
    osc_peak     <- max(finite_inside)
    y_lim_caustic <- osc_peak / GOLDEN_FILL
  } else {
    # Fallback: smooth envelope at orbit center
    q_center_idx <- which.min(abs(q_display - q_center))
    caustic_floor <- state$wkb_density_smooth[q_center_idx]
    if (!is.finite(caustic_floor) || caustic_floor <= 0) caustic_floor <- 1
    y_lim_caustic <- caustic_floor / (1 - GOLDEN_FILL)
  }

  dt_caustic <- data.table(q=q_display, wkb_density=state$wkb_density)
  dt_rho     <- data.table(q=q_display, rho_sympl=rho_sympl)

  # Y-scaling for the right column: include both symplectic and the
  # Airy overlay in the peak so neither curve clips at the panel top.
  rho_peak  <- max(c(rho_sympl, rho_airy), na.rm=TRUE)
  y_lim_rho <- rho_peak / (GOLDEN_FILL + 0.2)

  # QoA overlay (centered on orbit, same as Wigner figure).
  overlay_layers <- symplectic_overlay_layers(cov$Delta_q, cov$Delta_p,
                                              q_center=q_center)

  # Airy overlay for the right column. The Airy uniform density is the
  # textbook prior-art semiclassical resolution of the WKB caustic — the
  # established alternative against which the symplectic resolution is
  # contrasted. rho_airy was already computed above for y-scaling.
  #
  # This contrasts with the Wigner figure, where the right-column overlay
  # is the Husimi cross-section: there the alternative phase-space
  # resolution method is Husimi; here the alternative semiclassical
  # resolution method is Airy. The architectural parallel is intentional —
  # each figure contrasts the symplectic kernel against the established
  # alternative *of the same kind*.
  airy_overlay <- list(
    list(
      data       = data.frame(q = q_display, rho = rho_airy),
      color      = "gray70",
      linewidth  = 0.3 #,
      # fill       = "gray85",
      # fill_alpha = 0.5
    )
  )

  # Row label: quantum number. We index rows by n (not A_orbit/A_0) for
  # consistency across all bound-state figures in the manuscript — Wigner,
  # symplectic, Husimi, and Airy versions of the Morse and double-well
  # figures should all use the same row index so a reader comparing two
  # figures can match rows by n at a glance. The achieved A_orbit/A_0
  # values are reported to the console (above) and belong in the figure
  # caption rather than the row label itself.
  row_label <- sprintf("italic(n)==%d", n)

  list(
    row_label,
    plot_classical_orbit_phase_space(
      df_traj,
      q_lim=c(q_lo,q_hi), p_lim=c(p_lo,p_hi),
      custom_breaks_q=custom_breaks_q,
      custom_breaks_p=custom_breaks_p,
      label_format=label_format, base_font=base_font,
      overlay_layers=overlay_layers,
      orbit_color=HEATMAP_COLOR_HIGH,
      orbit_linewidth=0.3),
    plot_wkb_caustic_cross_section(
      dt_caustic, q_lim=c(q_lo,q_hi), y_lim=y_lim_caustic,
      custom_breaks=custom_breaks_q,
      label_format=label_format,
      q_minus=q_minus, q_plus=q_plus,
      base_font=base_font),
    plot_semiclassical_resolution(
      dt_rho, q_lim=c(q_lo,q_hi), y_lim=y_lim_rho,
      custom_breaks=custom_breaks_q,
      label_format=label_format, base_font=base_font,
      overlays=airy_overlay)
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
