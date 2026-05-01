# ==============================================================================
# plot_morse_semiclassical_airy.R
# Companion figure: Berry-Mount / Langer Airy uniform approximation for
# the Morse semiclassical caustics — the textbook approach to which our
# symplectic resolution is compared.
#
# Three rows indexed by quantum number n = 0, 8, 16. Same physical states
# as the symplectic figure, so a reader can lay the two figures side by
# side and compare the two methods row-by-row.
#
# Three columns:
#   Left:    classical orbit trajectory in phase space (no QoA overlay,
#            no heatmap — the Airy approach has no action-capacity
#            geometry, so the left column shows only the deterministic
#            classical motion that both methods build on).
#   Center:  analytical WKB caustic 1/sqrt(2*(E-V)) — diverges at turning
#            points. Identical to the corresponding column in the
#            symplectic figure: this is the divergence both methods are
#            trying to resolve.
#   Right:   1D position density rho_Airy(q) from the uniform Airy
#            replacement (Langer 1937, Berry & Mount 1972). Always
#            finite; matches the WKB envelope away from turning points
#            and replaces the divergence with an Airy-function lobe at
#            each turning point.
#
# Pipeline: quantum number -> Bohr-Sommerfeld Morse energy E_n -> turning
# points -> action integrals from each turning point -> Langer
# uniform density. No wavefunction, no orbit covariance, no symplectic
# kernel. The only inputs are V(q) and E.
# ==============================================================================

library(here)
library(patchwork)
library(data.table)

source(here("R", "plot_tools.R"))
source(here("R", "morse_system.R"))
source(here("R", "classical_action_tools.R"))   # classical_trajectory
source(here("R", "semiclassical_density.R"))    # build_semiclassical_state
                                                # (used only for wkb_density)
source(here("R", "airy_uniform.R"))             # airy_uniform_density
source(here("R", "schroedinger_solver.R"))      # solve_schroedinger,
                                                # schroedinger_density

GOLDEN_FILL <- 0.6

latex_font      <- "CMU Serif"
dir_figures     <- here("figures")
if (!dir.exists(dir_figures)) dir.create(dir_figures, recursive=TRUE)
file_output_pdf <- file.path(dir_figures, "morse_semiclassical_airy.pdf")

# Quantum numbers for the three rows. Same selection as the symplectic
# figure for cross-figure comparison.
target_quantum_numbers <- c(0, 8, 16)

# Solve Schroedinger once for all rows. Used for the dashed |psi_n|^2
# ground-truth overlay on the right column.
cat("Solving Schroedinger for Morse (used for ground-truth overlay)...\n")
morse_soln <- solve_schroedinger(morse_V,
                                 MORSE_Q_MIN, MORSE_Q_MAX, MORSE_DQ,
                                 n_states=MORSE_N_STATES)

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

  # Display window — same convention as the symplectic figure so the
  # left and middle columns line up between the two figures.
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

  # Classical orbit trajectory (left column) and WKB caustic (middle).
  # We reuse build_semiclassical_state purely for its wkb_density field
  # (the oscillating |psi_WKB|^2 displayed in the middle column); the
  # heatmap_dt and shell that come with it are unused in this figure.
  state    <- build_semiclassical_state(E_n, morse_V,
                                        q_lo, q_hi, p_lo, p_hi, q_display,
                                        q_minus=q_minus, q_plus=q_plus)
  df_traj  <- classical_trajectory(morse_V, E_n, tp)

  # Airy uniform density on the display grid. This is the right column.
  rho_airy <- airy_uniform_density(q_display, E_n, morse_V, tp)
  cat(sprintf("  rho_Airy: peak=%.4f, integral=%.4f\n",
              max(rho_airy, na.rm=TRUE),
              sum(rho_airy, na.rm=TRUE) * diff(q_display)[1]))

  # Schroedinger ground-truth density for the overlay (and for y-axis
  # peak inclusion so the overlay never gets clipped).
  psi_sq <- schroedinger_density(morse_soln, n, q_display)

  # Y-scaling.
  rho_peak  <- max(c(rho_airy, psi_sq), na.rm=TRUE)
  y_lim_rho <- rho_peak / (GOLDEN_FILL + 0.2)

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
    q_center_idx <- which.min(abs(q_display - q_center))
    caustic_floor <- state$wkb_density_smooth[q_center_idx]
    if (!is.finite(caustic_floor) || caustic_floor <= 0) caustic_floor <- 1
    y_lim_caustic <- caustic_floor / (1 - GOLDEN_FILL)
  }

  dt_caustic <- data.table(q=q_display, wkb_density=state$wkb_density)
  dt_rho     <- data.table(q=q_display, rho_airy=rho_airy)

  schroedinger_overlay <- list(
    list(
      data      = data.frame(q = q_display, rho = psi_sq),
      color     = "gray30",
      linetype  = 11,
      linewidth = 0.2
    )
  )

  # Row label: quantum number, matching the symplectic figure.
  row_label <- sprintf("italic(n)==%d", n)

  list(
    row_label,
    plot_classical_orbit_phase_space(
      df_traj,
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
    plot_airy_resolution(
      dt_rho, q_lim=c(q_lo,q_hi), y_lim=y_lim_rho,
      custom_breaks=custom_breaks_q,
      label_format=label_format, base_font=base_font,
      overlays=schroedinger_overlay)
  )
}

cat("Computing Morse semiclassical-Airy grid...\n")
rows    <- lapply(target_quantum_numbers,
                  function(n) build_morse_row(n, base_font=latex_font))
p_final <- assemble_grid(rows,
                         COLUMN_TITLE_CENTER_SEMICLASSICAL,
                         COLUMN_TITLE_RIGHT_AIRY,
                         base_font=latex_font,
                         title_left=COLUMN_TITLE_LEFT_ORBIT)

save_figure(p_final, file_output_pdf, length(target_quantum_numbers))
