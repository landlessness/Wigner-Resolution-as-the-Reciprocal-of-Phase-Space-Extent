# ==============================================================================
# plot_morse_semiclassical_symplectic.R
# Main paper: Symplectic resolution of semiclassical caustics — Morse potential
#
# Same three rows as the Wigner Morse figure: n=0, 8, 16.
# Three columns:
#   Left:    regularized 2D energy shell heatmap with QoA overlay
#   Center:  analytical WKB caustic 1/sqrt(2*(E-V)) — diverges at turning points
#   Right:   1D position density rho_{delta q}(q), the marginal of
#            (W_cl * G_{delta q})
#
# Architecture mirrors plot_morse_wigner_symplectic.R:
#   - solve_schrodinger() to get eigenstate energies (same energies the
#     Wigner figure uses, ensuring consistent level selection)
#   - build_semiclassical_state() builds the regularized 2D shell on
#     the same display window as Wigner figure
#   - apply_kernel_cross_section() with G_delta_q gives the 1D position
#     density via marginalization over p
# ==============================================================================

library(here)
library(patchwork)

source(here("R", "plot_tools.R"))
source(here("R", "morse_potential.R"))
source(here("R", "wigner_tools.R"))           # for solve_schrodinger
source(here("R", "semiclassical_state.R"))
source(here("R", "math_tools.R"))             # for fft_convolve_2d
source(here("R", "symplectic_kernel.R"))
source(here("R", "classical_action_tools.R"))

latex_font      <- "CMU Serif"
dir_figures     <- here("figures")
if (!dir.exists(dir_figures)) dir.create(dir_figures, recursive=TRUE)
file_output_pdf <- file.path(dir_figures, "morse_semiclassical_symplectic.pdf")

cat("Solving Morse Schrodinger equation...\n")
morse_soln <- solve_schrodinger(morse_V,
                                q_min=MORSE_Q_MIN,
                                q_max=MORSE_Q_MAX,
                                dq=MORSE_DQ,
                                n_states=MORSE_N_STATES)

target_n_levels <- c(0, 8, 16)

# ------------------------------------------------------------------------------

#' Apply symplectic kernel to a semiclassical state's W_matrix and return
#' the 1D position density rho(q) = integral (W * G) dp on q_display.
apply_symplectic_marginal <- function(state, kernel_fn, q_display) {
  K_mat <- kernel_fn(state$q_int, state$p_int)
  conv <- fft_convolve_2d(state$W_matrix, K_mat, state$dq_int, state$dp_int)
  rho_int <- rowSums(conv$P_mat) * state$dp_int
  rho <- approx(state$q_int, rho_int, xout=q_display, rule=1)$y
  rho[is.na(rho)] <- 0
  rho
}

# ------------------------------------------------------------------------------

build_morse_row <- function(n_val, soln, base_font="") {
  E_n <- soln$energies[n_val+1]
  tp  <- morse_turning_points(E_n)
  q_minus <- tp$q_minus
  q_plus  <- tp$q_plus

  cat(sprintf("\n== n=%d | E_n=%.4f | q-=%.4f | q+=%.4f ==\n",
              n_val, E_n, q_minus, q_plus))

  # Use the wavefunction's RS covariance to size the symplectic kernel.
  # This keeps the kernel sized by A (the quantum cell), consistent with
  # the Wigner figures.
  q_grid  <- soln$q_grid
  psi_vec <- soln$psi_matrix[, n_val+1]
  rs <- numerical_covariance(psi_vec, q_grid)
  cat(sprintf("  A_RS/A0=%.2f | Delta_q=%.3f Delta_p=%.3f\n",
              rs$A_over_A0, rs$Delta_q, rs$Delta_p))

  # Display window — same convention as Wigner Morse figure
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

  # Build the semiclassical state (regularized 2D shell + 1D WKB density)
  state <- build_semiclassical_state(E_n, morse_V,
                                     q_lo, q_hi, p_lo, p_hi, q_display)

  # Symplectic kernel for this state — closure binds Delta_q, Delta_p
  symplectic_kernel_for_state <- function(q_grid, p_grid) {
    G_delta_q_kernel_matrix(q_grid, p_grid, rs$Delta_q, rs$Delta_p)
  }
  rho_sympl <- apply_symplectic_marginal(state, symplectic_kernel_for_state, q_display)

  # Build cross-section data tables for the two right columns
  dt_caustic <- data.table(q=q_display, wkb_density=state$wkb_density)

  # Determine y-limit for caustic column: cap at a multiple of the
  # maximum *finite* density. The infinity arrows handle the divergence.
  finite_caustic <- state$wkb_density[is.finite(state$wkb_density) &
                                        state$wkb_density > 0]
  if (length(finite_caustic) > 0) {
    caustic_floor <- median(finite_caustic, na.rm=TRUE)
    y_lim_caustic <- caustic_floor * 6   # caustic value at midpoint × 6
  } else {
    y_lim_caustic <- 1
  }

  # Match right column scale to caustic column for visual comparability
  rho_max <- max(rho_sympl, na.rm=TRUE)
  y_lim_rho <- max(rho_max * 1.3, caustic_floor * 6)
  dt_rho <- data.table(q=q_display, rho_sympl=rho_sympl)

  # QoA overlay (centered on orbit, just like Wigner figure)
  overlay_layers <- symplectic_overlay_layers(rs$Delta_q, rs$Delta_p,
                                              q_center=q_center)

  list(
    sprintf("italic(n)==%d", n_val),
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
rows    <- lapply(target_n_levels,
                  function(n) build_morse_row(n, morse_soln, base_font=latex_font))
p_final <- assemble_grid(rows,
                         COLUMN_TITLE_CENTER_SEMICLASSICAL,
                         COLUMN_TITLE_RIGHT_SYMPLECTIC,
                         base_font=latex_font)
save_figure(p_final, file_output_pdf, length(target_n_levels))
