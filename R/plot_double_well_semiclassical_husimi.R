# ==============================================================================
# plot_double_well_semiclassical_husimi.R
# Appendix: Husimi resolution of semiclassical caustics — double-well
#
# Three rows: n=0, 2, 7. Three columns:
#   Left:    regularized 2D energy shell with Husimi unit circle overlay
#   Center:  analytical WKB caustic 1/sqrt(2*(E-V))
#   Right:   1D position density rho_Q(q), the marginal of (W_cl * Husimi)
# ==============================================================================

library(here)
library(patchwork)

source(here("R", "plot_tools.R"))
source(here("R", "double_well_potential.R"))
source(here("R", "wigner_tools.R"))
source(here("R", "semiclassical_state.R"))
source(here("R", "math_tools.R"))
source(here("R", "husimi_kernel.R"))
source(here("R", "classical_action_tools.R"))

GOLDEN_FILL <- 0.6

latex_font      <- "CMU Serif"
dir_figures     <- here("figures")
if (!dir.exists(dir_figures)) dir.create(dir_figures, recursive=TRUE)
file_output_pdf <- file.path(dir_figures, "double_well_semiclassical_husimi.pdf")

cat("Solving double-well Schrodinger equation...\n")
dw_soln <- solve_schrodinger(double_well_V,
                             q_min=DOUBLE_WELL_Q_MIN,
                             q_max=DOUBLE_WELL_Q_MAX,
                             dq=DOUBLE_WELL_DQ,
                             n_states=DOUBLE_WELL_N_STATES)

target_n_levels <- c(0, 2, 7)
V_BARRIER_TOP <- double_well_V(0)

# ------------------------------------------------------------------------------

apply_husimi_marginal <- function(state, q_display) {
  K_mat <- husimi_kernel_matrix(state$q_int, state$p_int)
  conv  <- fft_convolve_2d(state$W_matrix, K_mat, state$dq_int, state$dp_int)
  rho_int <- rowSums(conv$P_mat) * state$dp_int
  rho <- approx(state$q_int, rho_int, xout=q_display, rule=1)$y
  rho[is.na(rho)] <- 0
  rho
}

# ------------------------------------------------------------------------------

build_double_well_row <- function(n_val, soln, base_font="") {
  E_n  <- soln$energies[n_val+1]
  tp   <- double_well_turning_points(E_n)
  roots <- tp$roots
  is_sub_barrier <- (E_n < V_BARRIER_TOP) && length(roots) == 4

  if (is_sub_barrier) {
    q_minus <- roots[1]
    q_plus  <- roots[4]
    cat(sprintf("\n== n=%d | E_n=%.4f (sub-barrier) ==\n", n_val, E_n))
  } else {
    q_minus <- roots[1]
    q_plus  <- roots[length(roots)]
    cat(sprintf("\n== n=%d | E_n=%.4f (above-barrier) ==\n", n_val, E_n))
  }

  q_grid  <- soln$q_grid
  psi_vec <- soln$psi_matrix[, n_val+1]
  rs <- numerical_covariance(psi_vec, q_grid)
  cat(sprintf("  A_RS/A0=%.2f | Delta_q=%.3f Delta_p=%.3f\n",
              rs$A_over_A0, rs$Delta_q, rs$Delta_p))

  q_extent <- max(abs(q_minus), abs(q_plus))
  q_pad    <- q_extent * 0.3
  q_lo     <- -(q_extent + q_pad)
  q_hi     <-   q_extent + q_pad

  V_min <- min(double_well_V(q_grid))
  p_max <- sqrt(2 * (E_n - V_min))
  p_lo  <- -(p_max * 1.3)
  p_hi  <-   p_max * 1.3

  custom_breaks_q <- round(c(q_minus, 0, q_plus), 1)
  custom_breaks_p <- round(c(-p_max, 0, p_max), 1)
  label_format    <- function(x) sprintf("%.1f", x)
  q_display       <- seq(q_lo, q_hi, length.out=500)

  state <- build_semiclassical_state(E_n, double_well_V,
                                     q_lo, q_hi, p_lo, p_hi, q_display)

  rho_husimi <- apply_husimi_marginal(state, q_display)

  rho_peak  <- max(rho_husimi, na.rm=TRUE)
  y_lim_rho <- rho_peak / (GOLDEN_FILL + 0.2)

  if (is_sub_barrier) {
    finite_caustic <- state$wkb_density[is.finite(state$wkb_density) &
                                          state$wkb_density > 0]
    caustic_floor <- if (length(finite_caustic) > 0)
      min(finite_caustic, na.rm=TRUE) else 1
  } else {
    q0_idx <- which.min(abs(q_display))
    caustic_floor <- state$wkb_density[q0_idx]
    if (!is.finite(caustic_floor) || caustic_floor <= 0) {
      finite_caustic <- state$wkb_density[is.finite(state$wkb_density) &
                                            state$wkb_density > 0]
      caustic_floor <- if (length(finite_caustic) > 0)
        min(finite_caustic, na.rm=TRUE) else 1
    }
  }
  y_lim_caustic <- caustic_floor / (1 - GOLDEN_FILL)

  dt_caustic <- data.table(q=q_display, wkb_density=state$wkb_density)
  dt_rho     <- data.table(q=q_display, rho_husimi=rho_husimi)

  overlay_layers <- husimi_overlay_layers(q_center=0)

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
    plot_semiclassical_husimi_resolution(
      dt_rho, q_lim=c(q_lo,q_hi), y_lim=y_lim_rho,
      custom_breaks=custom_breaks_q,
      label_format=label_format, base_font=base_font)
  )
}

cat("Computing double-well semiclassical-Husimi grid...\n")
rows    <- lapply(target_n_levels,
                  function(n) build_double_well_row(n, dw_soln, base_font=latex_font))
p_final <- assemble_grid(rows,
                         title_center=COLUMN_TITLE_CENTER_SEMICLASSICAL,
                         title_right=COLUMN_TITLE_RIGHT_HUSIMI,
                         base_font=latex_font)
save_figure(p_final, file_output_pdf, length(target_n_levels))
