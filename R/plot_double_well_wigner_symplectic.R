# ==============================================================================
# plot_double_well_wigner_symplectic.R
# Main paper: Symplectic resolution of Wigner negativity — Double-well potential
#
# Same three rows as the double-well Husimi appendix figure:
#   n=0  ground-state symmetric tunneling doublet partner
#   n=2  upper symmetric doublet partner (still below barrier)
#   n=7  above-barrier state (single connected orbit)
#
# Same Wigner heatmap (left column) and Wigner cross-section (middle column)
# as the Husimi figure. The differences are entirely in the kernel:
#   - Left column overlay: three QoA ellipses (outer A, inner a_q, a_p)
#     instead of a single Husimi unit circle
#   - Right column: P_delta_q(q,0) = (W * G_delta_q)(q,0) in place of Q(q,0)
#
# Architecture:
#   build_wigner_state()              kernel-agnostic per-state computation
#   apply_kernel_cross_section()      kernel-specific (symplectic G_delta_q)
#   G_delta_q_kernel_matrix(),
#   symplectic_overlay_layers()       the two symplectic-specific pieces
# ==============================================================================

library(here)
library(patchwork)

source(here("R", "plot_tools.R"))
source(here("R", "double_well_potential.R"))
source(here("R", "wigner_tools.R"))
source(here("R", "wigner_state.R"))
source(here("R", "symplectic_kernel.R"))
source(here("R", "classical_action_tools.R"))

latex_font      <- "CMU Serif"
dir_figures     <- here("figures")
if (!dir.exists(dir_figures)) dir.create(dir_figures, recursive=TRUE)
file_output_pdf <- file.path(dir_figures, "double_well_wigner_symplectic.pdf")

cat("Solving double-well Schrodinger equation...\n")
dw_soln <- solve_schrodinger(double_well_V,
                             q_min=DOUBLE_WELL_Q_MIN,
                             q_max=DOUBLE_WELL_Q_MAX,
                             dq=DOUBLE_WELL_DQ,
                             n_states=DOUBLE_WELL_N_STATES)

target_n_levels <- c(0, 2, 7)

# ------------------------------------------------------------------------------

build_double_well_row <- function(n_val, soln, base_font="") {
  E_n     <- soln$energies[n_val+1]
  q_grid  <- soln$q_grid
  psi_vec <- soln$psi_matrix[, n_val+1]
  tp      <- double_well_turning_points(E_n)$roots

  cat(sprintf("\n== n=%d | E_n=%+.4f | turning points: %s ==\n",
              n_val, E_n,
              paste(sprintf("%+.3f", tp), collapse=", ")))

  rs   <- numerical_covariance(psi_vec, q_grid)
  A_BS <- classical_action(double_well_V, E_n, tp)

  cat(sprintf("  A_BS/A0=%.2f | A_RS/A0=%.2f | Delta_q=%.3f Delta_p=%.3f\n",
              A_BS, rs$A_over_A0, rs$Delta_q, rs$Delta_p))

  q_orbit_min <- min(tp)
  q_orbit_max <- max(tp)
  q_pad       <- (q_orbit_max - q_orbit_min) * 0.3
  q_lo        <- q_orbit_min - q_pad
  q_hi        <- q_orbit_max + q_pad
  p_max       <- sqrt(2*(E_n - min(double_well_V(seq(q_lo, q_hi, length.out=500)))))
  p_lo        <- -(p_max * 1.3)
  p_hi        <-   p_max * 1.3

  custom_breaks_q <- round(c(min(tp), 0, max(tp)), 1)
  custom_breaks_p <- round(c(-p_max, 0, p_max), 1)
  label_format    <- function(x) sprintf("%.1f", x)
  q_display       <- seq(q_lo, q_hi, length.out=500)

  state <- build_wigner_state(psi_vec, q_grid,
                              q_lo, q_hi, p_lo, p_hi, q_display)

  # Symplectic kernel for this state — closure binds Delta_q, Delta_p.
  symplectic_kernel_for_state <- function(q_grid, p_grid) {
    G_delta_q_kernel_matrix(q_grid, p_grid, rs$Delta_q, rs$Delta_p)
  }
  P_cross <- apply_kernel_cross_section(state, symplectic_kernel_for_state, q_display)

  w_max     <- max(abs(state$W_cross), na.rm=TRUE)
  y_lim     <- w_max * 1.3
  p_max_amp <- max(abs(P_cross), na.rm=TRUE)
  P_display <- if (p_max_amp > 0) P_cross/p_max_amp*w_max else P_cross
  dt_cross  <- data.table(q=q_display, W_raw=state$W_cross, P_sympl=P_display)

  df_traj        <- classical_trajectory(double_well_V, E_n, tp)
  overlay_layers <- symplectic_overlay_layers(rs$Delta_q, rs$Delta_p, q_center=0)

  list(
    sprintf("italic(n)==%d", n_val),
    plot_wigner_heatmap(
      state$heatmap_dt, overlay_layers, df_traj,
      q_lim=c(q_lo,q_hi), p_lim=c(p_lo,p_hi),
      custom_breaks_q=custom_breaks_q,
      custom_breaks_p=custom_breaks_p,
      label_format=label_format, base_font=base_font),
    plot_wigner_cross_section(
      dt_cross, q_lim=c(q_lo,q_hi), y_lim=y_lim,
      custom_breaks=custom_breaks_q,
      label_format=label_format, base_font=base_font),
    plot_symplectic_cross_section(
      dt_cross, q_lim=c(q_lo,q_hi), y_lim=y_lim,
      custom_breaks=custom_breaks_q,
      label_format=label_format, base_font=base_font)
  )
}

cat("Computing double-well Wigner-Symplectic grid...\n")
rows    <- lapply(target_n_levels,
                  function(n) build_double_well_row(n, dw_soln, base_font=latex_font))

p_final <- assemble_grid(rows,
                         title_center=COLUMN_TITLE_CENTER_WIGNER,
                         title_right=COLUMN_TITLE_RIGHT_SYMPLECTIC,
                         base_font=latex_font)

save_figure(p_final, file_output_pdf, length(target_n_levels))
