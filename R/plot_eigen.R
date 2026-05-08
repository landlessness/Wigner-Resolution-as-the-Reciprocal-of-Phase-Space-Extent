# ==============================================================================
# plot_wigner.R
# Symplectic resolution of Wigner negativity for bound-state systems.
# Each row shows one quantum eigenstate:
#   row 1: squeezed vacuum     (Gaussian; non-negative W; demonstrates that
#                               the symplectic kernel matches the state's
#                               anisotropic geometry while Husimi does not)
#   row 2: harmonic n=1        (Fock state; central negative dip resolved
#                               as positive density at delta_q)
#   row 3: Morse n=8           (anharmonic; high-n WKB-like ringing)
#   row 4: double-well n=5     (delocalized; central-barrier interference)
#
# Three columns:
#   left:    bare Wigner heatmap with symplectic overlay (Fermi blob plus
#            conjugate quantum blobs).
#   center:  raw Wigner cross-section W(q, 0) at p=0.
#   right:   symplectic resolution P_{delta q}(q, 0) at p=0.
#
# State-building (psi, RS covariance, Wigner via FFT, symplectic overlay)
# is delegated to build_eigenstate_state() in state_builder.R, so the
# tomography pipeline can consume the same per-state output.
#
# Author: Brian S. Mulloy
# ==============================================================================

library(here)
library(data.table)
library(patchwork)

source(here("R", "plot_tools.R"))

# Quantum-universe modules.
source(here("R", "harmonic_system.R"))           # harmonic_V, harmonic_psi
source(here("R", "morse_system.R"))              # morse_V, MORSE_*
source(here("R", "double_well_system.R"))        # double_well_V, DOUBLE_WELL_*
source(here("R", "schroedinger_solver.R"))       # solve_schroedinger
source(here("R", "wigner_density.R"))            # apply_kernel_cross_section
source(here("R", "symplectic_kernel.R"))         # G_delta_q_kernel_matrix
source(here("R", "husimi_kernel.R"))             # husimi_kernel_matrix
source(here("R", "state_builder.R"))             # build_eigenstate_state

latex_font  <- "CMU Serif"
dir_figures <- here("figures")
if (!dir.exists(dir_figures)) dir.create(dir_figures, recursive=TRUE)
file_output_pdf <- file.path(dir_figures, "wigner.pdf")

# ------------------------------------------------------------------------------
# WAVEFUNCTION SOURCES (quantum universe)
# ------------------------------------------------------------------------------

cat("Solving Schroedinger for Morse...\n")
suppressMessages({
  morse_soln <- solve_schroedinger(
    morse_V, MORSE_Q_MIN, MORSE_Q_MAX, MORSE_DQ,
    n_states=MORSE_N_STATES)
})

cat("Solving Schroedinger for double well...\n")
suppressMessages({
  dw_soln <- solve_schroedinger(
    double_well_V, DOUBLE_WELL_Q_MIN, DOUBLE_WELL_Q_MAX,
    DOUBLE_WELL_DQ, n_states=DOUBLE_WELL_N_STATES)
})

# ------------------------------------------------------------------------------
# SYSTEM DESCRIPTORS
# ------------------------------------------------------------------------------

# ---- Harmonic ----------------------------------------------------------------
harmonic_descriptor <- list(
  name = "harmonic",
  V    = harmonic_V,
  n_target = 1,
  E_fn = function(n) n + 0.5,
  psi_fn = function(n, q) harmonic_psi(n, q),
  q_window = function(E) {
    qt <- sqrt(2 * E); span <- 2 * qt
    list(q_lo = -qt - 0.3*span/2, q_hi = qt + 0.3*span/2)
  },
  p_window = function(E) {
    p_max <- sqrt(2 * E)
    list(p_lo = -1.3*p_max, p_hi = 1.3*p_max)
  },
  q_breaks_fn = function(E) {
    qt <- sqrt(2 * E); round(c(-qt, qt), 1)
  },
  p_breaks_fn = function(E) {
    p_max <- sqrt(2 * E); round(c(-p_max, 0, p_max), 1)
  },
  psi_q_grid = seq(-25, 25, by=0.02)
)

# ---- Squeezed vacuum ---------------------------------------------------------
SQUEEZED_R <- 0.5
squeezed_vacuum_descriptor <- list(
  name = "squeezed_vacuum",
  V    = harmonic_V,
  n_target = 0,
  E_fn = function(n) 0.5 * cosh(2*SQUEEZED_R),
  psi_fn = function(n, q) {
    sigma_q <- exp(-SQUEEZED_R) / sqrt(2)
    norm    <- (1 / (pi * sigma_q^2))^(1/4)
    norm * exp(-q^2 / (2 * sigma_q^2))
  },
  q_window = function(E) {
    Delta_q <- exp(-SQUEEZED_R)
    list(q_lo = -3 * Delta_q, q_hi = 3 * Delta_q)
  },
  p_window = function(E) {
    Delta_p <- exp(+SQUEEZED_R)
    list(p_lo = -1.5 * Delta_p, p_hi = 1.5 * Delta_p)
  },
  q_breaks_fn = function(E) {
    Delta_q <- exp(-SQUEEZED_R); round(c(-Delta_q, Delta_q), 2)
  },
  p_breaks_fn = function(E) {
    Delta_p <- exp(+SQUEEZED_R); round(c(-Delta_p, 0, Delta_p), 1)
  },
  psi_q_grid = seq(-15, 15, by=0.005)
)

# ---- Morse -------------------------------------------------------------------
morse_descriptor <- list(
  name = "morse",
  V    = morse_V,
  n_target = 8,
  E_fn = function(n) morse_soln$energies[n + 1],
  psi_fn = function(n, q) {
    psi_solver <- morse_soln$psi_matrix[, n + 1]
    psi_q      <- approx(morse_soln$q_grid, psi_solver, xout=q,
                         rule=2, yleft=0, yright=0)$y
    psi_q[is.na(psi_q)] <- 0
    psi_q
  },
  q_window = function(E) {
    qm <- -log(1 + sqrt(E/De)) / alpha
    qp <- -log(1 - sqrt(E/De)) / alpha
    span <- qp - qm
    list(q_lo = qm - 0.15*span, q_hi = qp + 0.15*span)
  },
  p_window = function(E) {
    p_max <- sqrt(2 * E)
    list(p_lo = -1.3*p_max, p_hi = 1.3*p_max)
  },
  q_breaks_fn = function(E) {
    qm <- -log(1 + sqrt(E/De)) / alpha
    qp <- -log(1 - sqrt(E/De)) / alpha
    round(c(qm, qp), 1)
  },
  p_breaks_fn = function(E) {
    p_max <- sqrt(2 * E); round(c(-p_max, 0, p_max), 1)
  },
  psi_q_grid = morse_soln$q_grid
)

# ---- Double well -------------------------------------------------------------
double_well_descriptor <- list(
  name = "double_well",
  V    = double_well_V,
  n_target = 5,
  E_fn = function(n) dw_soln$energies[n + 1],
  psi_fn = function(n, q) {
    psi_solver <- dw_soln$psi_matrix[, n + 1]
    psi_q      <- approx(dw_soln$q_grid, psi_solver, xout=q,
                         rule=2, yleft=0, yright=0)$y
    psi_q[is.na(psi_q)] <- 0
    psi_q
  },
  q_window = function(E) {
    roots <- polyroot(c(-E, 0, -mu2/2, 0, lambda/4))
    real_roots <- sort(Re(roots[abs(Im(roots)) < 1e-8]))
    q_lo <- min(real_roots); q_hi <- max(real_roots)
    span <- q_hi - q_lo
    list(q_lo = q_lo - 0.2*span, q_hi = q_hi + 0.2*span)
  },
  p_window = function(E) {
    V_min <- -double_well_barrier
    p_max <- sqrt(2 * (E - V_min))
    list(p_lo = -1.3*p_max, p_hi = 1.3*p_max)
  },
  q_breaks_fn = function(E) {
    roots <- polyroot(c(-E, 0, -mu2/2, 0, lambda/4))
    real_roots <- sort(Re(roots[abs(Im(roots)) < 1e-8]))
    round(c(min(real_roots), max(real_roots)), 1)
  },
  p_breaks_fn = function(E) {
    V_min <- -double_well_barrier
    p_max <- sqrt(2 * (E - V_min))
    round(c(-p_max, 0, p_max), 1)
  },
  psi_q_grid = dw_soln$q_grid
)

# ------------------------------------------------------------------------------
# ROW BUILDER
#
# Calls build_eigenstate_state() for the shared per-state pipeline, then
# adds the per-row Husimi cross-section and symplectic cross-section that
# this figure needs (and the tomography pipeline does not).
# ------------------------------------------------------------------------------

build_wigner_row <- function(descriptor, base_font="") {
  ps <- build_eigenstate_state(descriptor, base_font=base_font)

  # Symplectic kernel cross-section for the right column.
  symplectic_kernel_for_state <- function(qg, pg) {
    G_delta_q_kernel_matrix(qg, pg, ps$rs$Delta_q, ps$rs$Delta_p,
                            hbar=ps$hbar)
  }
  cat("  Convolving with symplectic kernel...\n")
  P_sympl_cross <- apply_kernel_cross_section(ps$state,
                                              symplectic_kernel_for_state,
                                              ps$q_display)

  # Husimi cross-section. (Computed even though the overlay is currently
  # disabled in plot_semiclassical_resolution(); cheap, leaves room to
  # toggle back on later.)
  husimi_kernel_for_state <- function(qg, pg) {
    husimi_kernel_matrix(qg, pg)
  }
  cat("  Computing Husimi cross-section at p=0...\n")
  Q_husimi_cross <- apply_kernel_cross_section(ps$state,
                                               husimi_kernel_for_state,
                                               ps$q_display)

  # Y-scaling.
  W_cross_peak <- max(abs(ps$state$W_cross), na.rm=TRUE)
  if (!is.finite(W_cross_peak) || W_cross_peak == 0) W_cross_peak <- 1
  y_lim_W <- W_cross_peak * 1.1

  P_peak_data <- max(P_sympl_cross, na.rm=TRUE)
  if (!is.finite(P_peak_data) || P_peak_data == 0) P_peak_data <- 1
  y_lim_P <- P_peak_data * 1.1

  # Cross-section data tables.
  dt_W       <- data.table(q=ps$q_display, W_raw=ps$state$W_cross)
  dt_P_sympl <- data.table(q=ps$q_display, rho_sympl=P_sympl_cross)

  # Three panels.
  list(
    plot_wigner_heatmap(
      ps$state$heatmap_dt, ps$overlay_layers, df_traj=NULL,
      q_lim=c(ps$q_lo, ps$q_hi), p_lim=c(ps$p_lo, ps$p_hi),
      custom_breaks_q=ps$custom_breaks_q,
      custom_breaks_p=ps$custom_breaks_p,
      label_format=ps$label_format, base_font=base_font),
    plot_wigner_cross_section(
      dt_W, q_lim=c(ps$q_lo, ps$q_hi), y_lim=y_lim_W,
      custom_breaks=ps$custom_breaks_q,
      label_format=ps$label_format, base_font=base_font),
    plot_semiclassical_resolution(
      dt_P_sympl, q_lim=c(ps$q_lo, ps$q_hi), y_lim=y_lim_P,
      custom_breaks=ps$custom_breaks_q,
      label_format=ps$label_format, base_font=base_font,
      overlays=NULL,
      y_label=expression(italic(P)[italic(delta*q)](italic(q)*","*0)))
  )
}

# ------------------------------------------------------------------------------
# DRIVE
# ------------------------------------------------------------------------------

cat("\nComputing Wigner figure (4 systems x 3 panels)...\n")

descriptors <- list(squeezed_vacuum_descriptor,
                    harmonic_descriptor,
                    morse_descriptor,
                    double_well_descriptor)

rows <- lapply(descriptors,
               function(d) build_wigner_row(d, base_font=latex_font))

p_final <- assemble_grid_unlabeled(rows,
                                   COLUMN_TITLE_CENTER_WIGNER,
                                   COLUMN_TITLE_RIGHT_SYMPLECTIC,
                                   base_font=latex_font)

save_figure(p_final, file_output_pdf, length(descriptors))
