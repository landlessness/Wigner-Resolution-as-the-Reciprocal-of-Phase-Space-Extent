# ==============================================================================
# plot_wigner.R
# Symplectic resolution of Wigner negativity across three systems.
# Each row shows one system at one chosen Schroedinger eigenstate:
#   row 1: harmonic    n=0   (ground state, analytical psi)
#   row 2: Morse       n=8   (highly excited bound state, multi-lobed Wigner)
#   row 3: double-well n=1   (antisymmetric lower-doublet partner with exact
#                             wavefunction node psi_1(0)=0; the deepest
#                             Wigner negativity coincides with the node)
#
# Three columns:
#   left:    bare Wigner heatmap with the symplectic quantum of action
#            (Fermi blob + conjugate quantum blobs) overlaid. Widths come
#            from the Robertson-Schroedinger covariance of psi_n itself,
#            not from any classical orbit.
#   center:  raw Wigner cross-section W_n(q, 0) at p=0 with negativity
#            visible as the curve dipping below zero.
#   right:   symplectic resolution P_{delta q}(q, 0) at p=0 (the Wigner
#            convolved with the symplectic kernel G_{delta q}, sliced at
#            p=0). Husimi cross-section Q(q, 0) overlaid as the prior-art
#            quantum comparator.
#
# QUANTUM UNIVERSE -- WHAT THIS PIPELINE TOUCHES (and does NOT touch):
#
#   This pipeline lives entirely in the quantum universe:
#     * potential V(q) and Schroedinger eigenstate psi_n(q)
#     * Wigner function W_n(q, p) via Fourier transform of psi_n
#     * Robertson-Schroedinger covariance of psi_n: <q^2>_psi, <p^2>_psi
#       compute Delta_q^RS, Delta_p^RS via numerical_covariance(psi)
#     * action capacity A_RS = pi * Delta_q^RS * Delta_p^RS (psi-derived)
#     * kernel widths delta_q^RS = hbar/Delta_p^RS, delta_p^RS = hbar/Delta_q^RS
#     * symplectic kernel G_{delta q} convolved with W_n
#     * Husimi kernel for the prior-art comparator overlay
#
#   This pipeline does NOT touch:
#     * classical orbit at energy E (no classical_trajectory)
#     * orbit time-averaged moments (no orbit_covariance)
#     * oscillating WKB densities (no oscillating_wkb_density,
#       no wkb_phase_space_lift, no build_semiclassical_state)
#     * Airy/Langer uniform construction (no airy_uniform_density)
#     * Bohr-Sommerfeld classical action (no classical_action)
#
#   The single shared object across the two universes is the symplectic
#   kernel G_{delta q} itself: a Gaussian of phase-space area h/2. Its
#   widths in this file come from psi_n's Robertson-Schroedinger
#   covariance -- NOT from any classical orbit. The semiclassical figure
#   uses the same Gaussian shape but derives its widths from the orbit.
#
# Author: Brian S. Mulloy
# ==============================================================================

library(here)
library(data.table)
library(patchwork)

# Plot helpers (shared across both universes; no quantum-state machinery here).
source(here("R", "plot_tools.R"))

# Quantum-universe modules ONLY. No classical_action_tools, no
# semiclassical_density, no airy_uniform.
source(here("R", "harmonic_system.R"))           # harmonic_V, harmonic_psi
                                                 # (analytical psi_n; turning
                                                 # points used only for axis
                                                 # extents, no classical action)
source(here("R", "morse_system.R"))              # morse_V, morse_E_BS,
                                                 # MORSE_Q_MIN/MAX/DQ
source(here("R", "double_well_system.R"))        # double_well_V, DOUBLE_WELL_*
source(here("R", "schroedinger_solver.R"))       # solve_schroedinger
source(here("R", "quantum_tools.R"))             # numerical_covariance
                                                 # (Robertson-Schroedinger
                                                 #  covariance of psi)
source(here("R", "wigner_density.R"))            # build_wigner_state,
                                                 # apply_kernel_cross_section,
                                                 # apply_kernel_density
source(here("R", "symplectic_kernel.R"))         # G_delta_q_kernel_matrix,
                                                 # symplectic_overlay_layers
source(here("R", "husimi_kernel.R"))             # husimi_kernel_matrix,
                                                 # husimi_marginal_density

latex_font  <- "CMU Serif"
dir_figures <- here("figures")
if (!dir.exists(dir_figures)) dir.create(dir_figures, recursive=TRUE)
file_output_pdf <- file.path(dir_figures, "wigner.pdf")

# ------------------------------------------------------------------------------
# WAVEFUNCTION SOURCES (quantum universe)
#
# Each system needs a function that returns psi_n(q) sampled on a grid.
# Harmonic uses the analytical Hermite-Gauss formula; Morse and double-
# well need numerical Schroedinger eigenfunctions. Solve once per system
# at module load and cache the soln object.
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
# SYSTEM DESCRIPTORS (quantum universe)
#
# Each descriptor exposes:
#   name        -- short identifier for console output
#   V           -- potential V(q) (universe-shared; just nature)
#   n_target    -- chosen Schroedinger quantum number for this row
#   E_fn        -- function(n) -> Schroedinger eigenvalue E_n
#   psi_fn      -- function(n, q) -> psi_n(q) sampled on q
#   q_window    -- function(E) -> list(q_lo, q_hi) for the display
#   p_window    -- function(E) -> list(p_lo, p_hi) for the display
#   q_breaks_fn -- function(E) -> custom q-axis breaks
#   p_breaks_fn -- function(E) -> custom p-axis breaks
#
# Note: q_window etc. take an energy E (or a wavefunction) so they can
# size the display window appropriately. For the quantum universe we
# use the wavefunction's classical turning points only as a heuristic
# for the axis extents -- the WAVEFUNCTION not the orbit drives the
# physics; the turning points are just a convenient axis-sizing rule.
# ------------------------------------------------------------------------------

# ---- Harmonic ----------------------------------------------------------------
harmonic_descriptor <- list(
  name = "harmonic",
  V    = harmonic_V,

  n_target = 0,

  E_fn = function(n) n + 0.5,   # exact for harmonic; analytical eigenvalue.

  psi_fn = function(n, q) {
    # harmonic_psi(n, q) returns the n-th Hermite-Gauss eigenstate sampled
    # on q. Real-valued; analytical.
    harmonic_psi(n, q)
  },

  q_window = function(E) {
    qt   <- sqrt(2 * E)
    span <- 2 * qt
    list(q_lo = -qt - 0.3*span/2, q_hi = qt + 0.3*span/2)
  },
  p_window = function(E) {
    p_max <- sqrt(2 * E)
    list(p_lo = -1.3*p_max, p_hi = 1.3*p_max)
  },
  q_breaks_fn = function(E) {
    qt <- sqrt(2 * E)
    round(c(-qt, qt), 1)
  },
  p_breaks_fn = function(E) {
    p_max <- sqrt(2 * E)
    round(c(-p_max, 0, p_max), 1)
  },

  # Wavefunction sampling grid for the FFT. Use a wide window so psi
  # decays to zero at the boundaries and the FFT doesn't ring.
  psi_q_grid = seq(-25, 25, by=0.02)
)

# ---- Morse -------------------------------------------------------------------
morse_descriptor <- list(
  name = "morse",
  V    = morse_V,

  n_target = 8,

  E_fn = function(n) morse_soln$energies[n + 1],

  psi_fn = function(n, q) {
    # Linearly interpolate the Schroedinger solver's psi onto the
    # requested grid. The solver was run on MORSE_Q_MIN..MORSE_Q_MAX
    # with spacing MORSE_DQ; we re-sample to whatever q is passed in.
    psi_solver <- morse_soln$psi_matrix[, n + 1]
    psi_q      <- approx(morse_soln$q_grid, psi_solver, xout=q,
                         rule=2, yleft=0, yright=0)$y
    psi_q[is.na(psi_q)] <- 0
    psi_q
  },

  q_window = function(E) {
    qm   <- -log(1 + sqrt(E/De)) / alpha
    qp   <- -log(1 - sqrt(E/De)) / alpha
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
    p_max <- sqrt(2 * E)
    round(c(-p_max, 0, p_max), 1)
  },

  psi_q_grid = morse_soln$q_grid
)

# ---- Double well -------------------------------------------------------------
double_well_descriptor <- list(
  name = "double_well",
  V    = double_well_V,

  # Schroedinger n=1: antisymmetric lower-doublet partner, E = -2.65175.
  # Has an exact node at q=0 (psi_1(0)=0); after symplectic resolution,
  # P_{delta q}(q=0, p=0) is positive (~5e-6), demonstrating the title's
  # "Schroedinger nodes resolved via the uncertainty principle" claim.
  # Aligns with the semiclassical figure's classical action target
  # A_BS/A_0 = 1.97593 at the same energy.
  n_target = 1,

  E_fn = function(n) dw_soln$energies[n + 1],

  psi_fn = function(n, q) {
    psi_solver <- dw_soln$psi_matrix[, n + 1]
    psi_q      <- approx(dw_soln$q_grid, psi_solver, xout=q,
                         rule=2, yleft=0, yright=0)$y
    psi_q[is.na(psi_q)] <- 0
    psi_q
  },

  q_window = function(E) {
    # Compute classical turning points (sized-for-display only; not used
    # to drive any physics in the quantum-universe pipeline).
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
# GENERIC WIGNER ROW BUILDER
#
# Pure quantum-universe pipeline. Takes a system descriptor, returns the
# three-panel Wigner row.
# ------------------------------------------------------------------------------

build_wigner_row <- function(descriptor, base_font="") {
  n   <- descriptor$n_target
  E_n <- descriptor$E_fn(n)

  cat(sprintf("\n== %s | n=%d | E_n=%.4f ==\n",
              descriptor$name, n, E_n))

  # 1. Sample psi_n on the wavefunction grid.
  psi_q_grid <- descriptor$psi_q_grid
  psi_vec    <- descriptor$psi_fn(n, psi_q_grid)

  # 2. Robertson-Schroedinger covariance from the wavefunction.
  #    This is the quantum-universe-honest source of kernel widths.
  rs <- numerical_covariance(psi_vec, psi_q_grid, hbar=1.0)
  cat(sprintf("  RS-covariance: A_RS/A0=%.4f Delta_q=%.3f Delta_p=%.3f\n",
              rs$A_over_A0, rs$Delta_q, rs$Delta_p))
  cat(sprintf("                 delta_q=%.3f delta_p=%.3f\n",
              rs$delta_q, rs$delta_p))

  # 3. Display windows.
  qw <- descriptor$q_window(E_n)
  pw <- descriptor$p_window(E_n)
  q_lo <- qw$q_lo; q_hi <- qw$q_hi
  p_lo <- pw$p_lo; p_hi <- pw$p_hi

  custom_breaks_q <- descriptor$q_breaks_fn(E_n)
  custom_breaks_p <- descriptor$p_breaks_fn(E_n)
  label_format    <- function(x) sprintf("%.1f", x)
  q_display       <- seq(q_lo, q_hi, length.out=500)

  # 4. Wigner state. build_wigner_state computes W via FFT of psi.
  cat("  Building Wigner state...\n")
  state <- build_wigner_state(
    psi_vec=psi_vec, psi_q_grid=psi_q_grid,
    q_lo=q_lo, q_hi=q_hi, p_lo=p_lo, p_hi=p_hi,
    q_display=q_display)

  # 5. Symplectic kernel: shape Gaussian, area h/2; widths from RS.
  symplectic_kernel_for_state <- function(qg, pg) {
    G_delta_q_kernel_matrix(qg, pg, rs$Delta_q, rs$Delta_p)
  }
  cat("  Convolving with symplectic kernel...\n")
  P_sympl_cross <- apply_kernel_cross_section(state, symplectic_kernel_for_state,
                                              q_display)

  # 6. Husimi: shape Gaussian, fixed width. Quantum-universe comparator.
  #    Use the Husimi cross-section at p=0 (NOT the marginal over p) so
  #    it is directly comparable to the symplectic cross-section -- both
  #    are 1D slices of their respective convolved 2D densities at p=0.
  husimi_kernel_for_state <- function(qg, pg) {
    husimi_kernel_matrix(qg, pg)
  }
  cat("  Computing Husimi cross-section at p=0...\n")
  Q_husimi_cross <- apply_kernel_cross_section(state, husimi_kernel_for_state,
                                               q_display)

  # 7. Y-scaling.
  W_cross_peak    <- max(abs(state$W_cross), na.rm=TRUE)
  if (!is.finite(W_cross_peak) || W_cross_peak == 0) W_cross_peak <- 1
  y_lim_W <- W_cross_peak * 1.1

  P_peak_data <- max(P_sympl_cross, na.rm=TRUE)
  if (!is.finite(P_peak_data) || P_peak_data == 0) P_peak_data <- 1
  y_lim_P <- P_peak_data * 1.1

  # 8. Build data.tables for the cross-section panels.
  dt_W       <- data.table(q=q_display, W_raw=state$W_cross)
  dt_P_sympl <- data.table(q=q_display, rho_sympl=P_sympl_cross)
  dt_Husimi  <- data.frame(q=q_display, rho=Q_husimi_cross)

  # 9. QoA overlay. RS-derived widths, centered at the wavefunction's
  #    <q>. For symmetric states <q> = 0; for asymmetric states like
  #    Morse n=8, <q> can be far from the well bottom.
  overlay_layers <- symplectic_overlay_layers(rs$Delta_q, rs$Delta_p,
                                              q_center=rs$q_mean)

  # 10. Husimi overlay on the symplectic cross-section panel.
  husimi_overlay <- list(
    list(
      data       = dt_Husimi,
      color      = "gray30",
      linewidth  = 0.35,
      fill       = "gray85",
      fill_alpha = 0.5
    )
  )

  # 11. Build the three panels.
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

cat("\nComputing Wigner figure (3 systems x 3 panels)...\n")

descriptors <- list(harmonic_descriptor,
                    morse_descriptor,
                    double_well_descriptor)

rows <- lapply(descriptors,
               function(d) build_wigner_row(d, base_font=latex_font))

p_final <- assemble_grid_unlabeled(rows,
                                   COLUMN_TITLE_CENTER_WIGNER,
                                   COLUMN_TITLE_RIGHT_SYMPLECTIC,
                                   base_font=latex_font)

save_figure(p_final, file_output_pdf, length(descriptors))
