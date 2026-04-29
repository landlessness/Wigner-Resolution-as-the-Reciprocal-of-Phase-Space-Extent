# ==============================================================================
# semiclassical_density.R
# Semiclassical phase-space densities and the per-state bundle that wraps them.
#
# The semiclassical regime takes V(q) and energy E_n as inputs and produces
# phase-space content WITHOUT solving Schroedinger:
#
#   wkb_caustic_density(q, E, V_fn): 1D analytical density
#     P_mc(q) = 1 / (T(E) * |p(q)|), where p(q) = sqrt(2*(E - V(q)))
#     and T(E) is the classical orbit period:
#       T(E) = oint dq / |p(q)|
#     Diverges at classical turning points where E = V(q) — the caustic.
#
#   semiclassical_shell_density(q, p, E, V_fn, epsilon): regularized 2D shell
#     W_cl(q,p) = (1/Z) * exp(-(H(q,p) - E)^2 / (2*epsilon^2))
#     where H(q,p) = p^2/2 + V(q). The regularization is purely numerical
#     (the discrete grid cannot represent a Dirac delta); epsilon is set
#     small relative to E_n and the kernel widths so that downstream
#     convolution is insensitive to it.
#
#   build_semiclassical_state(...): orchestrates the above into a per-state
#     bundle compatible in shape with build_wigner_state(), so kernel
#     application code can treat them uniformly.
#
# Kernel-sizing widths (Delta_q, Delta_p) for the symplectic kernel come
# from orbit_covariance() in classical_action_tools.R — also Schroedinger-free.
# The whole semiclassical pipeline is independent of the wavefunction.
#
# Reference: Berry 1977 J.Phys.A 10, 2083 (semiclassical phase space)
#            Berry & Mount 1972 Rep.Prog.Phys. 35, 315 (WKB caustics)
# Author: Brian S. Mulloy
# ==============================================================================

library(here)
library(data.table)
source(here("R", "math_tools.R"))

# Default regularization for the 2D shell:
# epsilon = SEMICLASSICAL_REG_FRACTION * |E_n|. Small relative to the energy
# scale, but large enough to span several grid cells of the integration grid.
SEMICLASSICAL_REG_FRACTION <- 0.02

# ------------------------------------------------------------------------------
# 1D WKB CAUSTIC DENSITY
# ------------------------------------------------------------------------------

#' Analytical WKB density at energy E for potential V_fn on q_grid.
#'
#' Computes P_mc(q) = 1 / (T(E) * |p(q)|) where p(q) = sqrt(2*(E - V(q))).
#' Diverges at classical turning points where E = V(q); these points are
#' rendered as Inf in the returned vector. Caller is expected to handle the
#' infinities (e.g., clip in plot, render with arrows).
#'
#' @param q_grid Position grid (numeric vector).
#' @param E Energy.
#' @param V_fn Function(q) returning potential values.
#' @return Numeric vector of P_mc values, with Inf at turning points and 0
#'         in classically forbidden region E < V(q).
wkb_caustic_density <- function(q_grid, E, V_fn) {
  V_vals <- V_fn(q_grid)
  arg    <- 2 * (E - V_vals)
  p_q    <- sqrt(pmax(arg, 0))   # |p(q)| where defined; 0 in forbidden region

  # Compute period via fine-grid integration over the classically-allowed region.
  q_fine  <- seq(min(q_grid), max(q_grid), length.out=10001)
  V_fine  <- V_fn(q_fine)
  arg_f   <- 2 * (E - V_fine)
  allowed <- arg_f > 0
  # Tiny floor avoids 0/0; integrand falls off as 1/sqrt away from the
  # turning point, integrable on a fine grid with trapezoidal quadrature.
  p_fine    <- sqrt(pmax(arg_f, 1e-15))
  dq_fine   <- diff(q_fine)[1]
  integrand <- ifelse(allowed, 1/p_fine, 0)
  T_E       <- sum(integrand) * dq_fine

  # Caustic density: 1 / (T * |p|), with Inf at turning points (p=0).
  caustic <- ifelse(p_q > 0, 1/(T_E * p_q), Inf)
  caustic[!is.finite(caustic) & p_q == 0 & arg == 0] <- Inf  # turning points
  caustic[arg < 0] <- 0  # forbidden region: zero, not infinity

  caustic
}

# ------------------------------------------------------------------------------
# 2D REGULARIZED ENERGY SHELL
# ------------------------------------------------------------------------------

#' Regularized classical energy shell on a (q,p) grid.
#'
#' W_cl(q,p) = (1/Z) * exp(-(H(q,p) - E)^2 / (2*epsilon^2))
#' where H(q,p) = p^2/2 + V(q). Z normalizes so that integral W_cl dq dp = 1.
#'
#' Used as input to the symplectic and Husimi convolutions on the
#' semiclassical side. The regularization width epsilon should be small
#' relative to E (so the band is thin) but large enough that the band
#' spans at least 2-3 grid cells (so numerical FFT convolution is
#' well-behaved).
#'
#' @param q_grid Position grid.
#' @param p_grid Momentum grid.
#' @param E Energy.
#' @param V_fn Function(q) returning potential values.
#' @param epsilon Regularization width (energy units).
#' @return Matrix W_cl[iq, ip].
semiclassical_shell_density <- function(q_grid, p_grid, E, V_fn, epsilon) {
  V_vals <- V_fn(q_grid)
  H_mat  <- outer(V_vals, p_grid^2/2, FUN="+")
  W_mat  <- exp(-(H_mat - E)^2 / (2*epsilon^2))
  dq     <- diff(q_grid)[1]
  dp     <- diff(p_grid)[1]
  Z      <- sum(W_mat) * dq * dp
  if (Z > 0) W_mat <- W_mat / Z
  W_mat
}

# ------------------------------------------------------------------------------
# BUILD SEMICLASSICAL STATE
# ------------------------------------------------------------------------------

#' Build a per-state semiclassical bundle.
#'
#' Drop-in replacement for build_wigner_state(): returns a bundle with the
#' same field names (q_int, p_int, dq_int, dp_int, W_matrix, W_cross,
#' heatmap_dt, norm) so apply_kernel_cross_section() works without
#' modification.
#'
#' In addition, the bundle carries:
#'   wkb_density   — analytical 1D WKB density on q_display
#'                   (Inf at turning points, 0 in forbidden region)
#'
#' @param E_n Energy of the eigenstate.
#' @param V_fn Function(q) returning potential values.
#' @param q_lo,q_hi Display window in q.
#' @param p_lo,p_hi Display window in p.
#' @param q_display Position grid for cross-section / projection rendering.
#' @param epsilon Regularization width for 2D shell. If NULL, defaults to
#'                SEMICLASSICAL_REG_FRACTION * |E_n|.
#' @param nq_int,np_int Integration grid resolution (default 801 x 601).
#' @return List with fields matching build_wigner_state() bundle plus
#'         wkb_density.
build_semiclassical_state <- function(E_n, V_fn,
                                      q_lo, q_hi, p_lo, p_hi, q_display,
                                      epsilon=NULL,
                                      nq_int=801, np_int=601) {
  if (is.null(epsilon)) epsilon <- SEMICLASSICAL_REG_FRACTION * abs(E_n)

  q_int  <- seq(q_lo, q_hi, length.out=nq_int)
  p_int  <- seq(p_lo, p_hi, length.out=np_int)
  dq_int <- diff(q_int)[1]
  dp_int <- diff(p_int)[1]

  cat(sprintf("    Building semiclassical shell on %d x %d grid (epsilon=%.4f)...\n",
              nq_int, np_int, epsilon))

  W_mat <- semiclassical_shell_density(q_int, p_int, E_n, V_fn, epsilon)
  norm  <- sum(W_mat) * dq_int * dp_int
  cat(sprintf("    Shell norm: %.6f\n", norm))

  # Extract p=0 cross-section (interpolated to q_display)
  W_cross <- extract_p0_cross_section(W_mat, q_int, p_int, q_display)

  # 1D WKB caustic density on q_display
  wkb <- wkb_caustic_density(q_display, E_n, V_fn)

  # Heatmap data for ggplot
  heatmap_dt <- as.data.table(expand.grid(q=q_int, p=p_int))
  heatmap_dt[, w := as.vector(W_mat)]
  # Apply soft contrast curve: sqrt for visual emphasis of the shell band
  max_w <- max(heatmap_dt$w, na.rm=TRUE)
  if (max_w > 0) {
    heatmap_dt[, w_plot := sqrt(pmax(w/max_w, 0))]
  } else {
    heatmap_dt[, w_plot := 0]
  }

  list(
    q_int       = q_int,
    p_int       = p_int,
    dq_int      = dq_int,
    dp_int      = dp_int,
    W_matrix    = W_mat,
    W_cross     = W_cross,
    wkb_density = wkb,
    heatmap_dt  = heatmap_dt,
    norm        = norm,
    epsilon     = epsilon
  )
}
