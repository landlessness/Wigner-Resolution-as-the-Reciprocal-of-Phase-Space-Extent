# ==============================================================================
# semiclassical_density.R
# Semiclassical phase-space densities and caustic curves.
#
# Two related quantities:
#
#   wkb_caustic_density(q, E, V_fn): 1D analytical density
#     P_mc(q) = 1 / (T(E) * |p(q)|), where p(q) = sqrt(2*(E - V(q)))
#     and T(E) is the classical orbit period:
#       T(E) = integral_{q_-}^{q_+} dq / |p(q)|
#     This is the WKB probability density. Diverges at classical turning
#     points where E = V(q) — the caustic.
#
#   semiclassical_shell_density(q, p, E, V_fn, epsilon): regularized 2D shell
#     W_cl(q,p) = (1/Z) * exp(-(H(q,p) - E)^2 / (2*epsilon^2))
#     where H(q,p) = p^2/2 + V(q). The regularization is purely numerical
#     (the discrete grid cannot represent a Dirac delta); epsilon is set
#     small relative to E_n and the kernel widths so that the symplectic
#     convolution result is insensitive to it.
#
# Reference: Berry 1977 J.Phys.A 10, 2083 (semiclassical phase space)
#            Berry & Mount 1972 Rep.Prog.Phys. 35, 315 (WKB caustics)
# Author: Brian S. Mulloy
# ==============================================================================

# ------------------------------------------------------------------------------
# 1D WKB CAUSTIC DENSITY
# ------------------------------------------------------------------------------

#' Analytical WKB density at energy E for potential V_fn on q_grid.
#'
#' Computes P_mc(q) = 1 / (T(E) * |p(q)|) where p(q) = sqrt(2*(E - V(q))).
#' Diverges at classical turning points where E = V(q); these points are
#' rendered as Inf in the returned vector. Caller is expected to handle
#' the infinities (e.g., clip in plot, render with arrows).
#'
#' @param q_grid Position grid (numeric vector).
#' @param E Energy.
#' @param V_fn Function(q) returning potential values.
#' @return Numeric vector of P_mc values, with Inf at turning points and
#'         0 in classically forbidden region E < V(q).
wkb_caustic_density <- function(q_grid, E, V_fn) {
  V_vals <- V_fn(q_grid)
  arg    <- 2 * (E - V_vals)
  p_q    <- sqrt(pmax(arg, 0))   # |p(q)| where defined; 0 in forbidden region

  # Compute period via fine-grid integration over the classically-allowed region.
  q_fine <- seq(min(q_grid), max(q_grid), length.out=10001)
  V_fine <- V_fn(q_fine)
  arg_f  <- 2 * (E - V_fine)
  allowed <- arg_f > 0
  p_fine <- sqrt(pmax(arg_f, 1e-15))   # tiny floor avoids 0/0; integrand falls
  # off as 1/sqrt away from turning point
  dq_fine <- diff(q_fine)[1]
  # Integrate 1/|p| over allowed region. The integrand diverges at turning
  # points but is integrable (1/sqrt singularity); use trapezoidal as a
  # reasonable approximation.
  integrand <- ifelse(allowed, 1/p_fine, 0)
  T_E <- sum(integrand) * dq_fine

  # Caustic density: 1 / (T * |p|), with Inf at turning points (p=0).
  caustic <- ifelse(p_q > 0, 1/(T_E * p_q), Inf)
  caustic[!is.finite(caustic) & p_q == 0 & arg == 0] <- Inf  # mark turning points
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
#' Used as input to the symplectic convolution. The regularization width
#' epsilon should be small relative to E (so the band is thin) but large
#' enough that the band spans at least 2-3 grid cells (so numerical FFT
#' convolution is well-behaved).
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
