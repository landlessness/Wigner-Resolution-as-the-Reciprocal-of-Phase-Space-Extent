# ==============================================================================
# quantum_tools.R
# Quantum-state properties of psi: Robertson-Schroedinger covariance and
# the action-capacity geometry derived from it.
#
# This file is the quantum analogue of classical_action_tools.R. The orbit
# covariance there depends only on V(q) and E; the covariance here depends
# on the wavefunction psi itself.
#
# Two routes to the same kinematic envelope:
#   robertson_schroedinger(sigma_qq, sigma_pp, sigma_qp)  — from moments
#   numerical_covariance(psi_vec, q_grid)                 — from sampled psi
#
# Both return Delta_q, Delta_p, delta_q, delta_p, and the Fermi-blob ratio
# A/A_0, in the same structure used by orbit_covariance() so callers can
# interchange them.
#
# Reference: Robertson Phys. Rev. 34, 163 (1929)
#            Schroedinger Sitzungsber. Preuss. Akad. Wiss. Berlin 24, 296
#            (1930)
#            Zurek Nature 412, 712 (2001) [reciprocal scales]
# Author: Brian S. Mulloy
# ==============================================================================

# Tolerance for RS-bound saturation check: relative deviation that we
# accept as numerically equivalent. Analytic ground states saturate RS
# exactly; the discrete finite-difference computation may sit fractionally
# below the continuous bound by O(dq^2).
RS_TOLERANCE <- sqrt(.Machine$double.eps)

# ------------------------------------------------------------------------------
# RS GEOMETRY FROM COVARIANCE MATRIX
# ------------------------------------------------------------------------------

#' Robertson-Schroedinger geometry from covariance matrix elements.
#'
#' Given the covariance matrix entries, returns the kinematic widths
#' (Delta_q, Delta_p), Zurek's reciprocal scales (delta_q, delta_p),
#' the Fermi-blob area in units of A_0 = h/2, and diagnostic flags.
#'
#' Note: the inequality sigma_qq * sigma_pp - sigma_qp^2 >= (hbar/2)^2 is
#' the Robertson-Schroedinger generalization of Heisenberg; this routine
#' issues warnings if the input covariance violates the bound, which would
#' indicate a non-physical psi or numerical-precision issues.
#'
#' @param sigma_qq Variance in q.
#' @param sigma_pp Variance in p.
#' @param sigma_qp Symmetric covariance (default 0).
#' @param hbar Planck constant in chosen units.
#' @return Named list with kinematic geometry and diagnostic flags.
robertson_schroedinger <- function(sigma_qq, sigma_pp, sigma_qp=0, hbar=1.0) {
  rs_lhs       <- sigma_qq*sigma_pp - sigma_qp^2
  rs_bound     <- (hbar/2)^2
  rs_satisfied <- rs_lhs >= rs_bound * (1 - RS_TOLERANCE)
  if (!rs_satisfied) warning(sprintf(
    "RS inequality violated: %.6e < %.6e", rs_lhs, rs_bound))

  Delta_q      <- sqrt(2*sigma_qq)
  Delta_p      <- sqrt(2*sigma_pp)
  delta_q      <- hbar/Delta_p
  delta_p      <- hbar/Delta_q

  sp_product   <- delta_q*Delta_p
  sp_satisfied <- abs(sp_product-hbar) < RS_TOLERANCE*hbar
  if (!sp_satisfied) warning(sprintf(
    "SP not saturated: %.6e != %.6e", sp_product, hbar))

  list(
    Delta_q      = Delta_q,
    Delta_p      = Delta_p,
    delta_q      = delta_q,
    delta_p      = delta_p,
    A_over_A0    = (Delta_q*Delta_p)/hbar,
    rs_satisfied = rs_satisfied,
    sp_satisfied = sp_satisfied
  )
}

# ------------------------------------------------------------------------------
# NUMERICAL COVARIANCE FROM SAMPLED WAVEFUNCTION
# ------------------------------------------------------------------------------

#' RS covariance from a numerically sampled wavefunction.
#'
#' Computes sigma_qq and sigma_pp from psi(q) on a uniform q grid, then
#' returns the kinematic geometry via robertson_schroedinger().
#'
#' sigma_pp is computed from <p^2> = hbar^2 * integral |dpsi/dq|^2 dq,
#' valid for normalizable psi where boundary terms vanish. Finite-difference
#' derivative; for high quantum numbers use a fine grid to avoid Nyquist
#' aliasing of psi's oscillations.
#'
#' @param psi_vec Wavefunction sampled on q_grid (will be re-normalized).
#' @param q_grid Uniform position grid.
#' @param hbar Planck constant in chosen units.
#' @return Same structure as robertson_schroedinger().
numerical_covariance <- function(psi_vec, q_grid, hbar=1.0) {
  dq       <- diff(q_grid)[1]
  norm     <- sqrt(sum(psi_vec^2)*dq)
  psi_n    <- psi_vec/norm
  prob     <- psi_n^2
  q_mean   <- sum(q_grid*prob)*dq
  sigma_qq <- sum((q_grid-q_mean)^2*prob)*dq
  dpsi     <- diff(psi_n)/dq
  sigma_pp <- hbar^2*sum(dpsi^2)*dq
  robertson_schroedinger(sigma_qq, sigma_pp, sigma_qp=0, hbar=hbar)
}
