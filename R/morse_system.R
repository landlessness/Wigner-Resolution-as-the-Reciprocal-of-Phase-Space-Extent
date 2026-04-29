# ==============================================================================
# morse_system.R
# Morse potential parameters, energy levels, and turning points.
#
# V(q) = D_e*(1-exp(-alpha*q))^2
#
# Dimensionless parameters chosen to span the textbook Morse phase-space
# regime: from compact harmonic-like ground state to dramatic horseshoe
# orbits near dissociation.
#
# Reference: Morse 1929 Phys. Rev. 34, 57
#            Dahl & Springborg J.Chem.Phys. 88, 4535 (1988)
#            Le Roy, J. Quant. Spectrosc. Radiat. Transfer 186, 158 (2017)
#            (RKR1 program: semiclassical inversion of diatomic spectra)
#
# Bound states: n = 0..N_max with N_max = floor(sqrt(2*D_e)/alpha - 0.5).
# For D_e=12.5, alpha=0.2: N_max = floor(24.5) = 24
#
# Schroedinger grid: q in [-2, 50], dq=0.01
# Author: Brian S. Mulloy
# ==============================================================================

De    <- 12.5
alpha <- 0.2

morse_V <- function(q) De * (1 - exp(-alpha*q))^2

morse_turning_points <- function(E_n) {
  list(q_minus = -log(1 + sqrt(E_n/De)) / alpha,
       q_plus  = -log(1 - sqrt(E_n/De)) / alpha)
}

# Schroedinger solver grid parameters
MORSE_Q_MIN    <- -2.0
MORSE_Q_MAX    <- 50.0
MORSE_DQ       <-  0.01
MORSE_N_STATES <- 18

# ------------------------------------------------------------------------------
# ENERGY FROM ACTION CAPACITY
#
# For Morse, the orbit covariance is asymmetric and the Fermi blob A
# differs from the Bohr-Sommerfeld classical action A_cl by O(1/n)
# corrections. We invert numerically: bisect on E in [0, D_e] to find
# the energy at which orbit_covariance(morse_V, E, tp) returns the
# requested A/A_0.
#
# Requires classical_action_tools.R to be sourced for orbit_covariance.
# ------------------------------------------------------------------------------

#' Energy at a given Fermi-blob action-capacity level for Morse.
#'
#' Numerically inverts A(E) = pi * Delta_q^orbit(E) * Delta_p^orbit(E) for E.
#' Brackets E in (E_lo, E_hi) where E_lo is just above zero (so the orbit
#' exists) and E_hi is just below D_e (the dissociation limit).
#'
#' @param A_over_A0 Target action capacity in units of A_0 = h/2.
#' @param tol Bisection tolerance on |A_obtained - A_target|.
#' @return Energy E in (0, D_e).
morse_E_at_action_capacity <- function(A_over_A0, tol=1e-4) {
  if (A_over_A0 <= 0) stop("A/A_0 must be positive")

  A_at <- function(E) {
    tp  <- morse_turning_points(E)
    cov <- suppressWarnings(orbit_covariance(morse_V, E, tp))
    cov$A_over_A0
  }

  E_lo <- 1e-6 * De
  E_hi <- (1 - 1e-6) * De
  A_lo <- A_at(E_lo)
  A_hi <- A_at(E_hi)
  if (A_over_A0 < A_lo || A_over_A0 > A_hi) {
    stop(sprintf(paste("A/A_0 = %.3f outside Morse bound-state range",
                       "[%.3f, %.3f]"), A_over_A0, A_lo, A_hi))
  }

  for (iter in 1:100) {
    E_mid <- 0.5 * (E_lo + E_hi)
    A_mid <- A_at(E_mid)
    if (abs(A_mid - A_over_A0) < tol) return(E_mid)
    if (A_mid < A_over_A0) E_lo <- E_mid else E_hi <- E_mid
  }
  warning(sprintf("morse_E_at_action_capacity: bisection did not converge to %.1e",
                  tol))
  0.5 * (E_lo + E_hi)
}
