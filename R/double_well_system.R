# ==============================================================================
# double_well_system.R
# Double-well potential parameters, barrier, minima, and turning points.
#
# V(q) = -mu2/2 * q^2 + lambda/4 * q^4
#
# Standard dimensionless parameters:
#   mu2=4, lambda=1
#   Barrier height: V_0 = mu2^2/(4*lambda) = 4
#   Well minima at: q = +/- sqrt(mu2/lambda) = +/- 2
#   Harmonic frequency at minima: omega_well = sqrt(2)*mu
#
# Reference: Griffiths Introduction to Quantum Mechanics Problem 2.47
#            Razavy, Quantum Theory of Tunneling (World Scientific 2003) Ch.2
#
# Schroedinger grid: q in [-6, 6], dq=0.01
#   V(-6) = V(6) = 216 >> all bound state energies
# Author: Brian S. Mulloy
# ==============================================================================

mu2    <- 4.0
lambda <- 1.0

double_well_V <- function(q) -mu2/2 * q^2 + lambda/4 * q^4

double_well_barrier <- mu2^2 / (4*lambda)       # = 4.0
double_well_minima  <- sqrt(mu2/lambda)          # = 2.0
double_well_omega   <- sqrt(2) * sqrt(mu2)       # harmonic freq at minima

double_well_turning_points <- function(E_n) {
  # Solve V(q) = E_n numerically for turning points.
  # Returns all real roots in ascending order.
  # For below-barrier states: 4 turning points
  # For above-barrier states: 2 turning points
  roots <- polyroot(c(-E_n, 0, -mu2/2, 0, lambda/4))
  real_roots <- Re(roots[abs(Im(roots)) < 1e-8])
  real_roots <- sort(real_roots)
  list(roots=real_roots)
}

# Schroedinger solver grid parameters
DOUBLE_WELL_Q_MIN    <- -6.0
DOUBLE_WELL_Q_MAX    <-  6.0
DOUBLE_WELL_DQ       <-  0.01
DOUBLE_WELL_N_STATES <-  8

# ------------------------------------------------------------------------------
# ENERGY FROM ACTION CAPACITY
#
# Numerical inversion of A(E) = pi * Delta_q^orbit(E) * Delta_p^orbit(E)
# for the double-well. Below the barrier (E < V_barrier = 4), the
# classical orbit consists of two disconnected segments and orbit_covariance
# integrates over both; above the barrier, single connected orbit. The
# combined moments give a continuous A(E) across the barrier-crossing
# regime change, monotonically increasing in E.
#
# Bracket: E in (V_min_strictly_above_minimum, E_max). The well minima
# are at V = -mu2^2 / (4*lambda) = -4 (for the standard parameters);
# E_max is set well above the barrier to cover above-barrier states.
#
# Requires classical_action_tools.R to be sourced for orbit_covariance.
# ------------------------------------------------------------------------------

#' Energy at a given Fermi-blob action-capacity level for the double well.
#'
#' Numerically inverts A(E) = pi * Delta_q^orbit(E) * Delta_p^orbit(E) for
#' E. Handles both sub-barrier (4 turning points) and above-barrier (2
#' turning points) regimes within one continuous monotone inversion.
#'
#' @param A_over_A0 Target action capacity in units of A_0 = h/2.
#' @param tol Bisection tolerance on |A_obtained - A_target|.
#' @return Energy E.
double_well_E_at_action_capacity <- function(A_over_A0, tol=1e-4) {
  if (A_over_A0 <= 0) stop("A/A_0 must be positive")

  V_min <- -mu2^2 / (4 * lambda)  # = -4

  A_at <- function(E) {
    tp_obj <- double_well_turning_points(E)
    cov    <- suppressWarnings(orbit_covariance(double_well_V, E, tp_obj$roots))
    cov$A_over_A0
  }

  # Search range: just above the well bottom to a value well above the barrier.
  E_lo <- V_min + 1e-3
  E_hi <- 50.0
  A_lo <- A_at(E_lo)
  A_hi <- A_at(E_hi)
  if (A_over_A0 < A_lo || A_over_A0 > A_hi) {
    stop(sprintf(paste("A/A_0 = %.3f outside double-well range",
                       "[%.3f, %.3f] for E in [%.3f, %.3f]"),
                 A_over_A0, A_lo, A_hi, E_lo, E_hi))
  }

  for (iter in 1:100) {
    E_mid <- 0.5 * (E_lo + E_hi)
    A_mid <- A_at(E_mid)
    if (abs(A_mid - A_over_A0) < tol) return(E_mid)
    if (A_mid < A_over_A0) E_lo <- E_mid else E_hi <- E_mid
  }
  warning(sprintf("double_well_E_at_action_capacity: bisection did not converge to %.1e",
                  tol))
  0.5 * (E_lo + E_hi)
}
