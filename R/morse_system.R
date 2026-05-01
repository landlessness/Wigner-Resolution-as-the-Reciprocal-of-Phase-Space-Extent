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

# Harmonic frequency at the well bottom (where V'' = 2*De*alpha^2 = m*omega^2,
# m = 1). For De=12.5, alpha=0.2 this gives omega = alpha*sqrt(2*De) = 1.0.
morse_omega <- alpha * sqrt(2 * De)

# ------------------------------------------------------------------------------
# BOHR-SOMMERFELD ENERGY LEVELS (analytic Morse)
#
# E_n = omega*(n + 1/2) - [omega*(n + 1/2)]^2 / (4*De)
#
# This is the exact Morse spectrum (and also the WKB / Bohr-Sommerfeld
# result, which is exact for Morse — a peculiarity of the potential).
# Used by indexing rows by quantum number rather than by action capacity.
#
# Reference: Morse 1929 Phys. Rev. 34, 57 (Eq. 23).
# ------------------------------------------------------------------------------

#' Bohr-Sommerfeld (= exact) Morse energy at quantum number n.
#'
#' @param n Non-negative integer; must satisfy n <= floor(sqrt(2*De)/alpha-0.5).
#' @return Energy E_n.
morse_E_BS <- function(n) {
  N_max <- floor(sqrt(2*De)/alpha - 0.5)
  if (n < 0 || n > N_max)
    stop(sprintf("n=%d outside Morse bound range [0, %d] for D_e=%.3f, alpha=%.3f",
                 n, N_max, De, alpha))
  morse_omega*(n + 0.5) - (morse_omega*(n + 0.5))^2 / (4*De)
}

# ------------------------------------------------------------------------------
# ENERGY FROM ACTION CAPACITY
#
# For Morse, the orbit covariance is asymmetric and the Fermi blob A
# differs from the Bohr-Sommerfeld classical action A_cl by O(1/n)
# corrections. We invert numerically: bisect on E in [0, D_e] to find
# the energy at which orbit_covariance(morse_V, E, tp) returns the
# requested A/A_0.
#
# Note: the Morse A_orbit(E) curve is non-monotone; see the implementation
# of morse_E_at_action_capacity() below for how the rising branch is
# isolated.
#
# Requires classical_action_tools.R to be sourced for orbit_covariance.
# ------------------------------------------------------------------------------

#' Energy at a given Fermi-blob action-capacity level for Morse.
#'
#' Numerically inverts A(E) = pi * Delta_q^orbit(E) * Delta_p^orbit(E) for E.
#' The Morse A(E) curve is non-monotone: it rises from zero through the bound
#' spectrum, peaks at some E_peak < D_e, then falls back toward zero as the
#' orbit becomes infinitely extended in q with vanishingly small typical |p|.
#' This routine therefore (a) locates E_peak by golden-section search, (b)
#' bisects on the strictly-rising branch (E_lo, E_peak), and (c) errors out
#' with a useful message if A_over_A0 exceeds the achievable peak.
#'
#' @param A_over_A0 Target action capacity in units of A_0 = h/2.
#' @param tol Bisection tolerance on |A_obtained - A_target|.
#' @return Energy E in (0, E_peak).
morse_E_at_action_capacity <- function(A_over_A0, tol=1e-4) {
  if (A_over_A0 <= 0) stop("A/A_0 must be positive")

  A_at <- function(E) {
    tp  <- morse_turning_points(E)
    cov <- suppressWarnings(orbit_covariance(morse_V, E, tp))
    cov$A_over_A0
  }

  E_lo <- 1e-6 * De
  E_hi <- (1 - 1e-6) * De

  # Locate the peak of A(E) by golden-section search on (E_lo, E_hi). The
  # function is unimodal in this interval (monotone rise to a single maximum,
  # then monotone fall), so golden section is the textbook tool.
  phi <- (sqrt(5) - 1) / 2
  a <- E_lo; b <- E_hi
  c_pt <- b - phi*(b - a); d_pt <- a + phi*(b - a)
  for (iter in 1:80) {
    if (abs(b - a) < 1e-5) break
    if (A_at(c_pt) < A_at(d_pt)) {
      a <- c_pt
    } else {
      b <- d_pt
    }
    c_pt <- b - phi*(b - a); d_pt <- a + phi*(b - a)
  }
  E_peak <- 0.5*(a + b)
  A_peak <- A_at(E_peak)

  A_at_lo <- A_at(E_lo)
  if (A_over_A0 < A_at_lo) {
    stop(sprintf("A/A_0 = %.3f below the lower-energy floor A/A_0 = %.3f",
                 A_over_A0, A_at_lo))
  }
  if (A_over_A0 > A_peak) {
    stop(sprintf(paste("A/A_0 = %.3f exceeds the maximum action capacity",
                       "A_peak/A_0 = %.3f reached at E = %.3f for this Morse",
                       "potential (D_e = %.3f, alpha = %.3f). The Morse A(E)",
                       "curve is non-monotone and bounded above; choose a",
                       "smaller target or re-parameterize the potential."),
                 A_over_A0, A_peak, E_peak, De, alpha))
  }

  # Bisect on the rising branch [E_lo, E_peak], which is strictly monotone.
  bracket_lo <- E_lo
  bracket_hi <- E_peak
  for (iter in 1:100) {
    E_mid <- 0.5 * (bracket_lo + bracket_hi)
    A_mid <- A_at(E_mid)
    if (abs(A_mid - A_over_A0) < tol) return(E_mid)
    if (A_mid < A_over_A0) bracket_lo <- E_mid else bracket_hi <- E_mid
  }
  warning(sprintf("morse_E_at_action_capacity: bisection did not converge to %.1e",
                  tol))
  0.5 * (bracket_lo + bracket_hi)
}
