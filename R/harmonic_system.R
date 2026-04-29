# ==============================================================================
# harmonic_system.R
# Quantum harmonic oscillator: V(q) = q^2/2
#
# Eigenstates are analytically known (Hermite-Gaussian functions), so no
# Schroedinger solver is needed. We expose a faux "soln" object with the
# same interface as solve_schroedinger() output (energies, psi_matrix,
# q_grid, dq) so the existing Wigner pipeline can be reused unchanged.
#
# Wavefunction:
#   psi_n(q) = (1/sqrt(2^n n! sqrt(pi))) * H_n(q) * exp(-q^2/2)
# Energies:
#   E_n = n + 1/2  (in units where hbar = m = omega = 1)
# Bohr-Sommerfeld action:
#   A_BS/A_0 = 2n+1  (exact)
#
# Reference: Griffiths Introduction to Quantum Mechanics Ch.2
#            Wigner 1932 Phys. Rev. 40, 749
# Author: Brian S. Mulloy
# ==============================================================================

# ------------------------------------------------------------------------------
# HERMITE POLYNOMIALS (physicists' convention)
# H_0 = 1, H_1 = 2q, H_{n+1} = 2q*H_n - 2n*H_{n-1}
# ------------------------------------------------------------------------------

#' Evaluate the n-th Hermite polynomial at q (vectorized in q).
hermite_poly <- function(n, q) {
  if (n == 0) return(rep(1, length(q)))
  if (n == 1) return(2*q)
  H_prev <- rep(1, length(q))
  H_curr <- 2*q
  for (k in 1:(n-1)) {
    H_next <- 2*q*H_curr - 2*k*H_prev
    H_prev <- H_curr
    H_curr <- H_next
  }
  H_curr
}

# ------------------------------------------------------------------------------
# HERMITE-GAUSSIAN WAVEFUNCTION
# Uses log-space normalization to avoid n! overflow at high n.
# log(2^n * n! * sqrt(pi)) = n*log(2) + lgamma(n+1) + (1/2)*log(pi)
# ------------------------------------------------------------------------------

#' Normalized harmonic oscillator eigenstate psi_n(q).
#' Uses log-space to handle high n (e.g. n=100) without factorial overflow.
harmonic_psi <- function(n, q) {
  log_norm <- 0.5*(n*log(2) + lgamma(n+1) + 0.5*log(pi))
  H_n      <- hermite_poly(n, q)
  exp(-q^2/2 - log_norm) * H_n
}

# ------------------------------------------------------------------------------
# HARMONIC POTENTIAL AND TURNING POINTS
# ------------------------------------------------------------------------------

harmonic_V <- function(q) 0.5 * q^2

#' Classical turning points at energy E for V = q^2/2: q_+- = +-sqrt(2E)
harmonic_turning_points <- function(E_n) {
  qt <- sqrt(2*E_n)
  list(q_minus = -qt, q_plus = qt)
}

# ------------------------------------------------------------------------------
# FAUX SOLN OBJECT
# Builds the same structure that solve_schroedinger() returns, but populated
# from analytic formulas. Compatible with build_*_row consumers.
# ------------------------------------------------------------------------------

#' Build a soln object containing the first n_states harmonic eigenstates
#' sampled on a uniform q grid. Same interface as solve_schroedinger output.
harmonic_soln <- function(n_states, q_min=-25, q_max=25, dq=0.02) {
  q_grid <- seq(q_min, q_max, by=dq)
  nq     <- length(q_grid)
  energies   <- (0:(n_states-1)) + 0.5
  psi_matrix <- matrix(0, nrow=nq, ncol=n_states)

  for (j in seq_len(n_states)) {
    n_val <- j - 1
    psi   <- harmonic_psi(n_val, q_grid)
    norm  <- sqrt(sum(psi^2)*dq)
    if (norm > 0) psi <- psi/norm
    psi_matrix[, j] <- psi
  }

  cat(sprintf("  Harmonic analytic: %d grid points, %d states\n", nq, n_states))
  for (j in seq_len(n_states)) {
    cat(sprintf("    n=%d: E=%.6f norm=%.6f\n",
                j-1, energies[j], sum(psi_matrix[,j]^2)*dq))
  }

  list(energies=energies, psi_matrix=psi_matrix, q_grid=q_grid, dq=dq)
}

# ------------------------------------------------------------------------------
# Grid parameters used by the plot files.
# n=100 has orbit radius sqrt(201) ~ 14.2; grid extends well beyond.
# dq=0.02 gives ~50 samples per node spacing for n=100 (oversampled).
# ------------------------------------------------------------------------------
HARMONIC_Q_MIN    <- -25.0
HARMONIC_Q_MAX    <-  25.0
HARMONIC_DQ       <-   0.02
HARMONIC_N_STATES <- 101    # need n=100 to be available

# ------------------------------------------------------------------------------
# ENERGY FROM ACTION CAPACITY
#
# For the harmonic oscillator the orbit covariance is exact:
#   Delta_q = sqrt(2E),  Delta_p = sqrt(2E)  =>  A/A_0 = 2E
# Inverse: E = (A/A_0) / 2.
#
# Equivalently in eigenstate terms: A_RS/A_0 = 2n+1, so A/A_0 = 1
# corresponds to n=0 (the ground state), A/A_0 = 3 to n=1, A/A_0 = 9 to
# n=4, etc. Continuous A/A_0 between integer levels gives a classical
# orbit at the corresponding energy with no quantum eigenstate associated.
# ------------------------------------------------------------------------------

#' Energy at a given action-capacity level, in dimensionless units.
#'
#' @param A_over_A0 Target action capacity in units of A_0 = h/2 (positive real).
#' @return Energy E such that orbit_covariance(harmonic_V, E, tp) yields
#'         A/A_0 = A_over_A0.
harmonic_E_at_action_capacity <- function(A_over_A0) {
  if (A_over_A0 <= 0) stop("A/A_0 must be positive")
  A_over_A0 / 2
}
