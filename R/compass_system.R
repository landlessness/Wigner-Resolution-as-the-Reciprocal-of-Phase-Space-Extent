# ==============================================================================
# compass_system.R
# Zurek's compass state: a free-space coherent superposition with no
# governing potential.
#
# psi_compass(q) is a normalized superposition of four coherent states
# placed at the cardinal compass points in phase space:
#   |N> centered at (q, p) = (0,  L/2)
#   |S> centered at (q, p) = (0, -L/2)
#   |E> centered at (q, p) = (L/2, 0)
#   |W> centered at (q, p) = (-L/2, 0)
#
# Each lobe is a Gaussian coherent state with width xi in q. For real-valued
# psi(q), the |N> + |S> superposition produces the cos(L q / 2) factor:
#
#   psi(q) = N * [ exp(-(q - L/2)^2 / (2 xi^2))   |E> lobe
#                + exp(-(q + L/2)^2 / (2 xi^2))   |W> lobe
#                + 2 cos(L q / 2) * exp(-q^2 / (2 xi^2)) ]   |N>+|S> superposition
#
# where the |N> + |S> term comes from
#   exp(+i (L/2) q) exp(-q^2/(2 xi^2)) + exp(-i (L/2) q) exp(-q^2/(2 xi^2))
#   = 2 cos(L q / 2) exp(-q^2 / (2 xi^2))
#
# This file contains NO Hamiltonian. The compass state is a free-space
# coherent superposition. Its action capacity comes from the Robertson-
# Schroedinger covariance of psi (via numerical_covariance in
# quantum_tools.R), not from any classical orbit.
#
# The default parameters match Zurek (2001) Fig. 2: L = 10, xi = 1, hbar = 1.
# In these units the chessboard tile half-width is pi / L ~= 0.31, the
# coherent-state width is 1, and the Fermi-blob action capacity is
# A/A_0 ~= 25 with reciprocal scales delta_q = delta_p ~= 0.19.
#
# Reference: Zurek 2001 Nature 412, 712 (compass state and sub-Planck scales)
#            Praxmeyer et al. 2016 PRA 93, 053835 (optical realization)
#            Toscano et al. 2006 PRA 73, 023803 (compass-state metrology)
# Author: Brian S. Mulloy
# ==============================================================================

# ------------------------------------------------------------------------------
# DEFAULT PARAMETERS (Zurek 2001 Fig. 2)
# ------------------------------------------------------------------------------

COMPASS_L   <- 10.0   # half-distance from origin to each lobe
COMPASS_XI  <-  1.0   # coherent-state width

# ------------------------------------------------------------------------------
# COMPASS WAVEFUNCTION
# ------------------------------------------------------------------------------

#' Real-valued compass-state wavefunction psi(q) sampled on a grid.
#'
#' Returns the normalized superposition of |W>, |E>, and (|N> + |S>) lobes.
#' The |N> + |S> contribution collapses to a real cosine modulation by
#' parity (cos(L q / 2) * Gaussian centered at 0), which keeps the whole
#' psi real-valued and compatible with the wigner_fft pipeline.
#'
#' @param q     Position grid (numeric vector).
#' @param L     Half-distance from origin to each lobe (default COMPASS_L).
#' @param xi    Coherent-state width (default COMPASS_XI).
#' @return Numeric vector psi(q), L2-normalized on q.
compass_psi <- function(q, L=COMPASS_L, xi=COMPASS_XI) {
  lobe_E <- exp(-(q - L/2)^2 / (2 * xi^2))
  lobe_W <- exp(-(q + L/2)^2 / (2 * xi^2))
  lobe_N_plus_S <- 2 * cos(L * q / 2) * exp(-q^2 / (2 * xi^2))

  psi   <- lobe_W + lobe_E + lobe_N_plus_S
  dq    <- diff(q)[1]
  norm  <- sqrt(sum(psi^2) * dq)
  if (norm > 0) psi <- psi / norm
  psi
}

# ------------------------------------------------------------------------------
# DISPLAY GRID PARAMETERS
# Span enough of phase space to enclose all four lobes plus the chessboard
# region between them. With L = 10 the lobes sit at +/- 5 with width xi = 1,
# so a window of [-10, 10] gives a full lobe-and-a-half of padding.
# ------------------------------------------------------------------------------

COMPASS_Q_MIN    <- -10.0
COMPASS_Q_MAX    <-  10.0
COMPASS_DQ       <-   0.02
