# ==============================================================================
# double_well_potential.R
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
# Schrodinger grid: q in [-6, 6], dq=0.01
#   V(-6) = V(6) = 216 >> all bound state energies ✓
# ==============================================================================

mu2    <- 4.0
lambda <- 1.0

double_well_V <- function(q) -mu2/2 * q^2 + lambda/4 * q^4

double_well_barrier <- mu2^2 / (4*lambda)       # = 4.0
double_well_minima  <- sqrt(mu2/lambda)          # = 2.0
double_well_omega   <- sqrt(2) * sqrt(mu2)       # harmonic freq at minima

double_well_turning_points <- function(E_n) {
  # Solve V(q) = E_n numerically for turning points
  # Returns all real roots in ascending order
  roots <- polyroot(c(-E_n, 0, -mu2/2, 0, lambda/4))
  real_roots <- Re(roots[abs(Im(roots)) < 1e-8])
  real_roots <- sort(real_roots)
  # For below-barrier states: 4 turning points
  # For above-barrier states: 2 turning points
  list(roots=real_roots)
}

# Schrodinger solver grid parameters
DOUBLE_WELL_Q_MIN    <- -6.0
DOUBLE_WELL_Q_MAX    <-  6.0
DOUBLE_WELL_DQ       <-  0.01
DOUBLE_WELL_N_STATES <-  8
