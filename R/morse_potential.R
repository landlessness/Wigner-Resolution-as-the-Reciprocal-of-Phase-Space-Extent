# ==============================================================================
# morse_potential.R
# Morse potential parameters, energy levels, and turning points.
#
# V(q) = D_e*(1-exp(-alpha*q))^2
#
# Dimensionless parameters (D_e=30, alpha=1) chosen to span the textbook
# Morse phase-space regime: from compact harmonic-like ground state to
# dramatic horseshoe orbits near dissociation.
#
# Reference: Morse 1929 Phys. Rev. 34, 57
#            Dahl & Springborg J.Chem.Phys. 88, 4535 (1988)
#
# Bound states: n = 0..7
# N_max = floor(sqrt(2*D_e)/alpha - 0.5) = floor(7.245) = 7
#
# Schrodinger grid: q in [-2, 25], dq=0.01
#   V(-2) = 30*(1-exp(2))^2 = 1234 >> E_7=29.97  ✓
#   V(25) ≈ 30 ≈ E_7=29.97 (n=7 has very long tail near dissociation)
#   q_max=25 captures n=0..6 fully and n=7 to ~99% of its support.
# ==============================================================================
De    <- 12.5
alpha <- 0.2
morse_V <- function(q) De * (1 - exp(-alpha*q))^2
morse_turning_points <- function(E_n) {
  list(q_minus = -log(1 + sqrt(E_n/De)) / alpha,
       q_plus  = -log(1 - sqrt(E_n/De)) / alpha)
}
# Schrodinger solver grid parameters
MORSE_Q_MIN    <- -2.0
MORSE_Q_MAX    <- 50
MORSE_DQ       <- 0.01
MORSE_N_STATES <- 18
