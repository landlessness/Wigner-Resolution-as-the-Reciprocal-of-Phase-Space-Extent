# ==============================================================================
# table_marginals.R
# Driver: loop over the eight states, compute marginals (exact and convolved),
# and emit a table of residuals.
#
# Eight states (matching plot_wigner.R and plot_cats.R rows):
#   1. squeezed vacuum (r = 0.5)
#   2. harmonic Fock n = 1
#   3. Morse n = 8
#   4. double-well n = 5
#   5. 2-cat
#   6. 3-cat
#   7. 4-cat square
#   8. 4-cat Zurek compass
#
# For each state:
#   - sample psi
#   - compute Delta_q, Delta_p from psi (numerical_covariance)
#   - build Wigner state on integration grid
#   - compute symplectic_q_marginal and exact_position_marginal on q_display
#   - compute symplectic_p_marginal and exact_momentum_marginal on p_display
#   - compute L1, L2, Linf residuals for both directions
#
# Output:
#   - per-state list of (q grids, p grids, all four marginals, residuals)
#     saved to data/marginals_data.rds
#   - residual table written to data/marginals_residuals.csv
#   - residual summary printed to console
#
# Author: Brian S. Mulloy
# ==============================================================================

library(here)
library(data.table)

source(here("R", "math_tools.R"))
source(here("R", "quantum_tools.R"))         # numerical_covariance
source(here("R", "wigner_density.R"))         # build_wigner_state
source(here("R", "symplectic_kernel.R"))      # G_delta_q_kernel_matrix
source(here("R", "marginal_tools.R"))         # the new helpers
source(here("R", "schroedinger_solver.R"))    # solve_schroedinger

source(here("R", "harmonic_system.R"))
source(here("R", "morse_system.R"))
source(here("R", "double_well_system.R"))
source(here("R", "cat_system.R"))

# ------------------------------------------------------------------------------
# STATE DESCRIPTORS
#
# Each descriptor returns:
#   $name        short label used in output
#   $psi_fn      function() -> list(psi_vec, psi_q_grid)
#   $window      function(Delta_q, Delta_p) -> list(q_lo, q_hi, p_lo, p_hi)
#                Defaults provided where states have a natural scale; may be
#                overridden per state to widen the integration domain.
# ------------------------------------------------------------------------------

# Default windows: 4 Delta in each direction. Wide enough to capture support
# without truncating the Wigner integrand.
default_window <- function(Delta_q, Delta_p) {
  list(q_lo = -4*Delta_q, q_hi = 4*Delta_q,
       p_lo = -4*Delta_p, p_hi = 4*Delta_p)
}

# ---- Squeezed vacuum (analytical Gaussian) ----
# psi(q) = (pi sigma_q^2)^{-1/4} exp(-q^2 / (2 sigma_q^2))
# with sigma_q = exp(-r), so Delta_q = exp(-r), Delta_p = exp(+r) (hbar=1).
SQUEEZED_R <- 0.5

squeezed_vacuum_descriptor <- list(
  name = "squeezed_vacuum",
  psi_fn = function() {
    psi_q_grid <- seq(-15, 15, by=0.005)
    sigma_q    <- exp(-SQUEEZED_R)
    psi_vec    <- (pi * sigma_q^2)^(-1/4) *
      exp(-psi_q_grid^2 / (2 * sigma_q^2))
    list(psi_vec=psi_vec, psi_q_grid=psi_q_grid)
  },
  window = default_window
)

# ---- Harmonic Fock n=1 ----
harmonic_n1_descriptor <- list(
  name = "harmonic_n1",
  psi_fn = function() {
    psi_q_grid <- seq(-25, 25, by=0.02)
    psi_vec    <- harmonic_psi(1, psi_q_grid)
    list(psi_vec=psi_vec, psi_q_grid=psi_q_grid)
  },
  window = default_window
)

# ---- Morse n=8 ----
# Solve once; capture into closure so repeated calls reuse the solution.
.morse_soln <- NULL
ensure_morse <- function() {
  if (is.null(.morse_soln)) {
    cat("  Solving Schroedinger for Morse...\n")
    .morse_soln <<- suppressMessages(solve_schroedinger(
      morse_V, MORSE_Q_MIN, MORSE_Q_MAX, MORSE_DQ,
      n_states=MORSE_N_STATES))
  }
  .morse_soln
}

morse_n8_descriptor <- list(
  name = "morse_n8",
  psi_fn = function() {
    soln       <- ensure_morse()
    psi_q_grid <- soln$q_grid
    psi_vec    <- soln$psi_matrix[, 8 + 1]   # n=8 is column 9
    list(psi_vec=psi_vec, psi_q_grid=psi_q_grid)
  },
  window = function(Delta_q, Delta_p) {
    # Morse is asymmetric and Delta_q underestimates the tail; widen.
    list(q_lo = -6*Delta_q, q_hi = 6*Delta_q,
         p_lo = -4*Delta_p, p_hi = 4*Delta_p)
  }
)

# ---- Double-well n=5 ----
.dw_soln <- NULL
ensure_dw <- function() {
  if (is.null(.dw_soln)) {
    cat("  Solving Schroedinger for double well...\n")
    .dw_soln <<- suppressMessages(solve_schroedinger(
      double_well_V, DOUBLE_WELL_Q_MIN, DOUBLE_WELL_Q_MAX,
      DOUBLE_WELL_DQ, n_states=DOUBLE_WELL_N_STATES))
  }
  .dw_soln
}

double_well_n5_descriptor <- list(
  name = "double_well_n5",
  psi_fn = function() {
    soln       <- ensure_dw()
    psi_q_grid <- soln$q_grid
    psi_vec    <- soln$psi_matrix[, 5 + 1]
    list(psi_vec=psi_vec, psi_q_grid=psi_q_grid)
  },
  window = default_window
)

# ---- Cat states ----
# Use cat_psi() directly. Naming conventions match plot_cats.R.
make_cat_descriptor <- function(name, n_cats, variant="diag") {
  list(
    name = name,
    psi_fn = function() {
      psi_q_grid <- seq(CAT_Q_MIN, CAT_Q_MAX, by=CAT_DQ)
      psi_vec    <- cat_psi(psi_q_grid, n_cats, variant=variant)
      list(psi_vec=psi_vec, psi_q_grid=psi_q_grid)
    },
    window = function(Delta_q, Delta_p) {
      # Cats: support extends to the lobe centers; widen to 5 Delta.
      list(q_lo = -5*Delta_q, q_hi = 5*Delta_q,
           p_lo = -5*Delta_p, p_hi = 5*Delta_p)
    }
  )
}

# 2-cat, 3-cat, 4-cat square (diag), Zurek compass (axis variant).
cat2_descriptor    <- make_cat_descriptor("cat_2",        n_cats=2)
cat3_descriptor    <- make_cat_descriptor("cat_3",        n_cats=3)
cat4sq_descriptor  <- make_cat_descriptor("cat_4_square", n_cats=4, variant="diag")
compass_descriptor <- make_cat_descriptor("cat_compass",  n_cats=4, variant="axis")

ALL_DESCRIPTORS <- list(
  squeezed_vacuum_descriptor,
  harmonic_n1_descriptor,
  morse_n8_descriptor,
  double_well_n5_descriptor,
  cat2_descriptor,
  cat3_descriptor,
  cat4sq_descriptor,
  compass_descriptor
)

# ------------------------------------------------------------------------------
# PER-STATE PROCESSING
#
# Returns a list with all marginals and residuals for one state.
# ------------------------------------------------------------------------------

process_state <- function(descriptor, n_q_display=1001, n_p_display=1001) {
  cat(sprintf("\n== %s ==\n", descriptor$name))

  # 1. Sample psi.
  ps         <- descriptor$psi_fn()
  psi_vec    <- ps$psi_vec
  psi_q_grid <- ps$psi_q_grid

  # 2. RS covariance from psi.
  rs <- numerical_covariance(psi_vec, psi_q_grid, hbar=1.0)
  cat(sprintf("  Delta_q=%.4f Delta_p=%.4f delta_q=%.4f delta_p=%.4f\n",
              rs$Delta_q, rs$Delta_p, rs$delta_q, rs$delta_p))

  # 3. Display windows (also used as integration extents).
  w  <- descriptor$window(rs$Delta_q, rs$Delta_p)
  q_lo <- w$q_lo; q_hi <- w$q_hi
  p_lo <- w$p_lo; p_hi <- w$p_hi

  q_display <- seq(q_lo, q_hi, length.out=n_q_display)
  p_display <- seq(p_lo, p_hi, length.out=n_p_display)
  dq_display <- diff(q_display)[1]
  dp_display <- diff(p_display)[1]

  # 4. Wigner state on the integration grid.
  cat("  Building Wigner state...\n")
  state <- build_wigner_state(
    psi_vec=psi_vec, psi_q_grid=psi_q_grid,
    q_lo=q_lo, q_hi=q_hi, p_lo=p_lo, p_hi=p_hi,
    q_display=q_display)

  # 5. Convolved marginals.
  cat("  Computing P_{delta q} marginal over p...\n")
  pdq_q_marg <- symplectic_q_marginal(state, rs$Delta_q, rs$Delta_p,
                                      q_display, hbar=1.0)
  cat("  Computing P_{delta p} marginal over q...\n")
  pdp_p_marg <- symplectic_p_marginal(state, rs$Delta_q, rs$Delta_p,
                                      p_display, hbar=1.0)

  # 6. Exact marginals.
  cat("  Computing |psi(q)|^2 and |psi-hat(p)|^2...\n")
  exact_q <- exact_position_marginal(psi_vec, psi_q_grid, q_display)
  exact_p <- exact_momentum_marginal(psi_vec, psi_q_grid, p_display, hbar=1.0)

  # 7. Residuals.
  res_q <- compute_residuals(exact_q, pdq_q_marg, dq_display)
  res_p <- compute_residuals(exact_p, pdp_p_marg, dp_display)

  cat(sprintf("  Position marginal: L1=%.3f%% L2=%.3f%% Linf=%.3f%%\n",
              100*res_q$l1, 100*res_q$l2, 100*res_q$linf))
  cat(sprintf("  Momentum marginal: L1=%.3f%% L2=%.3f%% Linf=%.3f%%\n",
              100*res_p$l1, 100*res_p$l2, 100*res_p$linf))

  list(
    name       = descriptor$name,
    Delta_q    = rs$Delta_q,
    Delta_p    = rs$Delta_p,
    delta_q    = rs$delta_q,
    delta_p    = rs$delta_p,
    q_display  = q_display,
    p_display  = p_display,
    exact_q    = exact_q,
    exact_p    = exact_p,
    pdq_q_marg = pdq_q_marg,
    pdp_p_marg = pdp_p_marg,
    res_q      = res_q,
    res_p      = res_p
  )
}

# ------------------------------------------------------------------------------
# DRIVE
# ------------------------------------------------------------------------------

cat("Computing marginal residuals across eight states...\n")

results <- lapply(ALL_DESCRIPTORS, process_state)
names(results) <- sapply(ALL_DESCRIPTORS, function(d) d$name)

# Build the residual table.
residual_table <- data.table(
  state    = sapply(results, function(r) r$name),
  Delta_q  = sapply(results, function(r) r$Delta_q),
  Delta_p  = sapply(results, function(r) r$Delta_p),
  delta_q  = sapply(results, function(r) r$delta_q),
  delta_p  = sapply(results, function(r) r$delta_p),
  q_L1_pct   = 100*sapply(results, function(r) r$res_q$l1),
  q_L2_pct   = 100*sapply(results, function(r) r$res_q$l2),
  q_Linf_pct = 100*sapply(results, function(r) r$res_q$linf),
  p_L1_pct   = 100*sapply(results, function(r) r$res_p$l1),
  p_L2_pct   = 100*sapply(results, function(r) r$res_p$l2),
  p_Linf_pct = 100*sapply(results, function(r) r$res_p$linf)
)

cat("\n==== Residual summary ====\n")
print(residual_table)

# Persist for the figure renderer.
dir.create(here("data"), showWarnings=FALSE)
saveRDS(results, here("data", "marginals_data.rds"))
fwrite(residual_table, here("data", "marginals_residuals.csv"))

cat(sprintf("\nWrote %s\n", here("data", "marginals_data.rds")))
cat(sprintf("Wrote %s\n", here("data", "marginals_residuals.csv")))
