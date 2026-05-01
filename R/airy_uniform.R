# ==============================================================================
# airy_uniform.R
# Berry-Mount / Langer uniform Airy approximation for the WKB caustic.
#
# Canonical Griffiths-Ch.8 construction: the Langer (1937) uniform
# wavefunction
#
#   psi_Langer(q) = N * |p(q)|^{-1/2} * z(q)^{1/4} * Ai(-z(q))
#
# (with appropriate sign flip and Bi tail in the forbidden region) is
# the textbook fix to the WKB caustic. It approximates the n-th
# eigenstate's wavefunction across the entire orbit, oscillating n+1
# times between the turning points (because Ai(-z) for large z is
# oscillatory, with phase advance pi(n+1/2) over the orbit by Bohr-
# Sommerfeld), and finite at the turning points themselves (where
# |p|^{-1/2} diverges and z^{1/4} -> 0 in compensation).
#
# We return the modulus squared:
#
#   rho_Airy(q) = |psi_Langer(q)|^2 = |p(q)|^{-1} * z(q)^{1/2} * Ai^2(-z(q))
#
# in the classically allowed region, and Ai^2(z) (decaying tail) in
# the forbidden region. This is the canonical comparator to a
# wavefunction-squared density: oscillating with n+1 lobes inside the
# orbit, finite at turning points (replacing the WKB caustic
# divergence with an Airy lobe), exponentially decaying outside.
#
# Apples-to-apples comparator with the symplectic-resolved oscillating
# WKB density: both methods take the diverging WKB |psi_WKB|^2 input
# and produce a finite-at-turning-points oscillating density.
#
# Reference: Langer, Phys. Rev. 51, 669 (1937)
#            Berry & Mount, Rep. Prog. Phys. 35, 315 (1972), Sec. 4
#            Griffiths, Introduction to Quantum Mechanics, Ch. 8
# Author: Brian S. Mulloy
# ==============================================================================

library(here)

if (!requireNamespace("gsl", quietly=TRUE))
  stop("airy_uniform.R requires package 'gsl' (already a project dependency)")

# ------------------------------------------------------------------------------
# ACTION INTEGRAL FROM A TURNING POINT
# ------------------------------------------------------------------------------

#' Trapezoidal-rule action integral from q to a turning point q_t.
#'
#' Returns the absolute action |int_q^{q_t} sqrt(2|V-E|) dq'|, positive
#' in BOTH allowed and forbidden regions. The sign is folded into the
#' Airy argument by the caller.
.action_from_turning_point <- function(q_grid, q_t, V_fn, E, n_int=2001) {
  S <- numeric(length(q_grid))
  for (k in seq_along(q_grid)) {
    q <- q_grid[k]
    if (q == q_t) { S[k] <- 0; next }
    a <- min(q, q_t); b <- max(q, q_t)
    qq <- seq(a, b, length.out=n_int)
    integrand <- sqrt(2 * abs(V_fn(qq) - E))
    dq_int <- diff(qq)[1]
    S[k] <- sum(integrand) * dq_int
  }
  S
}

# ------------------------------------------------------------------------------
# LANGER DENSITY ON ONE SIDE OF A TURNING POINT
#
# rho(q) = |p(q)|^{-1} * z(q)^{1/2} * Ai^2(arg)
# arg = -z in allowed (V<E), giving oscillating Ai
# arg = +z in forbidden (V>E), giving decaying Ai
#
# Cells with abs_p == 0 (literally on a turning point) are interpolated
# from neighbours; the analytic limit is finite.
# ------------------------------------------------------------------------------

.langer_density_one_side <- function(q_grid, q_t, V_fn, E, hbar=1.0) {
  V_vals <- V_fn(q_grid)
  S      <- .action_from_turning_point(q_grid, q_t, V_fn, E)
  z      <- ((3 * S) / (2 * hbar))^(2/3)

  airy_arg <- ifelse(V_vals < E, -z, z)
  Ai_val   <- gsl::airy_Ai(airy_arg)

  abs_p <- sqrt(2 * abs(V_vals - E))
  prefactor <- ifelse(abs_p > 0, sqrt(z) / abs_p, NA_real_)

  rho <- prefactor * Ai_val^2

  if (any(is.na(rho))) {
    good_idx <- which(!is.na(rho))
    if (length(good_idx) >= 2) {
      rho <- approx(x=q_grid[good_idx], y=rho[good_idx],
                    xout=q_grid, rule=2)$y
    } else {
      rho[is.na(rho)] <- 0
    }
  }
  rho
}

# ------------------------------------------------------------------------------
# DOUBLE-TURNING-POINT DENSITY (one orbit)
# Standard "double-Airy patched at the orbit midpoint" construction.
# ------------------------------------------------------------------------------

.langer_density_one_orbit <- function(q_grid, q_minus, q_plus, V_fn, E,
                                      hbar=1.0) {
  q_mid <- (q_minus + q_plus) / 2
  left_mask  <- q_grid <= q_mid
  right_mask <- q_grid >  q_mid
  rho <- numeric(length(q_grid))

  if (any(left_mask)) {
    rho[left_mask]  <- .langer_density_one_side(
      q_grid[left_mask], q_minus, V_fn, E, hbar=hbar)
  }
  if (any(right_mask)) {
    rho[right_mask] <- .langer_density_one_side(
      q_grid[right_mask], q_plus, V_fn, E, hbar=hbar)
  }
  rho
}

# ------------------------------------------------------------------------------
# PUBLIC API
# ------------------------------------------------------------------------------

#' Berry-Mount / Langer uniform Airy density on a q grid.
#'
#' Approximates |psi_n(q)|^2 across the whole orbit using the Langer
#' uniform formula. Oscillating with n+1 lobes; finite at turning points
#' (Airy lobe replaces the WKB divergence); exponentially decaying in
#' the forbidden region.
#'
#' Apples-to-apples with the symplectic figure's right column: both
#' figures' right columns approximate the oscillating |psi_n|^2 with
#' the WKB caustic divergence resolved.
#'
#' @param q_grid          Output q grid (uniform spacing).
#' @param E               Energy.
#' @param V_fn            Function(q) returning potential.
#' @param turning_points  list(q_minus, q_plus) OR numeric vector of
#'                        turning points in ascending order (length 2,4,...).
#' @param hbar            Planck constant in chosen units.
#' @return Numeric vector of |psi_Langer|^2 on q_grid, normalized so
#'         integral = 1 over the grid.
airy_uniform_density <- function(q_grid, E, V_fn, turning_points, hbar=1.0) {
  if (is.list(turning_points)) {
    tp <- c(turning_points$q_minus, turning_points$q_plus)
  } else {
    tp <- sort(as.numeric(turning_points))
  }
  if (length(tp) %% 2 != 0)
    stop("airy_uniform_density: need an even number of turning points")

  rho <- numeric(length(q_grid))
  for (k in seq(1, length(tp), by=2)) {
    rho <- rho + .langer_density_one_orbit(q_grid, tp[k], tp[k+1],
                                           V_fn, E, hbar=hbar)
  }

  dq <- diff(q_grid)[1]
  norm <- sum(rho) * dq
  if (norm > 0) rho <- rho / norm

  rho
}
