# ==============================================================================
# classical_action_tools.R
# Classical phase-space geometry: orbits, Bohr-Sommerfeld action, and
# orbit-derived covariance (the semiclassical kinematic envelope).
#
# This file is purely classical: no wavefunctions, no Hamiltonian operators,
# no Schroedinger machinery. It computes orbit-level quantities at energy E
# given a potential V(q) and turning points.
#
# The orbit covariance below is the semiclassical analogue of the Robertson-
# Schroedinger covariance. For a quantum eigenstate the RS covariance is
# computed from psi (see quantum_tools.R / wigner-side machinery). For the
# semiclassical microcanonical state on the energy shell, the covariance
# is given by time-averages over the classical orbit.
#
# These two routes agree exactly for the harmonic oscillator and approach
# each other as O(1/n) for anharmonic potentials.
#
# Reference: Goldstein, Poole & Safko Classical Mechanics 3rd ed. Ch.10
#            Landau & Lifshitz Mechanics Ch.7, Sec. 49
#            de Gosson, The Principles of Newtonian and Quantum Mechanics
#            (Imperial College Press 2001) Ch.1
# Author: Brian S. Mulloy
# ==============================================================================

# ------------------------------------------------------------------------------
# TURNING POINT NORMALIZATION
# Internal helper: accept either list(q_minus, q_plus) or numeric vector.
# Returns sorted numeric vector.
# ------------------------------------------------------------------------------

.normalize_turning_points <- function(turning_points) {
  if (is.list(turning_points) &&
      !is.null(turning_points$q_minus) && !is.null(turning_points$q_plus)) {
    tp <- c(turning_points$q_minus, turning_points$q_plus)
  } else if (is.numeric(turning_points)) {
    tp <- sort(turning_points)
  } else {
    stop("turning_points must be list(q_minus, q_plus) or numeric vector")
  }
  if (length(tp) %% 2 != 0)
    stop(sprintf("Need even number of turning points, got %d", length(tp)))
  tp
}

# ------------------------------------------------------------------------------
# CLASSICAL TRAJECTORY
# ------------------------------------------------------------------------------

#' Classical phase-space trajectory at energy E.
#'
#' Accepts either list(q_minus, q_plus) or numeric vector of turning points.
#' For 2 turning points: returns one closed orbit.
#' For 4 turning points: returns two closed orbits separated by an NA row
#'   so geom_path breaks the line between them.
classical_trajectory <- function(V_fn, E, turning_points, n_pts=500) {
  tp <- .normalize_turning_points(turning_points)

  build_orbit <- function(q_lo, q_hi, group_id) {
    q_seq <- seq(q_lo, q_hi, length.out=n_pts)
    p_sq  <- pmax(2*(E - V_fn(q_seq)), 0)
    p_pos <- sqrt(p_sq)
    data.frame(
      q     = c(q_seq, rev(q_seq), q_seq[1]),
      p     = c(p_pos, -rev(p_pos), p_pos[1]),
      group = group_id
    )
  }

  na_separator <- data.frame(q=NA_real_, p=NA_real_, group=NA_integer_)

  orbits <- list()
  for (k in seq(1, length(tp), by=2)) {
    if (length(orbits) > 0) orbits[[length(orbits)+1]] <- na_separator
    orbits[[length(orbits)+1]] <- build_orbit(tp[k], tp[k+1], (k+1)/2)
  }
  do.call(rbind, orbits)
}

# ------------------------------------------------------------------------------
# BOHR-SOMMERFELD ACTION
# ------------------------------------------------------------------------------

#' Classical Bohr-Sommerfeld action A_BS/A_0 in units of A_0 = h/2.
classical_action <- function(V_fn, E, turning_points, hbar=1.0, n_pts=1001) {
  tp <- .normalize_turning_points(turning_points)

  total_area <- 0
  for (k in seq(1, length(tp), by=2)) {
    q_seq <- seq(tp[k], tp[k+1], length.out=n_pts)
    dq    <- diff(q_seq)[1]
    p_sq  <- pmax(2*(E - V_fn(q_seq)), 0)
    p_pos <- sqrt(p_sq)
    total_area <- total_area + 2 * sum(p_pos) * dq
  }
  total_area / (pi * hbar)
}

# ------------------------------------------------------------------------------
# ORBIT COVARIANCE (semiclassical kinematic envelope)
#
# Time-averaged second moments of the microcanonical phase-space distribution
# rho_cl(q,p) ~ delta(H(q,p) - E) on a closed orbit. The natural time measure
# is dt = dq/|p(q)|, so for any orbital observable f(q,p):
#
#   <f>_orbit = oint f(q,p) dt / oint dt = oint f(q,p) dq/|p| / oint dq/|p|
#
# For a single orbit between turning points q_- and q_+:
#   T(E)        = 2 * integral_{q_-}^{q_+} dq / |p(q)|       (period)
#   <q>_orbit   = (2/T) * integral_{q_-}^{q_+} q dq / |p(q)|
#   <q^2>_orbit = (2/T) * integral_{q_-}^{q_+} q^2 dq / |p(q)|
#   <p^2>_orbit = (2/T) * integral_{q_-}^{q_+} |p(q)| dq    (since |p| dt = dq)
#   <p>_orbit   = 0  (time-symmetry of bound orbit)
#
# Variances:
#   sigma_qq = <q^2>_orbit - <q>_orbit^2
#   sigma_pp = <p^2>_orbit
#
# Widths:
#   Delta_q = sqrt(2 * sigma_qq), Delta_p = sqrt(2 * sigma_pp)
#
# Numerical note: |p(q)| -> 0 at turning points produces an integrable 1/sqrt
# singularity in 1/|p|. Trapezoidal integration on a fine grid converges.
# The integrand for <p^2> contains |p| (no singularity) and converges trivially.
# ------------------------------------------------------------------------------

#' Semiclassical orbit covariance at energy E.
#'
#' Computes (Delta_q, Delta_p) and the Fermi-blob area A from time-averaged
#' second moments on the classical orbit. Returns the same structure as
#' robertson_schroedinger() in the quantum pipeline so call sites can be
#' interchanged.
#'
#' For multi-orbit potentials (e.g., sub-barrier double-well with 4 turning
#' points), the function combines moments across all orbital segments
#' weighted by their respective periods, treating the union as a single
#' microcanonical distribution at energy E.
#'
#' @param V_fn Function(q) returning potential values.
#' @param E Energy.
#' @param turning_points list(q_minus, q_plus) or sorted numeric vector.
#' @param hbar Planck constant in chosen units.
#' @param n_pts Number of grid points per orbit segment for quadrature.
#' @return Named list with the same fields as robertson_schroedinger():
#'         Delta_q, Delta_p, delta_q, delta_p, A_over_A0, plus an extra
#'         flag heisenberg_satisfied for the diagnostic A >= h/2.
#'         A warning is issued if the semiclassical envelope falls below
#'         the Heisenberg bound; this can happen at low quantum numbers
#'         in anharmonic potentials and indicates the regime where the
#'         semiclassical approximation breaks down. The returned widths
#'         remain operationally usable.
orbit_covariance <- function(V_fn, E, turning_points, hbar=1.0, n_pts=10001) {
  tp <- .normalize_turning_points(turning_points)

  # Accumulate weighted moments across all orbital segments.
  # T_total: total orbital period (sum across segments)
  # M_q, M_q2, M_p2: numerators of <q>, <q^2>, <p^2> (sum across segments)
  T_total <- 0
  M_q     <- 0
  M_q2    <- 0
  M_p2    <- 0

  for (k in seq(1, length(tp), by=2)) {
    q_seq   <- seq(tp[k], tp[k+1], length.out=n_pts)
    dq      <- diff(q_seq)[1]
    p_sq    <- pmax(2*(E - V_fn(q_seq)), 0)
    # Floor avoids 1/0 at turning points; the 1/sqrt singularity is
    # integrable so trapezoidal converges as floor -> 0.
    p_pos   <- sqrt(pmax(p_sq, 1e-15))
    inv_p   <- 1 / p_pos

    # Period of this segment: T_seg = 2 * int dq/|p|
    T_seg   <- 2 * sum(inv_p) * dq
    # Moments on this segment, weighted by 2*dq/|p| (forward + backward leg)
    M_q_seg  <- 2 * sum(q_seq      * inv_p) * dq
    M_q2_seg <- 2 * sum(q_seq^2    * inv_p) * dq
    # <p^2>: integrand p^2 / |p| = |p|, weighted by 2*dq for the two legs
    M_p2_seg <- 2 * sum(p_pos)              * dq

    T_total <- T_total + T_seg
    M_q     <- M_q     + M_q_seg
    M_q2    <- M_q2    + M_q2_seg
    M_p2    <- M_p2    + M_p2_seg
  }

  # Time-averages
  q_mean   <- M_q  / T_total
  q2_mean  <- M_q2 / T_total
  p2_mean  <- M_p2 / T_total

  sigma_qq <- q2_mean - q_mean^2
  sigma_pp <- p2_mean
  # By time-symmetry of bound orbits: <p> = 0 and sigma_qp = 0

  Delta_q  <- sqrt(2 * sigma_qq)
  Delta_p  <- sqrt(2 * sigma_pp)
  delta_q  <- hbar / Delta_p
  delta_p  <- hbar / Delta_q

  # A_0 = pi * hbar (= h/2). Fermi-blob area A = pi * Delta_q * Delta_p.
  # Diagnostic: at low quantum numbers in anharmonic potentials, A may
  # come out fractionally below h/2 (the orbit is tighter than the
  # quantum minimum-uncertainty state). Issue a warning but do not cap.
  A_over_A0 <- (Delta_q * Delta_p) / hbar
  heisenberg_satisfied <- A_over_A0 >= 1
  if (!heisenberg_satisfied) warning(sprintf(
    paste("orbit_covariance: A/A0 = %.4f < 1 at E = %.4f.",
          "Semiclassical regime breaks down at low quantum numbers in",
          "anharmonic potentials; widths remain usable but the symplectic",
          "interpretation as conjugate quantum blobs is operative only",
          "when A >= h/2."), A_over_A0, E))

  list(
    Delta_q              = Delta_q,
    Delta_p              = Delta_p,
    delta_q              = delta_q,
    delta_p              = delta_p,
    A_over_A0            = A_over_A0,
    q_mean               = q_mean,
    heisenberg_satisfied = heisenberg_satisfied
  )
}
