# ==============================================================================
# classical_action_tools.R
# Classical action computation and phase-space geometry.
#
# Reference: Goldstein, Poole & Safko Classical Mechanics 3rd ed. Ch.10
#            Landau & Lifshitz Mechanics Ch.7
#            de Gosson The Principles of Newtonian and Quantum Mechanics
#            (Imperial College Press 2001) Ch.1
# Author: Brian S. Mulloy
# ==============================================================================

#' Classical phase-space trajectory at energy E.
#' Accepts either a single (q_minus, q_plus) pair or a vector of turning points.
#' For 2 turning points: returns one closed orbit.
#' For 4 turning points: returns two closed orbits separated by an NA row
#'   so geom_path breaks the line between them.
classical_trajectory <- function(V_fn, E, turning_points, n_pts=500) {
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

#' Classical Bohr-Sommerfeld action A_BS/A_0 in units of A_0 = h/2.
classical_action <- function(V_fn, E, turning_points, hbar=1.0, n_pts=1001) {
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
# ROBERTSON-SCHRODINGER COVARIANCE
# RS check tolerance is sqrt(machine_eps) ~ 1e-8 to accommodate finite-grid
# discretization of the momentum derivative; analytic ground states saturate
# RS exactly, so the discrete computation may sit fractionally below the
# continuous bound by O(dq^2).
# ------------------------------------------------------------------------------

# Tolerance for RS-bound check: relative deviation from saturation that we
# accept as numerically equivalent.
RS_TOLERANCE <- sqrt(.Machine$double.eps)

#' RS geometry from covariance matrix elements.
robertson_schrodinger <- function(sigma_qq, sigma_pp, sigma_qp=0, hbar=1.0) {
  rs_lhs       <- sigma_qq*sigma_pp - sigma_qp^2
  rs_bound     <- (hbar/2)^2
  rs_satisfied <- rs_lhs >= rs_bound * (1 - RS_TOLERANCE)
  if (!rs_satisfied) warning(sprintf(
    "RS inequality violated: %.6e < %.6e", rs_lhs, rs_bound))
  Delta_q      <- sqrt(2*sigma_qq)
  Delta_p      <- sqrt(2*sigma_pp)
  delta_q      <- hbar/Delta_p
  delta_p      <- hbar/Delta_q
  sp_product   <- delta_q*Delta_p
  sp_satisfied <- abs(sp_product-hbar) < RS_TOLERANCE*hbar
  if (!sp_satisfied) warning(sprintf(
    "SP not saturated: %.6e != %.6e", sp_product, hbar))
  list(
    Delta_q      = Delta_q,
    Delta_p      = Delta_p,
    delta_q      = delta_q,
    delta_p      = delta_p,
    A_over_A0    = (Delta_q*Delta_p)/hbar,
    rs_satisfied = rs_satisfied,
    sp_satisfied = sp_satisfied
  )
}

#' RS covariance from numerically sampled wavefunction.
numerical_covariance <- function(psi_vec, q_grid, hbar=1.0) {
  dq       <- diff(q_grid)[1]
  norm     <- sqrt(sum(psi_vec^2)*dq)
  psi_n    <- psi_vec/norm
  prob     <- psi_n^2
  q_mean   <- sum(q_grid*prob)*dq
  sigma_qq <- sum((q_grid-q_mean)^2*prob)*dq
  dpsi     <- diff(psi_n)/dq
  sigma_pp <- hbar^2*sum(dpsi^2)*dq
  robertson_schrodinger(sigma_qq, sigma_pp, sigma_qp=0, hbar=hbar)
}
