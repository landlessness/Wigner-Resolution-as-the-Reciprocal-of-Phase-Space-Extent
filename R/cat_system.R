# ==============================================================================
# cat_system.R
# n-component cat states: coherent superposition of n coherent states placed
# in phase space such that the cross-section at p=0 cuts through canonical
# interference structure with visible Wigner negativity.
#
# Four configurations supported, all anchored so that the extreme p
# coordinates of the lobes are at +/- P_MAX. The q coordinates adjust to
# maintain equally-spaced rotational symmetry around the origin:
#   n_cats=2:                       lobes at (0, +5), (0, -5).
#   n_cats=3:                       lobes at (0, +5), (+10/sqrt(3), -5),
#                                          (-10/sqrt(3), -5).
#                                   Triangle pointing up with apex at
#                                   (0, +5) and base at p = -5.
#   n_cats=4 orientation="cardinal" (default):
#                                   lobes at (0, +5), (+5, 0), (0, -5), (-5, 0).
#                                   Zurek cardinal compass.
#   n_cats=4 orientation="diagonal":
#                                   lobes at (+5, +5), (-5, +5), (-5, -5), (+5, -5).
#                                   Zurek compass rotated 45 deg
#                                   (lobes on diagonals, extreme p at +/-5).
#
# The p=0 horizontal cross-section cuts through:
#   n_cats=2: central interference fringe between the two lobes.
#   n_cats=3: two pair-fringe envelopes between top lobe and bottom pair.
#   n_cats=4 cardinal: W and E lobes plus the central chessboard.
#   n_cats=4 diagonal: central chessboard plus two pair-fringe envelopes
#                      at (+/- 5, 0) (NE+SE and NW+SW pair midpoints).
#
# Each lobe is a coherent state at phase-space location (q_k, p_k) with
# Gaussian width xi in q. The cat wavefunction in position space is
#   psi(q) = sum_k exp(-(q - q_k)^2 / (2 xi^2)) * exp(i p_k q / hbar)
# normalized so integral |psi|^2 dq = 1.
#
# In general psi(q) is complex-valued (the wigner_fft pipeline handles this).
# It is real-valued only when lobes are mirror-symmetric across the q-axis
# (n_cats=2, n_cats=4 cardinal, n_cats=4 diagonal). For the n_cats=3
# triangle-pointing-up configuration, psi has a non-trivial imaginary part.
#
# Reference: Schleich, "Quantum Optics in Phase Space" (Wiley, 2001) Ch. 7
#            for n-cat states; Wikipedia "Cat state" for the canonical
#            three-row visualization (n=2,3,4); Zurek 2001 Nature 412, 712
#            for the original compass state (cardinal orientation).
# Author: Brian S. Mulloy
# ==============================================================================

# ------------------------------------------------------------------------------
# DEFAULT PARAMETERS
# Lobes are anchored such that the extreme p coordinates are at +/- CAT_P_MAX.
# This guarantees the p=0 horizontal cross-section sits halfway between the
# extreme-p lobes in every row.
# ------------------------------------------------------------------------------

CAT_P_MAX <- 5.0   # extreme p coordinate of lobes (always +/- P_MAX)
CAT_XI    <- 1.0   # coherent-state width

# ------------------------------------------------------------------------------
# LOBE GEOMETRY
# ------------------------------------------------------------------------------

#' Compute lobe positions (q_k, p_k) for an n-cat configuration.
#'
#' Returns a list with components q (vector of length n) and p (vector of
#' length n). Lobes are equally spaced rotationally around the origin and
#' anchored so that the extreme p coordinates are at +/- p_max.
#'
#' @param n_cats     Number of coherent states (2, 3, or 4).
#' @param p_max      Extreme p coordinate of lobes (default CAT_P_MAX).
#' @param orientation For n_cats=4: "cardinal" (default, lobes on cardinal
#'                   axes) or "diagonal" (lobes on diagonals, extreme p at
#'                   +/- p_max). Ignored for n_cats=2 and n_cats=3.
#' @return list(q=numeric, p=numeric).
cat_lobe_positions <- function(n_cats, p_max=CAT_P_MAX, orientation="cardinal") {
  if (n_cats == 2) {
    # N, S: (0, +p_max), (0, -p_max).
    list(q=c(0, 0), p=c(p_max, -p_max))
  } else if (n_cats == 3) {
    # Triangle pointing up: apex at (0, +p_max), base lobes at (+/- q_b, -p_max).
    # Equal spacing 120 deg around the origin requires base q such that
    # adjacent-lobe distance is the same on all three sides:
    #   apex-to-base  = sqrt(q_b^2 + (2 p_max)^2)
    #   base-to-base  = 2 q_b
    # Setting equal -> q_b = 2 p_max / sqrt(3).
    q_b <- 2 * p_max / sqrt(3)
    list(q=c(0, q_b, -q_b), p=c(p_max, -p_max, -p_max))
  } else if (n_cats == 4) {
    if (orientation == "cardinal") {
      # Zurek cardinal compass: lobes on cardinal axes at radius p_max.
      list(q=c(0, p_max, 0, -p_max), p=c(p_max, 0, -p_max, 0))
    } else if (orientation == "diagonal") {
      # Diagonal compass: lobes at (+/- p_max, +/- p_max). Extreme p still
      # at +/- p_max; extreme q also at +/- p_max.
      list(q=c(p_max, -p_max, -p_max,  p_max),
           p=c(p_max,  p_max, -p_max, -p_max))
    } else {
      stop(sprintf(
        "cat_lobe_positions: orientation must be 'cardinal' or 'diagonal' (got '%s')",
        orientation))
    }
  } else {
    stop("cat_lobe_positions: n_cats must be 2, 3, or 4")
  }
}

# ------------------------------------------------------------------------------
# CAT WAVEFUNCTION
# ------------------------------------------------------------------------------

#' Complex-valued n-cat wavefunction psi(q) sampled on a grid.
#'
#' psi(q) = sum_k exp(-(q - q_k)^2 / (2 xi^2)) * exp(i p_k q / hbar)
#' with normalization sum |psi|^2 dq = 1.
#'
#' @param q           Position grid (numeric vector).
#' @param n_cats      Number of coherent states (2, 3, or 4).
#' @param p_max       Extreme p coordinate of lobes (default CAT_P_MAX).
#' @param xi          Coherent-state width (default CAT_XI).
#' @param hbar        Planck constant (default 1.0).
#' @param orientation For n_cats=4: "cardinal" or "diagonal".
#' @return Complex vector psi(q), L2-normalized on q.
cat_psi <- function(q, n_cats, p_max=CAT_P_MAX, xi=CAT_XI, hbar=1.0,
                    orientation="cardinal") {
  lobes <- cat_lobe_positions(n_cats, p_max, orientation)
  psi <- complex(real=rep(0, length(q)), imaginary=rep(0, length(q)))
  for (k in seq_along(lobes$q)) {
    q_k <- lobes$q[k]
    p_k <- lobes$p[k]
    psi <- psi + exp(-(q - q_k)^2 / (2*xi^2)) * exp(1i * p_k * q / hbar)
  }
  dq   <- diff(q)[1]
  norm <- sqrt(sum(abs(psi)^2) * dq)
  if (norm > 0) psi <- psi / norm
  psi
}

# ------------------------------------------------------------------------------
# DISPLAY GRID PARAMETERS
# Wide enough to enclose all lobes (extending to about q ~ 5.77 for n=3) plus
# their interference envelopes, with room to show the cross-section structure.
# ------------------------------------------------------------------------------

CAT_Q_MIN    <- -10.0
CAT_Q_MAX    <-  10.0
CAT_DQ       <-   0.02
