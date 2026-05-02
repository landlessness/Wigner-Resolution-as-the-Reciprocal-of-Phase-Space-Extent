# ==============================================================================
# cat_system.R
# n-component cat states: coherent superposition of n coherent states placed
# in phase space. Cat states are non-eigenstate quantum objects with no
# classical orbit and no Bohr-Sommerfeld condition; they exist purely in
# the quantum universe of this paper.
#
# Three configurations supported, all anchored so adjacent cats sit at
# distance 2*p_max from each other:
#   n_cats=2:           lobes at (0, +5), (0, -5)
#                       -- N/S
#   n_cats=3:           lobes at (0, +5), (+10/sqrt(3), -5),
#                                         (-10/sqrt(3), -5)
#                       -- triangle pointing up; equilateral with side
#                          length 10
#   n_cats=4 (diag):    lobes at (+5, +5), (-5, +5), (-5, -5), (+5, -5)
#                       -- compass with lobes on the diagonals; side
#                          length 10
#   n_cats=4 (axis):    lobes at (+5sqrt(2), 0), (0, +5sqrt(2)),
#                                (-5sqrt(2), 0), (0, -5sqrt(2))
#                       -- the diag compass rotated 45 deg, with the
#                          same adjacent-lobe distance of 10. Lobes
#                          reach further out in q and p (extreme
#                          coordinate p_max*sqrt(2) ~ 7.07 instead of
#                          p_max=5), but the underlying physics is
#                          identical to diag up to phase-space
#                          rotation: same A_RS/A_0, same kernel widths,
#                          same QoA size.
#
# The p=0 horizontal cross-section cuts through:
#   n_cats=2:        central interference fringe between the two lobes.
#   n_cats=3:        two pair-fringe envelopes between top lobe and
#                    bottom pair.
#   n_cats=4 (diag): central chessboard plus two pair-fringe envelopes
#                    at (+/-5, 0) -- the NE+SE and NW+SW pair midpoints.
#   n_cats=4 (axis): the E and W lobes themselves (which sit on p=0),
#                    plus three pair-fringe envelopes.
#
# The cat wavefunction is a coherent superposition of Gaussian wavepackets:
#   psi(q) = sum_k exp(-(q - q_k)^2 / (2 xi^2)) * exp(i p_k q / hbar)
# normalized so int |psi|^2 dq = 1.
#
# Reference: Schleich, "Quantum Optics in Phase Space" (Wiley, 2001) Ch. 7;
#            Zurek 2001 Nature 412, 712 for the original compass state.
# Author: Brian S. Mulloy
# ==============================================================================

# ------------------------------------------------------------------------------
# DEFAULT PARAMETERS
# ------------------------------------------------------------------------------

CAT_P_MAX  <- 5.0   # extreme p coordinate of any lobe (always +/- P_MAX)
CAT_XI     <- 1.0   # coherent-state position width
CAT_HBAR   <- 1.0   # natural units

# Display window for cat figures. Wide enough to enclose all lobes plus
# their interference envelopes, with room to show cross-section structure.
CAT_Q_MIN  <- -13.0
CAT_Q_MAX  <-  13.0
CAT_DQ     <-   0.02
CAT_Q_DISPLAY <- 10.0   # display half-window for figures

# ------------------------------------------------------------------------------
# LOBE GEOMETRY
# ------------------------------------------------------------------------------

#' Compute lobe positions (q_k, p_k) for an n-cat configuration.
#'
#' Lobes are anchored so that adjacent cats sit at distance 2*p_max
#' from each other. For n_cats=2 and n_cats=3 the extreme p coordinate
#' is +/- p_max. For n_cats=4 (diag) the extreme coordinate is also
#' +/- p_max (lobes on the +/-p_max square corners). For n_cats=4
#' (axis) the extreme coordinate is +/- p_max*sqrt(2) (lobes on the
#' axes) so that the adjacent-cat distance still equals 2*p_max --
#' axis is the literal 45 deg rotation of diag.
#'
#' @param n_cats  Number of coherent states (2, 3, or 4).
#' @param variant For n_cats=4 only: "diag" (compass with lobes on
#'                the diagonals) or "axis" (compass rotated 45 deg
#'                with lobes on the q and p axes). Ignored otherwise.
#' @param p_max   Sets the adjacent-cat distance to 2*p_max.
#' @return list(q=numeric, p=numeric).
cat_lobe_positions <- function(n_cats, variant="diag", p_max=CAT_P_MAX) {
  if (n_cats == 2) {
    # N, S
    list(q = c(0,    0),
         p = c(p_max, -p_max))
  } else if (n_cats == 3) {
    # Triangle pointing up: apex at (0, +p_max), base lobes at (+/- q_b, -p_max)
    # with q_b = 2 p_max / sqrt(3) so all three sides of the triangle are equal.
    q_b <- 2 * p_max / sqrt(3)
    list(q = c(0,    q_b,   -q_b),
         p = c(p_max, -p_max, -p_max))
  } else if (n_cats == 4) {
    if (variant == "diag") {
      # Compass with lobes at the corners of a 2 p_max x 2 p_max square.
      list(q = c(p_max, -p_max, -p_max,  p_max),
           p = c(p_max,  p_max, -p_max, -p_max))
    } else if (variant == "axis") {
      # 45 deg rotation of diag. Adjacent-cat distance is still 2*p_max,
      # but achieved by placing lobes at +/- p_max*sqrt(2) on the axes.
      r <- p_max * sqrt(2)
      list(q = c(0,  r,  0, -r),
           p = c(r,  0, -r,  0))
    } else {
      stop(sprintf("cat_lobe_positions: variant must be 'axis' or 'diag', got '%s'",
                   variant))
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
#' normalized so sum |psi|^2 dq = 1.
#'
#' @param q       Position grid (numeric vector).
#' @param n_cats  Number of coherent states.
#' @param p_max   Extreme p coordinate of lobes.
#' @param xi      Coherent-state position width.
#' @param hbar    Natural units.
#' @return Complex vector psi(q), L2-normalized on q.
cat_psi <- function(q, n_cats, variant="diag", p_max=CAT_P_MAX, xi=CAT_XI, hbar=CAT_HBAR) {
  lobes <- cat_lobe_positions(n_cats, variant=variant, p_max=p_max)
  psi   <- complex(real=rep(0, length(q)), imaginary=rep(0, length(q)))
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
