# ==============================================================================
# symplectic_tools.R
# Shared symplectic geometry tools for resolving Wigner negativity and
# semiclassical caustics via finite classical action.
#
# Units: all positions in units of q_0, momenta in units of p_0.
# This gives hbar = q_0 * p_0 = 1 and A_0 = pi * q_0 * p_0 = h/2.
# Action ratios A/A_0 = 2n+1 for the QHO, derived from Robertson-Schrödinger.
#
# Reference: de Gosson (2009), Zurek (2001), Robertson (1929), Schrödinger (1930)
# Author: Brian S. Mulloy
# ==============================================================================

library(gsl)
library(data.table)

# ------------------------------------------------------------------------------
# UNIT SYSTEM
# ------------------------------------------------------------------------------
# Position : q measured in units of q_0  (ground state position uncertainty)
# Momentum : p measured in units of p_0  (ground state momentum uncertainty)
# Action   : A measured in units of A_0 = pi * q_0 * p_0 = h/2
# In these units: hbar = q_0 * p_0 = 1, h = 2*pi, A_0 = pi
# All functions below assume this unit system unless hbar is passed explicitly.
# ------------------------------------------------------------------------------

#' Robertson-Schrödinger symplectic geometry from a covariance matrix.
#'
#' Computes classical kinematic boundaries Delta_q, Delta_p and Zurek
#' reciprocal scales delta_q, delta_p from the quantum covariance matrix.
#' Verifies the Robertson-Schrödinger inequality and the kernel symplectic
#' positivity condition delta_q * Delta_p = hbar.
#'
#' The Robertson-Schrödinger uncertainty relation is:
#'   sigma_qq * sigma_pp - sigma_qp^2 >= (hbar/2)^2
#'
#' For the QHO eigenstate n: sigma_qq = sigma_pp = (2n+1)/2, sigma_qp = 0.
#' For Morse and double-well: sigma_qp != 0 and Delta_q != Delta_p.
#'
#' @param sigma_qq Position variance <q^2> - <q>^2
#' @param sigma_pp Momentum variance <p^2> - <p>^2
#' @param sigma_qp Position-momentum covariance (<qp+pq>/2 - <q><p>)
#' @param hbar Reduced Planck constant (default 1.0 in q_0*p_0 units)
#' @return Named list: Delta_q, Delta_p, delta_q, delta_p, A_over_A0,
#'         rs_satisfied, sp_satisfied (symplectic positivity check)
robertson_schrodinger <- function(sigma_qq, sigma_pp, sigma_qp = 0, hbar = 1.0) {

  # Robertson-Schrödinger inequality check
  rs_lhs       <- sigma_qq * sigma_pp - sigma_qp^2
  rs_bound     <- (hbar / 2)^2
  rs_satisfied <- rs_lhs >= rs_bound - .Machine$double.eps * abs(rs_bound)

  if (!rs_satisfied) {
    warning(sprintf(
      "RS inequality violated: sigma_qq*sigma_pp - sigma_qp^2 = %.6e < (hbar/2)^2 = %.6e",
      rs_lhs, rs_bound
    ))
  }

  # Classical semi-axes from 1-sigma Fermi blob (de Gosson convention)
  Delta_q <- sqrt(2 * sigma_qq)
  Delta_p <- sqrt(2 * sigma_pp)

  # Zurek reciprocal scales
  delta_q <- hbar / Delta_p
  delta_p <- hbar / Delta_q

  # Symplectic positivity check: kernel saturates Hudson bound iff delta_q * Delta_p = hbar
  sp_product   <- delta_q * Delta_p
  sp_satisfied <- abs(sp_product - hbar) < .Machine$double.eps^0.5 * hbar

  if (!sp_satisfied) {
    warning(sprintf(
      "Kernel symplectic positivity not saturated: delta_q * Delta_p = %.6e (expected hbar = %.6e)",
      sp_product, hbar
    ))
  }

  # Classical action ratio A/A_0
  # A = pi * Delta_q * Delta_p, A_0 = pi * hbar, so A/A_0 = Delta_q * Delta_p / hbar
  A_over_A0 <- (Delta_q * Delta_p) / hbar

  list(
    Delta_q      = Delta_q,
    Delta_p      = Delta_p,
    delta_q      = delta_q,
    delta_p      = delta_p,
    A_over_A0    = A_over_A0,
    rs_satisfied = rs_satisfied,
    sp_satisfied = sp_satisfied
  )
}

#' QHO Robertson-Schrödinger geometry for eigenstate n.
#'
#' Exact covariance for harmonic oscillator |n>:
#'   sigma_qq = sigma_pp = (2n+1) * hbar / 2,  sigma_qp = 0
#' giving A/A_0 = 2n+1.
#'
#' @param n Non-negative integer quantum number
#' @param hbar Reduced Planck constant (default 1.0)
#' @return Output of robertson_schrodinger()
qho_covariance <- function(n, hbar = 1.0) {
  alpha    <- 2*n + 1
  sigma_qq <- alpha * hbar / 2
  sigma_pp <- alpha * hbar / 2
  robertson_schrodinger(sigma_qq, sigma_pp, sigma_qp = 0, hbar = hbar)
}

#' True QHO Wigner function W_n(q, p).
#'
#' Exact formula in q_0, p_0 units (hbar = 1):
#'   W_n(q,p) = (-1)^n / pi * exp(-(q^2 + p^2)) * L_n(2*(q^2 + p^2))
#'
#' Uses GSL laguerre_n() for numerical stability at large n.
#' Spot-check values: W_n(0,0) = (-1)^n / pi.
#'
#' @param n Quantum number
#' @param q Position in units of q_0 (scalar or vector)
#' @param p Momentum in units of p_0 (scalar or vector)
#' @return W_n(q, p)
qho_wigner <- function(n, q, p) {
  rho2 <- q^2 + p^2
  (-1)^n / pi * exp(-rho2) * laguerre_n(n, 0, 2 * rho2)
}

#' Squeezed kernel G_delta_q(q, p) in q_0, p_0 units.
#'
#' Squeezed in q (width delta_q), extended in p (width Delta_p).
#' Symplectic area = h/2, saturating Hudson's bound.
#' In q_0, p_0 units (hbar = 1): prefactor = 1/pi.
#'
#' @param q Position values
#' @param p Momentum values
#' @param rs Output of robertson_schrodinger() or qho_covariance()
#' @return G_delta_q(q, p)
squeezed_kernel_q <- function(q, p, rs) {
  (1/pi) * exp(-q^2 / rs$delta_q^2 - p^2 / rs$Delta_p^2)
}

#' Conjugate squeezed kernel G_delta_p(q, p) in q_0, p_0 units.
#'
#' Squeezed in p (width delta_p), extended in q (width Delta_q).
#' Symplectic area = h/2, saturating Hudson's bound.
#'
#' @param q Position values
#' @param p Momentum values
#' @param rs Output of robertson_schrodinger() or qho_covariance()
#' @return G_delta_p(q, p)
squeezed_kernel_p <- function(q, p, rs) {
  (1/pi) * exp(-q^2 / rs$Delta_q^2 - p^2 / rs$delta_p^2)
}

#' 2D FFT convolution P = W * K on a uniform grid.
#'
#' Uses ifftshift on the kernel so its center maps to the FFT origin [1,1].
#' Hudson's theorem guarantees P >= 0 when K has symplectic area >= h/2.
#' Negative values below the grid-resolution tolerance are floating point
#' artifacts and are zeroed. Values exceeding tolerance indicate a grid
#' resolution problem and trigger a warning.
#'
#' @param W_mat Wigner function matrix (nq x np)
#' @param K_mat Kernel matrix (nq x np, unnormalized)
#' @param dq Grid spacing in q
#' @param dp Grid spacing in p
#' @return Named list: P_mat, max_negative, tolerance
fft_convolve_2d <- function(W_mat, K_mat, dq, dp) {
  nq <- nrow(W_mat)
  np <- ncol(W_mat)

  # Normalize kernel to unit integral
  K_norm <- K_mat / (sum(K_mat) * dq * dp)

  # ifftshift: move kernel peak to FFT origin [1,1]
  ifftshift2d <- function(m) {
    nr <- nrow(m); nc <- ncol(m)
    sr <- floor(nr/2); sc <- floor(nc/2)
    rbind(
      cbind(m[(sr+1):nr, (sc+1):nc], m[(sr+1):nr, 1:sc]),
      cbind(m[1:sr,      (sc+1):nc], m[1:sr,      1:sc])
    )
  }
  K_shift <- ifftshift2d(K_norm)

  # FFT convolution theorem
  P_mat <- Re(fft(fft(W_mat) * fft(K_shift), inverse = TRUE)) / (nq * np)
  P_mat <- P_mat * dq * dp

  # Tolerance scaled to grid size — FFT floating point error is O(N * eps)
  peak_val     <- max(abs(P_mat))
  tol          <- peak_val * sqrt(.Machine$double.eps) * sqrt(nq * np)
  max_negative <- min(P_mat)

  if (max_negative < -tol) {
    warning(sprintf(
      "fft_convolve_2d: min(P_mat) = %.2e exceeds tolerance %.2e. Increase integ_res.",
      max_negative, -tol
    ))
  }

  # Zero sub-tolerance negatives — floating point artifacts, not physics
  P_mat[P_mat < 0 & P_mat >= -tol] <- 0

  list(P_mat = P_mat, max_negative = max_negative, tolerance = tol)
}

#' Compute P_delta_q(q): 2D Wigner convolution projected to position space.
#'
#' Pipeline:
#'   1. Build W_n on a wide integration grid to capture Gaussian tails
#'   2. Verify Wigner normalization: integral W dq dp = 1
#'   3. Build kernel G from Robertson-Schrödinger widths
#'   4. Verify kernel symplectic positivity: delta_q * Delta_p = hbar
#'   5. Convolve W * G via FFT (Hudson guarantees result >= 0)
#'   6. Project to 1D: P_q(q) = integral P(q,p) dp
#'   7. Normalize exactly once on the integration grid
#'   8. Interpolate onto display grid (view only, no renormalization)
#'   9. Log W_n(0,0) spot-check against analytic value (-1)^n / pi
#'
#' @param n Quantum number
#' @param wigner_fn Function(n, q, p) -> Wigner values
#' @param kernel_fn Function(q, p, rs) -> kernel values
#' @param rs Output of robertson_schrodinger() or qho_covariance()
#' @param q_display Display grid for interpolation
#' @param hbar Reduced Planck constant (default 1.0)
#' @return Named list: P_q, w_norm, w_spot_check, max_negative, tolerance
compute_symplectic_density <- function(n, wigner_fn, kernel_fn, rs,
                                       q_display, hbar = 1.0) {
  Delta_q <- rs$Delta_q

  # Integration grid: wide enough to capture Gaussian tails
  integ_lim <- max(Delta_q * 2.0, 4.0 * sqrt(hbar))
  integ_res <- max(601, 20 * ceiling(rs$A_over_A0))
  if (integ_res %% 2 == 0) integ_res <- integ_res + 1

  q_int <- seq(-integ_lim, integ_lim, length.out = integ_res)
  p_int <- seq(-integ_lim, integ_lim, length.out = integ_res)
  dq    <- diff(q_int)[1]
  dp    <- diff(p_int)[1]

  # Build Wigner matrix
  W_mat <- outer(q_int, p_int, FUN = function(q, p) wigner_fn(n, q, p))

  # Normalization check
  w_norm <- sum(W_mat) * dq * dp
  if (abs(w_norm - 1) > 1e-3) {
    warning(sprintf(
      "n=%d: Wigner norm = %.6f. Increase integ_lim or integ_res.",
      n, w_norm
    ))
  }

  # Spot-check W_n(0,0): analytic value is (-1)^n / pi
  w_at_origin  <- wigner_fn(n, 0, 0)
  w_expected   <- (-1)^n / pi
  w_spot_check <- abs(w_at_origin - w_expected)
  if (w_spot_check > 1e-6) {
    warning(sprintf(
      "n=%d: W_n(0,0) = %.6f, expected (-1)^n/pi = %.6f (error %.2e)",
      n, w_at_origin, w_expected, w_spot_check
    ))
  }

  # Build kernel and convolve
  K_mat <- outer(q_int, p_int, FUN = function(q, p) kernel_fn(q, p, rs))
  conv  <- fft_convolve_2d(W_mat, K_mat, dq, dp)
  P_mat <- conv$P_mat

  # Project to 1D
  P_q_int  <- rowSums(P_mat) * dp
  norm_int <- sum(P_q_int) * dq

  if (abs(norm_int) < .Machine$double.eps^0.5) {
    stop(sprintf("n=%d: P_delta_q integrates to zero. Check wigner_fn and kernel_fn.", n))
  }

  # Single normalization on integration grid
  P_q_int <- P_q_int / norm_int

  # Interpolate onto display grid — display is a view, not a physics operation
  P_q <- approx(q_int, P_q_int, xout = q_display, rule = 1)$y
  P_q[is.na(P_q)] <- 0

  list(
    P_q          = P_q,
    w_norm       = w_norm,
    w_spot_check = w_spot_check,
    max_negative = conv$max_negative,
    tolerance    = conv$tolerance
  )
}
