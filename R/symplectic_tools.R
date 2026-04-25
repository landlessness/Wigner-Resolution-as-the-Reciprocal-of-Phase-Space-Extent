# ==============================================================================
# symplectic_tools.R
# Shared symplectic geometry tools for resolving Wigner negativity and
# semiclassical caustics via finite classical action.
#
# Reference: de Gosson (2009), Zurek (2001), Robertson (1929), Schrödinger (1930)
# Author: Brian S. Mulloy
# ==============================================================================

# Required packages
library(gsl)        # Laguerre polynomials via laguerre_n() — numerically stable
library(data.table)

# ------------------------------------------------------------------------------
# PHYSICAL CONSTANTS
# ------------------------------------------------------------------------------
# A_0 = pi * q_0 * p_0 = h/2 is the symplectic quantum of action (ground state)
# In natural units q_0 = p_0 = 1, so A_0 = pi/2... but we work in A/A_0 ratios
# so A_0 cancels throughout. alpha = A/A_0 = 2n+1 for the QHO.
# ------------------------------------------------------------------------------

#' Compute the Robertson-Schrödinger covariance matrix for a given state.
#'
#' The Robertson-Schrödinger uncertainty relation generalizes Heisenberg:
#'   sigma_qq * sigma_pp - sigma_qp^2 >= (hbar/2)^2
#'
#' For the QHO eigenstate n in natural units (hbar=1, q0=p0=1):
#'   sigma_qq = (2n+1)/2, sigma_pp = (2n+1)/2, sigma_qp = 0
#'
#' For non-circular systems (Morse, double-well), sigma_qp != 0 and
#' Delta_q != Delta_p. This function is the single source of truth for
#' all kernel width computations across QHO, Morse, and double-well.
#'
#' @param sigma_qq Variance in position (q^2 expectation minus mean^2)
#' @param sigma_pp Variance in momentum (p^2 expectation minus mean^2)
#' @param sigma_qp Covariance between q and p (0 for QHO eigenstates)
#' @param hbar Reduced Planck constant (default 1.0 in natural units)
#'
#' @return A named list with:
#'   Delta_q  : classical position semi-axis (sqrt(2) * sigma_q)
#'   Delta_p  : classical momentum semi-axis (sqrt(2) * sigma_p)
#'   delta_q  : squeezed position width = hbar / Delta_p
#'   delta_p  : squeezed momentum width = hbar / Delta_q
#'   A_over_A0: classical action in units of A_0 = h/2
#'   rs_bound : RS bound = hbar^2/4, must satisfy sigma_qq*sigma_pp - sigma_qp^2 >= rs_bound
#'   rs_satisfied: logical, TRUE if RS inequality is satisfied
robertson_schrodinger <- function(sigma_qq, sigma_pp, sigma_qp = 0, hbar = 1.0) {
  # Verify Robertson-Schrödinger inequality
  rs_lhs   <- sigma_qq * sigma_pp - sigma_qp^2
  rs_bound <- (hbar / 2)^2
  rs_satisfied <- rs_lhs >= rs_bound - .Machine$double.eps * abs(rs_bound)

  if (!rs_satisfied) {
    warning(sprintf(
      "Robertson-Schrödinger inequality violated: %.6e < %.6e",
      rs_lhs, rs_bound
    ))
  }

  # Classical semi-axes from 1-sigma contour (de Gosson Fermi blob convention)
  Delta_q <- sqrt(2 * sigma_qq)
  Delta_p <- sqrt(2 * sigma_pp)

  # Zurek reciprocal scales — the irreducible quantum widths
  delta_q <- hbar / Delta_p
  delta_p <- hbar / Delta_q

  # Classical action in units of A_0 = h/2 = pi*hbar
  # A = pi * Delta_q * Delta_p (Fermi blob area)
  # A_0 = pi * hbar (ground state area, since Delta_q=Delta_p=sqrt(hbar) at n=0)
  A_over_A0 <- (Delta_q * Delta_p) / hbar

  list(
    Delta_q      = Delta_q,
    Delta_p      = Delta_p,
    delta_q      = delta_q,
    delta_p      = delta_p,
    A_over_A0    = A_over_A0,
    rs_bound     = rs_bound,
    rs_satisfied = rs_satisfied
  )
}

#' QHO covariance matrix for eigenstate n in natural units.
#'
#' For the harmonic oscillator eigenstate |n>, the exact covariance is:
#'   <q^2> = (2n+1)/2,  <p^2> = (2n+1)/2,  <qp+pq>/2 = 0
#'
#' @param n Quantum number (non-negative integer)
#' @param hbar Reduced Planck constant (default 1.0)
#' @return Output of robertson_schrodinger() for eigenstate n
qho_covariance <- function(n, hbar = 1.0) {
  alpha    <- 2*n + 1
  sigma_qq <- alpha * hbar / 2
  sigma_pp <- alpha * hbar / 2
  sigma_qp <- 0
  robertson_schrodinger(sigma_qq, sigma_pp, sigma_qp, hbar)
}

#' True QHO Wigner function W_n(q,p) in natural units.
#'
#' Uses GSL Laguerre polynomials for numerical stability at large n.
#' The formula is:
#'   W_n(q,p) = (-1)^n / pi * exp(-(q^2+p^2)) * L_n(2*(q^2+p^2))
#'
#' Valid in natural units where hbar=1, q_0=p_0=1.
#'
#' @param n Quantum number
#' @param q Position (scalar or vector, in units of q_0)
#' @param p Momentum (scalar or vector, in units of p_0)
#' @return Wigner function value(s)
qho_wigner <- function(n, q, p) {
  rho2 <- q^2 + p^2
  (-1)^n / pi * exp(-rho2) * laguerre_n(n, 0, 2 * rho2)
}

#' Squeezed Gaussian kernel G_delta_q(q,p).
#'
#' The minimum-uncertainty kernel consistent with the classical action A.
#' Squeezed in q (width delta_q) and extended in p (width Delta_p).
#' Symplectic area = h/2, saturating the Heisenberg limit.
#'
#' @param q Position grid values
#' @param p Momentum grid values
#' @param rs Output of robertson_schrodinger() or qho_covariance()
#' @return Kernel values G_delta_q(q,p)
squeezed_kernel_q <- function(q, p, rs) {
  (1/pi) * exp(-q^2 / rs$delta_q^2 - p^2 / rs$Delta_p^2)
}

#' Squeezed Gaussian kernel G_delta_p(q,p).
#'
#' Conjugate kernel: squeezed in p (width delta_p), extended in q (width Delta_q).
#' Symplectic area = h/2, saturating the Heisenberg limit.
#'
#' @param q Position grid values
#' @param p Momentum grid values
#' @param rs Output of robertson_schrodinger() or qho_covariance()
#' @return Kernel values G_delta_p(q,p)
squeezed_kernel_p <- function(q, p, rs) {
  (1/pi) * exp(-q^2 / rs$Delta_q^2 - p^2 / rs$delta_p^2)
}

#' 2D FFT convolution of Wigner function with squeezed kernel.
#'
#' Computes P_delta_q(q,p) = W * G_delta_q on a uniform grid.
#' Uses ifftshift on the kernel to respect FFT origin convention.
#'
#' Hudson's theorem guarantees the result is non-negative since the kernel
#' has symplectic area h/2. Any negative values in the output are numerical
#' noise from floating point arithmetic in the FFT, not physical negativity.
#' These are flagged if they exceed machine epsilon tolerance.
#'
#' @param W_mat 2D matrix of Wigner function values on (q_int x p_int) grid
#' @param K_mat 2D matrix of kernel values on same grid (unnormalized)
#' @param dq Grid spacing in q
#' @param dp Grid spacing in p
#' @return Named list with P_mat (2D convolution) and max_negative (diagnostic)
fft_convolve_2d <- function(W_mat, K_mat, dq, dp) {
  nq <- nrow(W_mat)
  np <- ncol(W_mat)

  # Normalize kernel to unit integral
  K_norm <- K_mat / (sum(K_mat) * dq * dp)

  # ifftshift: move kernel center to FFT origin [1,1]
  ifftshift2d <- function(m) {
    nr <- nrow(m); nc <- ncol(m)
    sr <- floor(nr/2); sc <- floor(nc/2)
    rbind(
      cbind(m[(sr+1):nr, (sc+1):nc], m[(sr+1):nr, 1:sc]),
      cbind(m[1:sr,      (sc+1):nc], m[1:sr,      1:sc])
    )
  }
  K_shift <- ifftshift2d(K_norm)

  # FFT convolution
  P_mat <- Re(fft(fft(W_mat) * fft(K_shift), inverse = TRUE)) / (nq * np)
  P_mat <- P_mat * dq * dp

  # Diagnose numerical negatives — Hudson guarantees none should exist
  # physically. Flag if any exceed machine epsilon tolerance.
  peak_val    <- max(abs(P_mat))
  tol <- max(abs(P_mat)) * sqrt(.Machine$double.eps) * sqrt(nq * np)
  max_negative <- min(P_mat)

  if (max_negative < -tol) {
    warning(sprintf(
      "fft_convolve_2d: negative values %.2e exceed tolerance %.2e. Check grid resolution.",
      max_negative, -tol
    ))
  }

  # Zero negatives below tolerance — these are floating point artifacts
  P_mat[P_mat < 0 & P_mat >= -tol] <- 0

  list(P_mat = P_mat, max_negative = max_negative, tolerance = tol)
}

#' Compute the 1D symplectic density P_delta_q(q) by 2D convolution and projection.
#'
#' Full pipeline:
#'   1. Build W_n on a wide integration grid (captures Gaussian tails)
#'   2. Build G_delta_q kernel from Robertson-Schrödinger widths
#'   3. Convolve via FFT
#'   4. Project to 1D by integrating over p
#'   5. Normalize exactly once on the integration grid
#'   6. Interpolate onto the display grid (no renormalization)
#'
#' @param n Quantum number
#' @param wigner_fn Function(n, q, p) returning Wigner values
#' @param kernel_fn Function(q, p, rs) returning kernel values
#' @param rs Output of robertson_schrodinger() or qho_covariance()
#' @param q_display Display grid (interpolation target)
#' @param hbar Reduced Planck constant (default 1.0)
#' @return Named list with P_q (density on display grid) and diagnostics
compute_symplectic_density <- function(n, wigner_fn, kernel_fn, rs,
                                       q_display, hbar = 1.0) {
  Delta_q <- rs$Delta_q

  # Wide integration grid — must capture Gaussian tails beyond display range
  integ_lim <- max(Delta_q * 2.0, 4.0 * sqrt(hbar))
  integ_res <- max(601, 20 * ceiling(rs$A_over_A0))
  if (integ_res %% 2 == 0) integ_res <- integ_res + 1  # force odd for exact center

  q_int <- seq(-integ_lim, integ_lim, length.out = integ_res)
  p_int <- seq(-integ_lim, integ_lim, length.out = integ_res)
  dq    <- diff(q_int)[1]
  dp    <- diff(p_int)[1]

  # Build Wigner matrix and verify normalization
  W_mat  <- outer(q_int, p_int, FUN = function(q, p) wigner_fn(n, q, p))
  w_norm <- sum(W_mat) * dq * dp
  if (abs(w_norm - 1) > 1e-3) {
    warning(sprintf(
      "n=%d: Wigner norm = %.6f (expected 1.0). Increase integ_lim or integ_res.",
      n, w_norm
    ))
  }

  # Build kernel matrix
  K_mat <- outer(q_int, p_int, FUN = function(q, p) kernel_fn(q, p, rs))

  # 2D FFT convolution
  conv  <- fft_convolve_2d(W_mat, K_mat, dq, dp)
  P_mat <- conv$P_mat

  # Project to 1D: integrate over p
  P_q_int <- rowSums(P_mat) * dp

  # Normalize exactly once on the integration grid
  norm_int <- sum(P_q_int) * dq
  if (abs(norm_int) < .Machine$double.eps^0.5) {
    stop(sprintf("n=%d: Symplectic density integrates to zero. Check inputs.", n))
  }
  P_q_int <- P_q_int / norm_int

  # Interpolate onto display grid — display is a view, not a renormalization
  P_q <- approx(q_int, P_q_int, xout = q_display, rule = 1)$y
  P_q[is.na(P_q)] <- 0

  list(
    P_q          = P_q,
    w_norm       = w_norm,
    max_negative = conv$max_negative,
    tolerance    = conv$tolerance,
    integ_lim    = integ_lim,
    integ_res    = integ_res
  )
}
