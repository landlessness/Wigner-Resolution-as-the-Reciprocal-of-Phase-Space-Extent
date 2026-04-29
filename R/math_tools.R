# ==============================================================================
# math_tools.R
# Pure numerical primitives — no physics, no display, no orchestration.
#
# Exports:
#   fft_convolve_2d                 2D convolution via the convolution theorem
#   extract_p0_cross_section        interpolate p=0 slice of a 2D matrix to a
#                                   display grid in q
#   extract_diagonal_cross_section  interpolate the SW-NE diagonal (q = p)
#                                   slice of a 2D matrix to an arc-length
#                                   display grid (bilinear interpolation)
#
# Reference: Bracewell, The Fourier Transform and Its Applications
#            (McGraw-Hill); Press et al., Numerical Recipes 3rd ed.
# Author: Brian S. Mulloy
# ==============================================================================

# ------------------------------------------------------------------------------
# 2D FFT CONVOLUTION
# P = W * K via the convolution theorem.
# Reference: Bracewell Ch.3, Press et al. Ch.12
# ------------------------------------------------------------------------------

#' 2D FFT convolution P = W * K on a uniform grid.
#'
#' @param W_mat Input matrix (nq x np)
#' @param K_mat Kernel matrix (nq x np), evaluated relative to the grid
#'   midpoint so that ifftshift places the kernel peak at the FFT origin.
#' @param dq Grid spacing in q
#' @param dp Grid spacing in p
#' @return Named list: P_mat, max_negative, tolerance.
fft_convolve_2d <- function(W_mat, K_mat, dq, dp) {
  nq <- nrow(W_mat); np <- ncol(W_mat)

  # Normalize kernel to unit integral.
  K_norm <- K_mat / (sum(K_mat) * dq * dp)

  # ifftshift: move kernel center to FFT origin [1,1].
  ifftshift2d <- function(m) {
    nr <- nrow(m); nc <- ncol(m)
    sr <- floor(nr/2); sc <- floor(nc/2)
    rbind(cbind(m[(sr+1):nr,(sc+1):nc], m[(sr+1):nr,1:sc]),
          cbind(m[1:sr,    (sc+1):nc],  m[1:sr,    1:sc]))
  }

  K_shift      <- ifftshift2d(K_norm)
  P_mat        <- Re(fft(fft(W_mat)*fft(K_shift), inverse=TRUE)) / (nq*np) * dq*dp
  peak_val     <- max(abs(P_mat))
  tol          <- peak_val * sqrt(.Machine$double.eps) * sqrt(nq*np)
  max_negative <- min(P_mat)

  if (max_negative < -tol) warning(sprintf(
    "fft_convolve_2d: min=%.2e exceeds tol %.2e.", max_negative, -tol))

  P_mat[P_mat < 0 & P_mat >= -tol] <- 0

  list(P_mat=P_mat, max_negative=max_negative, tolerance=tol)
}

# ------------------------------------------------------------------------------
# CROSS-SECTION EXTRACTION
# Extracts p=0 slice from a 2D matrix and interpolates onto display grid.
# ------------------------------------------------------------------------------

#' Extract the p=0 slice of a 2D matrix and interpolate to a display grid.
#'
#' @param mat        2D matrix [nq x np]
#' @param q_int      Integration grid in q
#' @param p_int      Integration grid in p
#' @param q_display  Display grid in q
#' @return Vector of values on q_display.
extract_p0_cross_section <- function(mat, q_int, p_int, q_display) {
  p0_idx <- which.min(abs(p_int))
  cross  <- approx(q_int, mat[, p0_idx], xout=q_display, rule=1)$y
  cross[is.na(cross)] <- 0
  cross
}

# ------------------------------------------------------------------------------
# DIAGONAL CROSS-SECTION EXTRACTION
# Extracts a slice along the SW-NE diagonal (q = p) of a 2D matrix and
# interpolates onto an arc-length display grid.
#
# Parametrization: (q, p) = (s/sqrt(2), s/sqrt(2)) where s is arc length
# along the diagonal. Endpoints s = +/- sqrt(2) * max(|q|, |p|).
#
# Bilinear interpolation over the (q_int, p_int) grid: for each target s,
# compute (q_s, p_s), find the four enclosing grid points, and interpolate.
# Out-of-grid points return zero.
# ------------------------------------------------------------------------------

#' Extract the SW-NE diagonal slice of a 2D matrix and interpolate onto an
#' arc-length grid.
#'
#' @param mat        2D matrix [nq x np]
#' @param q_int      Integration grid in q (uniform, sorted)
#' @param p_int      Integration grid in p (uniform, sorted)
#' @param s_display  Arc-length display grid along the diagonal
#' @return Vector of values on s_display.
extract_diagonal_cross_section <- function(mat, q_int, p_int, s_display) {
  q_vals <- s_display / sqrt(2)
  p_vals <- s_display / sqrt(2)

  dq <- diff(q_int)[1]
  dp <- diff(p_int)[1]
  q_min <- q_int[1]; q_max <- q_int[length(q_int)]
  p_min <- p_int[1]; p_max <- p_int[length(p_int)]

  # Bilinear interpolation per target point.
  cross <- numeric(length(s_display))
  for (k in seq_along(s_display)) {
    qk <- q_vals[k]; pk <- p_vals[k]
    if (qk < q_min || qk > q_max || pk < p_min || pk > p_max) {
      cross[k] <- 0
      next
    }
    # Real-valued grid indices; floor/ceiling brackets.
    iq_real <- (qk - q_min) / dq + 1
    ip_real <- (pk - p_min) / dp + 1
    iq_lo <- floor(iq_real); iq_hi <- min(iq_lo + 1, length(q_int))
    ip_lo <- floor(ip_real); ip_hi <- min(ip_lo + 1, length(p_int))
    fq <- iq_real - iq_lo
    fp <- ip_real - ip_lo

    cross[k] <- (1-fq)*(1-fp)*mat[iq_lo, ip_lo] +
                   fq *(1-fp)*mat[iq_hi, ip_lo] +
                (1-fq)*   fp *mat[iq_lo, ip_hi] +
                   fq *   fp *mat[iq_hi, ip_hi]
  }
  cross
}
