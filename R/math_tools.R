# ==============================================================================
# math_tools.R
# Pure numerical methods — no physics, no display.
# Used internally by wigner_tools.R and husimi_tools.R.
#
# Reference: Bracewell, The Fourier Transform and Its Applications (McGraw-Hill)
#            Press et al., Numerical Recipes 3rd ed.
# Author: Brian S. Mulloy
# ==============================================================================

# ------------------------------------------------------------------------------
# 2D FFT CONVOLUTION
# P = W * K via the convolution theorem.
# Reference: Bracewell Ch.3, Press et al. Ch.12
# ------------------------------------------------------------------------------

#' 2D FFT convolution P = W * K on a uniform grid.
#' @param W_mat Input matrix (nq x np)
#' @param K_mat Kernel matrix (nq x np, unnormalized)
#' @param dq Grid spacing in q
#' @param dp Grid spacing in p
#' @return Named list: P_mat, max_negative, tolerance
fft_convolve_2d <- function(W_mat, K_mat, dq, dp) {
  nq <- nrow(W_mat); np <- ncol(W_mat)

  # Normalize kernel to unit integral
  K_norm <- K_mat / (sum(K_mat) * dq * dp)

  # ifftshift: move kernel center to FFT origin [1,1]
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
# Used by both Wigner and Husimi pipelines.
# ------------------------------------------------------------------------------

#' Extract p=0 cross-section from 2D matrix and interpolate to display grid.
#' @param mat 2D matrix [nq x np]
#' @param q_int Integration grid in q
#' @param p_int Integration grid in p
#' @param q_display Display grid
#' @return Vector of values on q_display
extract_p0_cross_section <- function(mat, q_int, p_int, q_display) {
  p0_idx <- which.min(abs(p_int))
  cross  <- approx(q_int, mat[, p0_idx], xout=q_display, rule=1)$y
  cross[is.na(cross)] <- 0
  cross
}

# ------------------------------------------------------------------------------
# CANONICAL PIPELINE
# Computes W(q,0) and P_conv(q,0) from the same 2D grid.
# Called by wigner_tools.R and husimi_tools.R — not directly by plot files.
# ------------------------------------------------------------------------------

#' Compute raw and convolved p=0 cross-sections from the same 2D grid.
#' @param wigner_fn Function(q_grid, p_grid) -> W matrix [nq x np]
#' @param kernel_fn Function(q_grid, p_grid) -> K matrix [nq x np]
#' @param q_int Integration grid in q
#' @param p_int Integration grid in p
#' @param q_display Display grid
#' @return Named list: W_cross, P_cross, w_norm, max_negative, tolerance
compute_cross_sections <- function(wigner_fn, kernel_fn, q_int, p_int, q_display) {
  dq <- diff(q_int)[1]; dp <- diff(p_int)[1]
  nq <- length(q_int);  np <- length(p_int)

  cat("    Building W on", nq, "x", np, "grid...\n")
  W_mat  <- wigner_fn(q_int, p_int)
  w_norm <- sum(W_mat) * dq * dp
  if (abs(w_norm-1) > 1e-3) warning(sprintf("Wigner norm=%.6f (expected 1.0).", w_norm))
  cat(sprintf("    Wigner norm: %.6f\n", w_norm))

  cat("    Building K on same grid...\n")
  K_mat <- kernel_fn(q_int, p_int)

  cat("    Convolving via FFT...\n")
  conv  <- fft_convolve_2d(W_mat, K_mat, dq, dp)
  P_mat <- conv$P_mat

  cat(sprintf("    Hudson check: max_neg=%.2e tol=%.2e OK=%s\n",
              conv$max_negative, conv$tolerance,
              ifelse(conv$max_negative >= -conv$tolerance, "YES", "NO")))

  W_cross <- extract_p0_cross_section(W_mat, q_int, p_int, q_display)
  P_cross <- extract_p0_cross_section(P_mat, q_int, p_int, q_display)

  list(W_cross=W_cross, P_cross=P_cross,
       w_norm=w_norm, max_negative=conv$max_negative, tolerance=conv$tolerance)
}
