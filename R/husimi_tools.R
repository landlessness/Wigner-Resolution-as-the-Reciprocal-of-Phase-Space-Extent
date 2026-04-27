# ==============================================================================
# husimi_tools.R
# Husimi Q-function — fixed coherent-state kernel.
#
# Reference: Husimi 1940 Proc. Phys.-Math. Soc. Japan 22, 264
#            Takahashi & Saito PRL 55, 645 (1985)
#            Lee Phys. Rep. 259, 147 (1995)
#
# Units: positions in q_0, momenta in p_0, hbar = 1
# At the ground state A=A_0 the Husimi kernel coincides with our symplectic
# kernel. For A > A_0 our kernel sharpens anisotropically; Husimi does not.
# Author: Brian S. Mulloy
# ==============================================================================

library(here)
source(here("R", "math_tools.R"))

# Default heatmap resolution (per axis). Raise for high quantum numbers.
HUSIMI_DEFAULT_NHEAT <- 400

# ------------------------------------------------------------------------------
# HUSIMI KERNEL
# G_coh(q,p) = (1/pi)*exp(-q^2-p^2)
# Kernel must be evaluated relative to the integration grid's midpoint so
# that ifftshift will place the kernel peak at the FFT origin.
# ------------------------------------------------------------------------------

husimi_kernel <- function(q, p) {
  (1/pi) * exp(-q^2 - p^2)
}

husimi_ellipse_data <- function(q_center=0) {
  list(circle=data.frame(x0=q_center, y0=0, r=1.0))
}

#' Build the Husimi kernel on a grid, centered on the grid's midpoint.
husimi_kernel_matrix <- function(q_grid, p_grid) {
  q_mid <- (min(q_grid) + max(q_grid)) / 2
  p_mid <- (min(p_grid) + max(p_grid)) / 2
  outer(q_grid, p_grid,
        FUN=function(q,p) husimi_kernel(q - q_mid, p - p_mid))
}

# ------------------------------------------------------------------------------
# HUSIMI PIPELINE
# ------------------------------------------------------------------------------

#' Compute Q(q,0) Husimi cross-section from sampled wavefunction.
#' @param n_q_int / n_p_int Integration grid resolution (default 801 / 601).
compute_husimi_cross_section <- function(psi_vec, psi_q_grid,
                                         q_lo, q_hi, p_lo, p_hi,
                                         q_display, wigner_fft_fn,
                                         n_q_int=WIGNER_DEFAULT_NQ,
                                         n_p_int=WIGNER_DEFAULT_NP) {
  disp_width <- q_hi - q_lo
  q_int <- seq(q_lo - disp_width, q_hi + disp_width, length.out=n_q_int)
  p_int <- seq(p_lo - 2,          p_hi + 2,          length.out=n_p_int)

  psi_int <- approx(psi_q_grid, psi_vec, xout=q_int,
                    rule=1, yleft=0, yright=0)$y

  wigner_fn <- function(q_grid, p_grid) wigner_fft_fn(psi_int, q_int, p_grid)
  kernel_fn <- function(q_grid, p_grid) husimi_kernel_matrix(q_grid, p_grid)

  cat("    Computing Husimi cross-section...\n")
  result <- compute_cross_sections(wigner_fn, kernel_fn, q_int, p_int, q_display)
  result$P_cross
}

#' Compute 2D Husimi heatmap from sampled wavefunction.
#' @param n_heat Heatmap resolution per axis (default 400). Raise for high n.
compute_husimi_heatmap <- function(psi_vec, psi_q_grid,
                                   q_lo, q_hi, p_lo, p_hi,
                                   wigner_fft_fn,
                                   n_heat=HUSIMI_DEFAULT_NHEAT) {
  disp_width <- q_hi - q_lo
  q_heat <- seq(q_lo - disp_width, q_hi + disp_width, length.out=n_heat)
  p_heat <- seq(p_lo - 2,          p_hi + 2,          length.out=n_heat)

  psi_heat <- approx(psi_q_grid, psi_vec, xout=q_heat,
                     rule=1, yleft=0, yright=0)$y

  W_heat   <- wigner_fft_fn(psi_heat, q_heat, p_heat)
  K_heat   <- husimi_kernel_matrix(q_heat, p_heat)
  dq_h     <- diff(q_heat)[1]; dp_h <- diff(p_heat)[1]
  conv_h   <- fft_convolve_2d(W_heat, K_heat, dq_h, dp_h)

  q_idx <- q_heat >= q_lo & q_heat <= q_hi
  p_idx <- p_heat >= p_lo & p_heat <= p_hi

  dt <- as.data.table(expand.grid(q=q_heat[q_idx], p=p_heat[p_idx]))
  dt[, w := as.vector(conv_h$P_mat[q_idx, p_idx])]
  # Linear normalization to [-1, 1] for symmetric diverging colormap.
  # No power compression — keeps oscillation contrast visually honest.
  max_abs <- max(abs(dt$w), na.rm=TRUE)
  if (max_abs > 0) dt[, w_plot := w/max_abs] else dt[, w_plot := w]
  dt
}
