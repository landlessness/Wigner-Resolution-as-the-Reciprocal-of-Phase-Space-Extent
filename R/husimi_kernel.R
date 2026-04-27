# ==============================================================================
# husimi_kernel.R
# Husimi-specific data: the coherent-state kernel and its visual annotation.
#
# This file only contains things that are unique to the Husimi method:
#   - the kernel function itself (a fixed unit-width Gaussian)
#   - the kernel-matrix builder used by the convolution pipeline
#   - the ellipse annotation data used as overlay on phase-space plots
#
# Everything else — the Wigner computation, the convolution pipeline, the
# cross-section extraction — lives in kernel-agnostic files (wigner_tools.R,
# wigner_state.R, math_tools.R) and is shared with the symplectic kernel.
#
# Reference: Husimi 1940 Proc. Phys.-Math. Soc. Japan 22, 264
#            Takahashi & Saito PRL 55, 645 (1985)
#            Lee Phys. Rep. 259, 147 (1995)
#
# Author: Brian S. Mulloy
# ==============================================================================

# ------------------------------------------------------------------------------
# HUSIMI KERNEL
# G_coh(q,p) = (1/pi) * exp(-q^2 - p^2)
# A fixed unit-width Gaussian, independent of the system's classical action.
# ------------------------------------------------------------------------------

#' Husimi coherent-state kernel evaluated at (q, p).
husimi_kernel <- function(q, p) {
  (1/pi) * exp(-q^2 - p^2)
}

#' Build the Husimi kernel matrix on a (q_grid, p_grid) integration grid.
#'
#' The kernel is centered on the grid's midpoint so that ifftshift in
#' fft_convolve_2d places the kernel peak at the FFT origin.
#'
#' This function has the same signature as future symplectic_kernel_matrix(),
#' so the convolution pipeline can call either one interchangeably.
husimi_kernel_matrix <- function(q_grid, p_grid) {
  q_mid <- (min(q_grid) + max(q_grid)) / 2
  p_mid <- (min(p_grid) + max(p_grid)) / 2
  outer(q_grid, p_grid,
        FUN=function(q, p) husimi_kernel(q - q_mid, p - p_mid))
}

# ------------------------------------------------------------------------------
# HUSIMI ELLIPSE ANNOTATION
# Visual overlay for phase-space plots: a single dashed circle of radius 1
# centered at q_center, expressing the Husimi kernel's fixed coherent-state
# width.
# ------------------------------------------------------------------------------

#' Annotation data for the Husimi kernel overlay (a unit circle).
#' @param q_center Center of the circle in q (default 0). The kernel itself
#'   is translation-invariant during convolution; q_center is only the
#'   visual placement on the heatmap.
husimi_ellipse_data <- function(q_center=0) {
  list(circle=data.frame(x0=q_center, y0=0, r=1.0))
}
