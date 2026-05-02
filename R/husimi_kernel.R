# ==============================================================================
# husimi_kernel.R
# Husimi Q-function: fixed coherent-state convolution kernel.
#
# This file provides the Husimi-specific kernel matrix builder used by the
# convolution pipeline. Everything else — the Wigner / semiclassical
# computation, the convolution pipeline, the cross-section extraction, the
# heatmap rendering — lives in kernel-agnostic files (math_tools.R,
# wigner_density.R, semiclassical_density.R, plot_tools.R) and is shared
# with the symplectic kernel.
#
# Reference: Husimi 1940 Proc. Phys.-Math. Soc. Japan 22, 264;
#            Takahashi & Saito PRL 55, 645 (1985);
#            Lee Phys. Rep. 259, 147 (1995).
# Author: Brian S. Mulloy
# ==============================================================================

library(here)
source(here("R", "math_tools.R"))

# ------------------------------------------------------------------------------
# HUSIMI KERNEL
# G_coh(q,p) = (1/pi) * exp(-q^2 - p^2) — fixed unit-width Gaussian.
# ------------------------------------------------------------------------------

#' Husimi coherent-state kernel evaluated at offset (q, p) from grid midpoint.
husimi_kernel <- function(q, p) {
  (1/pi) * exp(-q^2 - p^2)
}

#' Build the Husimi kernel matrix on a (q_grid, p_grid) integration grid.
#'
#' Centered on the grid midpoint so ifftshift in fft_convolve_2d places the
#' kernel peak at the FFT origin.
husimi_kernel_matrix <- function(q_grid, p_grid) {
  q_mid <- (min(q_grid) + max(q_grid)) / 2
  p_mid <- (min(p_grid) + max(p_grid)) / 2
  outer(q_grid, p_grid,
        FUN = function(q, p) husimi_kernel(q - q_mid, p - p_mid))
}
