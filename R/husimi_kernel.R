# ==============================================================================
# husimi_kernel.R
# Husimi-specific data: the coherent-state kernel and its visual overlay.
#
# This file only contains things unique to the Husimi method:
#   - the kernel matrix builder used by the convolution pipeline
#   - the overlay-layers builder used as visual annotation on a heatmap
#
# Everything else — the Wigner computation, the convolution pipeline, the
# cross-section extraction, the heatmap rendering — lives in kernel-agnostic
# files (wigner_tools.R, wigner_state.R, math_tools.R, plot_tools.R) and is
# shared with the symplectic kernel.
#
# Reference: Husimi 1940 Proc. Phys.-Math. Soc. Japan 22, 264;
#            Takahashi & Saito PRL 55, 645 (1985);
#            Lee Phys. Rep. 259, 147 (1995).
# Author: Brian S. Mulloy
# ==============================================================================

library(ggplot2)
library(ggforce)

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

# ------------------------------------------------------------------------------
# HUSIMI OVERLAY
# Single dashed circle of unit radius — the Husimi kernel's fixed contour.
# ------------------------------------------------------------------------------

#' Husimi overlay layers for a phase-space heatmap.
#'
#' Returns a list of ggplot layers ready to be added to a heatmap plot.
#' Same signature pattern as symplectic_overlay_layers().
#'
#' @param q_center Center of the circle in q (orbit center).
#' @return A list of ggplot layers (here: a single geom_circle).
husimi_overlay_layers <- function(q_center=0) {
  circle_data <- data.frame(x0=q_center, y0=0, r=1.0)
  list(
    geom_circle(data=circle_data, aes(x0=x0, y0=y0, r=r),
                inherit.aes=FALSE,
                color="gray40",
                linewidth=0.4,
                linetype="22")
  )
}
