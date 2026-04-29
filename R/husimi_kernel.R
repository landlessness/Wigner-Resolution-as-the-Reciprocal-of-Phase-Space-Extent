# ==============================================================================
# husimi_kernel.R
# Husimi Q-function: fixed coherent-state convolution kernel.
#
# This file provides only the Husimi-specific pieces:
#   - the kernel matrix builder used by the convolution pipeline
#   - the overlay-layers builder used as visual annotation on a heatmap
#
# Everything else — the Wigner / semiclassical computation, the convolution
# pipeline, the cross-section extraction, the heatmap rendering — lives in
# kernel-agnostic files (math_tools.R, wigner_density.R,
# semiclassical_density.R, plot_tools.R) and is shared with the
# symplectic kernel.
#
# Reference: Husimi 1940 Proc. Phys.-Math. Soc. Japan 22, 264;
#            Takahashi & Saito PRL 55, 645 (1985);
#            Lee Phys. Rep. 259, 147 (1995).
# Author: Brian S. Mulloy
# ==============================================================================

library(here)
library(ggplot2)
library(ggforce)
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
                color="black",
                linewidth=0.25,
                linetype="solid")
  )
}

# ------------------------------------------------------------------------------
# HUSIMI MARGINAL DENSITY
#
# Used in the right column of semiclassical figures:
#   rho_Q(q) = integral (W_cl * G_husimi)(q, p) dp
#
# Mirror of symplectic_marginal_density() but with the fixed unit-width
# Husimi kernel — no state-specific widths to thread through.
# ------------------------------------------------------------------------------

#' Apply the Husimi kernel to a state's W_matrix and return the 1D
#' marginal rho_Q(q) on a display grid.
#'
#' @param state State bundle from build_semiclassical_state() or
#'              build_wigner_state() — must have q_int, p_int, dq_int,
#'              dp_int, W_matrix.
#' @param q_display Display grid in q for the output.
#' @return Numeric vector of rho values on q_display.
husimi_marginal_density <- function(state, q_display) {
  K_mat   <- husimi_kernel_matrix(state$q_int, state$p_int)
  conv    <- fft_convolve_2d(state$W_matrix, K_mat,
                             state$dq_int, state$dp_int)
  rho_int <- rowSums(conv$P_mat) * state$dp_int
  rho     <- approx(state$q_int, rho_int, xout=q_display, rule=1)$y
  rho[is.na(rho)] <- 0
  rho
}
