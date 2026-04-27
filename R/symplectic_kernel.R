# ==============================================================================
# symplectic_kernel.R
# Symplectic-specific data: action-scaled squeezed Gaussian kernels and the
# three-ellipse visual overlay expressing the quantum of action.
#
# The symplectic resolution uses two conjugate squeezed kernels (Mulloy 2025
# manuscript, Sec. "Resolution via Finite Classical Action"):
#
#   G_delta_q(q,p) = (1/pi) * exp(-q^2/delta_q^2 - p^2/Delta_p^2)
#                    squeezed in q (width delta_q), spans Delta_p in p
#                    area = pi * delta_q * Delta_p = h/2 (= a_q)
#
#   G_delta_p(q,p) = (1/pi) * exp(-q^2/Delta_q^2 - p^2/delta_p^2)
#                    squeezed in p (width delta_p), spans Delta_q in q
#                    area = pi * Delta_q * delta_p = h/2 (= a_p)
#
# where delta_q = hbar/Delta_p and delta_p = hbar/Delta_q. Both saturate
# the Heisenberg limit; their areas sum to the quantum of action h.
#
# Tomographic complementarity: a_q and a_p are accessible only on separate
# copies of the state — never simultaneously on a single copy. The overlay
# shows both contours plus the outer Fermi blob A.
#
# Geometry: a_q and a_p are always inscribed within A. a_q has the same
# Delta_p extent as A but is narrower in q; a_p has the same Delta_q extent
# as A but is narrower in p. The QoA's bounding box equals A's.
#
# This file provides only the symplectic-specific pieces:
#   - kernel-width computation from RS covariance widths
#   - kernel matrix builder for G_delta_q (used in cross-section convolution)
#   - overlay-layers builder for the three QoA ellipses
#
# Reference: de Gosson, Symplectic Methods in Harmonic Analysis (Birkhauser
#            2011); Zurek Nature 412, 712 (2001); Mulloy 2025 (this work).
# Author: Brian S. Mulloy
# ==============================================================================

library(ggplot2)
library(ggforce)

# ------------------------------------------------------------------------------
# KERNEL WIDTHS FROM RS COVARIANCE
# ------------------------------------------------------------------------------

#' Symplectic resolution scales (Zurek's reciprocal scales) from RS widths.
#'
#' @param Delta_q State width in q (= sqrt(2*sigma_qq)).
#' @param Delta_p State width in p (= sqrt(2*sigma_pp)).
#' @param hbar Planck constant in chosen units.
#' @return list(delta_q, delta_p) with delta_q = hbar/Delta_p,
#'         delta_p = hbar/Delta_q.
symplectic_kernel_widths <- function(Delta_q, Delta_p, hbar=1.0) {
  list(delta_q = hbar/Delta_p,
       delta_p = hbar/Delta_q)
}

# ------------------------------------------------------------------------------
# G_DELTA_Q KERNEL MATRIX
# Used in the right-column cross-section convolution: P_delta_q = W * G_delta_q.
# Built relative to the integration grid's midpoint so ifftshift in
# fft_convolve_2d places the kernel peak at the FFT origin.
# ------------------------------------------------------------------------------

#' G_delta_q kernel evaluated at offset (q,p) from grid midpoint.
G_delta_q_kernel <- function(q, p, Delta_q, Delta_p, hbar=1.0) {
  w <- symplectic_kernel_widths(Delta_q, Delta_p, hbar=hbar)
  (1/pi) * exp(-(q/w$delta_q)^2 - (p/Delta_p)^2)
}

#' Build the G_delta_q kernel matrix on a (q_grid, p_grid) integration grid.
#'
#' Same calling convention as husimi_kernel_matrix() — kernel is centered on
#' the grid midpoint. State-specific widths are supplied via closure when
#' this function is passed to apply_kernel_cross_section().
G_delta_q_kernel_matrix <- function(q_grid, p_grid, Delta_q, Delta_p, hbar=1.0) {
  q_mid <- (min(q_grid) + max(q_grid)) / 2
  p_mid <- (min(p_grid) + max(p_grid)) / 2
  outer(q_grid, p_grid,
        FUN = function(q, p) G_delta_q_kernel(q - q_mid, p - p_mid,
                                              Delta_q, Delta_p, hbar=hbar))
}

# ------------------------------------------------------------------------------
# SYMPLECTIC OVERLAY
# Three nested ellipses: outer A (solid), inner a_q (dashed), inner a_p
# (dashed). All centered at q_center.
# ------------------------------------------------------------------------------

#' Symplectic overlay layers for a phase-space heatmap.
#'
#' Returns a list of ggplot layers ready to be added to a heatmap plot.
#' Same signature pattern as husimi_overlay_layers().
#'
#' Drawing order: outer A first, then inner a_q and a_p so they appear on
#' top of the outer.
#'
#' @param Delta_q,Delta_p RS-covariance widths.
#' @param q_center Center of the overlay in q (orbit center).
#' @param hbar Planck constant.
#' @return A list of three ggplot layers (one per ellipse).
symplectic_overlay_layers <- function(Delta_q, Delta_p, q_center=0, hbar=1.0) {
  w <- symplectic_kernel_widths(Delta_q, Delta_p, hbar=hbar)
  ellipse_A   <- data.frame(x0=q_center, y0=0, a=Delta_q,   b=Delta_p,   angle=0)
  ellipse_a_q <- data.frame(x0=q_center, y0=0, a=w$delta_q, b=Delta_p,   angle=0)
  ellipse_a_p <- data.frame(x0=q_center, y0=0, a=Delta_q,   b=w$delta_p, angle=0)
  list(
    geom_ellipse(data=ellipse_A,
                 aes(x0=x0, y0=y0, a=a, b=b, angle=angle),
                 inherit.aes=FALSE,
                 color="gray80", linewidth=0.5, linetype="22"),
    geom_ellipse(data=ellipse_a_q,
                 aes(x0=x0, y0=y0, a=a, b=b, angle=angle),
                 inherit.aes=FALSE,
                 color="gray40", linewidth=0.4, linetype="22"),
    geom_ellipse(data=ellipse_a_p,
                 aes(x0=x0, y0=y0, a=a, b=b, angle=angle),
                 inherit.aes=FALSE,
                 color="gray40", linewidth=0.4, linetype="22")
  )
}
