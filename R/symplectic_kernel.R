# ==============================================================================
# symplectic_kernel.R
# Action-scaled squeezed Gaussian kernels and the conjugate-blob overlay
# expressing the quantum of action.
#
# The symplectic resolution uses two conjugate squeezed kernels (Mulloy 2025
# manuscript, Sec. "Resolution via Action Capacity"):
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
# This file provides only the symplectic-specific pieces:
#   - kernel-width computation from RS / orbit covariance widths
#   - kernel matrix builder for G_delta_q (used in convolution)
#   - overlay-layers builder for the three QoA ellipses
#   - symplectic_marginal_density: convenience function for the
#     semiclassical right-column rho_{delta q}(q) marginal
#
# Reference: de Gosson, Symplectic Methods in Harmonic Analysis (Birkhauser
#            2011); Zurek Nature 412, 712 (2001); Mulloy 2025 (this work).
# Author: Brian S. Mulloy
# ==============================================================================

library(here)
library(ggplot2)
library(ggforce)
source(here("R", "math_tools.R"))

# ------------------------------------------------------------------------------
# KERNEL WIDTHS FROM COVARIANCE WIDTHS
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
#
# Used in the right-column convolution: P_delta_q = W * G_delta_q.
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
# Up to three nested ellipses: outer A, inner a_q, inner a_p.
# All centered at q_center. The `cells` parameter selects which subset
# to render; default is all three.
# ------------------------------------------------------------------------------

#' Symplectic overlay layers for a phase-space heatmap.
#'
#' Returns a list of ggplot layers ready to be added to a heatmap plot.
#' Same signature pattern as husimi_overlay_layers().
#'
#' Drawing order: outer A first, then inner cells so they appear on top
#' of the outer.
#'
#' Visual hierarchy: outer Fermi blob A is drawn with a slightly thicker
#' line than the inner conjugate quantum blobs a_q, a_p, to emphasize A
#' as the primary kinematic envelope. All cells are solid black.
#'
#' @param Delta_q,Delta_p Covariance widths (RS or orbit).
#' @param q_center Center of the overlay in q (orbit center).
#' @param hbar Planck constant.
#' @param cells Character vector specifying which cells to render.
#'   Any subset of c("A", "a_q", "a_p"). Order in this vector does not
#'   matter; cells are always drawn in canonical order (A first, then
#'   a_q, then a_p) so inner cells appear on top of the outer envelope.
#'   Defaults to all three.
#' @return A list of ggplot layers, one per requested cell.
symplectic_overlay_layers <- function(Delta_q, Delta_p, q_center=0, hbar=1.0,
                                      cells=c("A", "a_q", "a_p")) {
  valid_cells <- c("A", "a_q", "a_p")
  unknown <- setdiff(cells, valid_cells)
  if (length(unknown) > 0) {
    stop(sprintf("Unknown cell(s): %s. Valid: %s",
                 paste(unknown, collapse=", "),
                 paste(valid_cells, collapse=", ")))
  }

  w <- symplectic_kernel_widths(Delta_q, Delta_p, hbar=hbar)
  # Build ellipse paths directly with geom_path. Older ggforce versions on
  # some R installs drop linewidth on geom_ellipse and warn about it; we
  # construct the parametric paths ourselves so the line weights are
  # respected exactly. Each path has its own group id so the cells are
  # rendered as separate closed curves.
  theta <- seq(0, 2*pi, length.out=361)
  ellipse_path <- function(a, b, group_id) {
    data.frame(q     = q_center + a * cos(theta),
               p     =            b * sin(theta),
               group = group_id)
  }

  layers <- list()
  if ("A" %in% cells) {
    path_A <- ellipse_path(Delta_q, Delta_p, "A")
    layers <- c(layers, list(
      geom_path(data=path_A, aes(x=q, y=p, group=group),
                inherit.aes=FALSE, color="black", linewidth=0.2)))
  }
  if ("a_q" %in% cells) {
    path_a_q <- ellipse_path(w$delta_q, Delta_p, "a_q")
    layers <- c(layers, list(
      geom_path(data=path_a_q, aes(x=q, y=p, group=group),
                inherit.aes=FALSE, color="black", linewidth=0.2)))
  }
  if ("a_p" %in% cells) {
    path_a_p <- ellipse_path(Delta_q, w$delta_p, "a_p")
    layers <- c(layers, list(
      geom_path(data=path_a_p, aes(x=q, y=p, group=group),
                inherit.aes=FALSE, color="black", linewidth=0.2)))
  }
  layers
}

# ------------------------------------------------------------------------------
# SYMPLECTIC MARGINAL DENSITY
#
# Used in the right column of semiclassical figures:
#   rho_{delta q}(q) = integral (W_cl * G_{delta q})(q, p) dp
#
# Convolves the regularized energy shell with the symplectic kernel and
# marginalizes over p. The kernel widths are state-specific and supplied
# via closure (kernel_fn).
# ------------------------------------------------------------------------------

#' Apply the symplectic kernel to a state's W_matrix and return the 1D
#' marginal rho(q) on a display grid.
#'
#' @param state State bundle from build_semiclassical_state() or
#'              build_wigner_state() — must have q_int, p_int, dq_int,
#'              dp_int, W_matrix.
#' @param kernel_fn Closure with signature kernel_fn(q_grid, p_grid),
#'                  typically built via:
#'                    function(q, p) G_delta_q_kernel_matrix(q, p,
#'                                                           Delta_q, Delta_p)
#' @param q_display Display grid in q for the output.
#' @return Numeric vector of rho values on q_display, with NAs replaced by 0.
symplectic_marginal_density <- function(state, kernel_fn, q_display) {
  K_mat   <- kernel_fn(state$q_int, state$p_int)
  conv    <- fft_convolve_2d(state$W_matrix, K_mat,
                             state$dq_int, state$dp_int)
  rho_int <- rowSums(conv$P_mat) * state$dp_int
  rho     <- approx(state$q_int, rho_int, xout=q_display, rule=1)$y
  rho[is.na(rho)] <- 0
  rho
}
