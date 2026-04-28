# ==============================================================================
# semiclassical_state.R
# Per-state bundle for semiclassical phase-space rendering.
#
# Drop-in replacement for build_wigner_state(): returns a bundle with the
# same field names (q_int, p_int, dq_int, dp_int, W_matrix, W_cross,
# heatmap_dt, norm) so apply_kernel_cross_section() works without
# modification.
#
# In addition, the bundle carries:
#   wkb_density   — analytical 1D WKB density on q_display
#                   (Inf at turning points, 0 in forbidden region)
#
# This file is kernel-agnostic. The kernel-specific bits live in
# symplectic_kernel.R / husimi_kernel.R as before.
#
# Author: Brian S. Mulloy
# ==============================================================================

library(here)
library(data.table)
source(here("R", "math_tools.R"))
source(here("R", "semiclassical_density.R"))

# ------------------------------------------------------------------------------
# DEFAULT REGULARIZATION
# epsilon = REG_FRACTION * E_n. Small relative to energy scale, large
# enough to span several grid cells of the integration grid.
# ------------------------------------------------------------------------------

SEMICLASSICAL_REG_FRACTION <- 0.02

# ------------------------------------------------------------------------------
# BUILD SEMICLASSICAL STATE
# ------------------------------------------------------------------------------

#' Build a per-state semiclassical bundle.
#'
#' @param E_n Energy of the eigenstate.
#' @param V_fn Function(q) returning potential values.
#' @param q_lo,q_hi Display window in q.
#' @param p_lo,p_hi Display window in p.
#' @param q_display Position grid for cross-section / projection rendering.
#' @param epsilon Regularization width for 2D shell. If NULL, defaults to
#'                SEMICLASSICAL_REG_FRACTION * E_n.
#' @param nq_int,np_int Integration grid resolution. Defaults parallel
#'                Wigner pipeline (801 x 601).
#' @return List with fields matching build_wigner_state() bundle plus
#'         wkb_density (1D analytical caustic density on q_display).
build_semiclassical_state <- function(E_n, V_fn,
                                      q_lo, q_hi, p_lo, p_hi, q_display,
                                      epsilon=NULL,
                                      nq_int=801, np_int=601) {
  if (is.null(epsilon)) epsilon <- SEMICLASSICAL_REG_FRACTION * abs(E_n)

  q_int <- seq(q_lo, q_hi, length.out=nq_int)
  p_int <- seq(p_lo, p_hi, length.out=np_int)
  dq_int <- diff(q_int)[1]
  dp_int <- diff(p_int)[1]

  cat(sprintf("    Building semiclassical shell on %d x %d grid (epsilon=%.4f)...\n",
              nq_int, np_int, epsilon))

  W_mat <- semiclassical_shell_density(q_int, p_int, E_n, V_fn, epsilon)
  norm  <- sum(W_mat) * dq_int * dp_int
  cat(sprintf("    Shell norm: %.6f\n", norm))

  # Extract p=0 cross-section (interpolated to q_display)
  p0_idx  <- which.min(abs(p_int))
  W_cross_int <- W_mat[, p0_idx]
  W_cross <- approx(q_int, W_cross_int, xout=q_display, rule=1)$y
  W_cross[is.na(W_cross)] <- 0

  # 1D WKB caustic density on q_display
  wkb <- wkb_caustic_density(q_display, E_n, V_fn)

  # Heatmap data for ggplot
  heatmap_dt <- as.data.table(expand.grid(q=q_int, p=p_int))
  heatmap_dt[, w := as.vector(W_mat)]
  # Apply soft contrast curve: sqrt for visual emphasis of the shell band
  max_w <- max(heatmap_dt$w, na.rm=TRUE)
  if (max_w > 0) {
    heatmap_dt[, w_plot := sqrt(pmax(w/max_w, 0))]
  } else {
    heatmap_dt[, w_plot := 0]
  }

  list(
    q_int=q_int, p_int=p_int,
    dq_int=dq_int, dp_int=dp_int,
    W_matrix=W_mat,
    W_cross=W_cross,
    wkb_density=wkb,
    heatmap_dt=heatmap_dt,
    norm=norm,
    epsilon=epsilon
  )
}
