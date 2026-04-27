# ==============================================================================
# wigner_state.R
# The per-state Wigner data bundle and the operations that build it.
#
# A WignerState is everything we know about a quantum state in phase space
# that does not depend on which kernel (Husimi, symplectic, ...) we plan
# to apply later:
#
#   q_int, p_int          uniform integration grids
#   dq_int, dp_int        grid spacings
#   W_matrix              Wigner function on the integration grid
#   W_cross               W(q, 0) on the display grid
#   heatmap_dt            data.table of W clipped to display window, ready
#                         for ggplot's geom_raster
#   norm                  diagnostic: integrated Wigner norm (should be ~1)
#
# Both Husimi and symplectic figures call build_wigner_state() once per
# (state, display window) and then call apply_kernel_cross_section() with
# their respective kernel-builder. This is the spine of the kernel-agnostic
# vs kernel-specific separation.
#
# Author: Brian S. Mulloy
# ==============================================================================

library(here)
library(data.table)
source(here("R", "math_tools.R"))
source(here("R", "wigner_tools.R"))

# ------------------------------------------------------------------------------
# BUILD WIGNER STATE
# Computes W on an integration grid that extends beyond the display window
# to suppress boundary artifacts. Returns the integration-grid W matrix
# along with the cross-section and heatmap data extracted from it.
# ------------------------------------------------------------------------------

#' Build the per-state Wigner data bundle.
#'
#' @param psi_vec      Wavefunction sampled on psi_q_grid.
#' @param psi_q_grid   The grid on which psi_vec is sampled.
#' @param q_lo,q_hi    Display window in q.
#' @param p_lo,p_hi    Display window in p.
#' @param q_display    The display q grid for cross-section output.
#' @param n_q_int      Integration grid size in q (default 801). Raise for
#'                     high-quantum-number states to avoid Nyquist aliasing.
#' @param n_p_int      Integration grid size in p (default 601).
#' @return A list with q_int, p_int, dq_int, dp_int, W_matrix, W_cross,
#'   heatmap_dt, norm. Heatmap data uses the colour-ready name `w_plot`.
build_wigner_state <- function(psi_vec, psi_q_grid,
                               q_lo, q_hi, p_lo, p_hi,
                               q_display,
                               n_q_int=WIGNER_DEFAULT_NQ,
                               n_p_int=WIGNER_DEFAULT_NP) {

  # Integration grid extends one display-width on each side in q to suppress
  # boundary artifacts; in p we extend by 2 units (a few coherent-state widths).
  disp_width <- q_hi - q_lo
  q_int <- seq(q_lo - disp_width, q_hi + disp_width, length.out=n_q_int)
  p_int <- seq(p_lo - 2,          p_hi + 2,          length.out=n_p_int)
  dq_int <- diff(q_int)[1]
  dp_int <- diff(p_int)[1]

  # Interpolate psi onto the integration q grid (zero outside support).
  psi_int <- approx(psi_q_grid, psi_vec, xout=q_int,
                    rule=1, yleft=0, yright=0)$y

  # Compute the Wigner matrix.
  cat(sprintf("    Building Wigner on %d x %d grid (dq=%.4f)...\n",
              n_q_int, n_p_int, dq_int))
  W_mat  <- wigner_fft(psi_int, q_int, p_int)
  w_norm <- sum(W_mat) * dq_int * dp_int
  cat(sprintf("    Wigner norm: %.6f\n", w_norm))
  if (abs(w_norm-1) > 1e-3) warning(sprintf("Wigner norm=%.6f.", w_norm))

  # Extract W(q, 0) cross-section on the display grid.
  W_cross <- extract_p0_cross_section(W_mat, q_int, p_int, q_display)

  # Build the heatmap data: clip W to display window and pack as data.table.
  q_mask <- q_int >= q_lo & q_int <= q_hi
  p_mask <- p_int >= p_lo & p_int <= p_hi
  heatmap_dt <- as.data.table(expand.grid(
    q = q_int[q_mask],
    p = p_int[p_mask]
  ))
  heatmap_dt[, w := as.vector(W_mat[q_mask, p_mask])]
  # Symmetric normalization for diverging colormap.
  max_abs <- max(abs(heatmap_dt$w), na.rm=TRUE)
  heatmap_dt[, w_plot := if (max_abs > 0) w/max_abs else w]

  list(
    q_int      = q_int,
    p_int      = p_int,
    dq_int     = dq_int,
    dp_int     = dp_int,
    W_matrix   = W_mat,
    W_cross    = W_cross,
    heatmap_dt = heatmap_dt,
    norm       = w_norm
  )
}

# ------------------------------------------------------------------------------
# APPLY KERNEL: COMPUTE CONVOLVED CROSS-SECTION
# Takes a precomputed WignerState and a kernel-builder function, returns
# the convolved P(q, 0) cross-section on the display grid.
#
# The kernel-builder must have signature  kernel_fn(q_grid, p_grid) -> matrix
# evaluated relative to the grid midpoint (so that ifftshift places the
# kernel peak at the FFT origin). Both husimi_kernel_matrix and the future
# symplectic_kernel_matrix follow this convention.
# ------------------------------------------------------------------------------

#' Apply a kernel to a precomputed WignerState and return P(q, 0).
#'
#' @param state      Output of build_wigner_state().
#' @param kernel_fn  Kernel-builder: function(q_grid, p_grid) -> matrix.
#' @param q_display  Display grid in q for the output cross-section.
#' @return P_cross vector on q_display.
apply_kernel_cross_section <- function(state, kernel_fn, q_display) {
  K_mat <- kernel_fn(state$q_int, state$p_int)
  cat("    Convolving via FFT...\n")
  conv  <- fft_convolve_2d(state$W_matrix, K_mat,
                           state$dq_int, state$dp_int)
  cat(sprintf("    Hudson check: max_neg=%.2e tol=%.2e OK=%s\n",
              conv$max_negative, conv$tolerance,
              ifelse(conv$max_negative >= -conv$tolerance, "YES", "NO")))
  extract_p0_cross_section(conv$P_mat, state$q_int, state$p_int, q_display)
}
