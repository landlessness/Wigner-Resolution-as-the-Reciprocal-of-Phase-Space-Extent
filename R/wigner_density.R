# ==============================================================================
# wigner_density.R
# Wigner phase-space density W(q,p) and the per-state bundle that wraps it.
#
# Definition: W(q, p) = (1/pi*hbar) * integral psi(q+x) psi(q-x) exp(2i*p*x/hbar) dx
#
# Two layers in this file:
#
#   Low level — primitive computation:
#     wigner_fft(psi_vec, q_grid, p_grid)
#                           Wigner via FFT of psi(q+x)*conj(psi(q-x)) along x;
#                           Leonhardt 1997 Ch.5
#
#   High level — per-state bundle:
#     build_wigner_state(...)     orchestrates wigner_fft + cross-section
#                                 extraction + heatmap data into a single
#                                 bundle for plotting
#     apply_kernel_cross_section  convolves the Wigner with a kernel and
#                                 returns the p=0 cross-section on display grid
#
# The bundle returned by build_wigner_state() carries:
#   q_int, p_int          uniform integration grids
#   dq_int, dp_int        grid spacings
#   W_matrix              Wigner function on the integration grid
#   W_cross               W(q, 0) on the display grid
#   heatmap_dt            data.table of W clipped to display window, ready
#                         for ggplot's geom_raster
#   norm                  diagnostic: integrated Wigner norm (should be ~1)
#
# Both Husimi and symplectic Wigner figures call build_wigner_state() once
# per (state, display window) and then apply_kernel_cross_section() with
# their respective kernel-builder. This is the spine of the kernel-agnostic
# vs kernel-specific separation.
#
# Units: positions in q_0, momenta in p_0, hbar = q_0*p_0 = 1
# Reference: Wigner Phys. Rev. 40, 749 (1932)
#            Leonhardt, Measuring the Quantum State of Light Ch.5
# Author: Brian S. Mulloy
# ==============================================================================

library(here)
library(data.table)
source(here("R", "math_tools.R"))

# Default integration grid sizes.
# For high quantum numbers the wavefunction oscillates fast; pass larger
# values via build_wigner_state() to avoid Nyquist aliasing.
WIGNER_DEFAULT_NQ <- 801
WIGNER_DEFAULT_NP <- 601

# ------------------------------------------------------------------------------
# WIGNER FFT METHOD
#
# Discrete form: for each iq, define
#   rho[iq, k] = psi(q_iq + k*dq) * psi(q_iq - k*dq), k = 0..nq-1
# Then for real psi:
#   W(q_iq, p) = (dq/pi) * [ rho[iq,0] + 2 * sum_{k>=1} rho[iq,k] * cos(2*p*k*dq) ]
#              = (dq/pi) * ( 2 * Re[FFT(rho_padded)] - rho[iq,0] )
#
# Output frequencies: p_native[m] = pi * m / (M * dq), m = 0..M/2.
# Wigner is even in p for real psi, so we interpolate at |p|.
#
# Reference: Leonhardt 1997 Ch.5
# ------------------------------------------------------------------------------

#' Wigner function via FFT method, vectorized over k and over p.
#'
#' @param psi_vec Wavefunction sampled on q_grid (will be re-normalized).
#' @param q_grid Uniform position grid.
#' @param p_grid Momentum grid at which W is evaluated.
#' @return Matrix W[nq x np] with W[iq, ip] = W(q_grid[iq], p_grid[ip]).
wigner_fft <- function(psi_vec, q_grid, p_grid) {
  nq    <- length(q_grid)
  np    <- length(p_grid)
  dq    <- diff(q_grid)[1]
  # Use |psi|^2 so the norm is correct for complex psi as well as real.
  norm  <- sqrt(sum(abs(psi_vec)^2)*dq)
  if (norm > 0) psi_vec <- psi_vec/norm

  # Build rho[iq, k_idx] = psi(q_iq + k*dq) * conj(psi(q_iq - k*dq)) for
  # k = -(nq-1), ..., -1, 0, 1, ..., nq-1. Both signs of k must be filled
  # explicitly when psi is complex; the Hermitian symmetry rho(q, -k) =
  # conj(rho(q, k)) is what makes W real, but the FFT must see both sides.
  # For real psi, rho is real and symmetric in k, and this reduces to the
  # standard real-symmetric case.
  k_max  <- nq - 1
  k_vals <- -k_max:k_max                          # length 2*nq - 1
  nk     <- length(k_vals)
  iq_grid <- matrix(seq_len(nq), nrow=nq, ncol=nk)
  k_grid  <- matrix(k_vals,      nrow=nq, ncol=nk, byrow=TRUE)
  ip_idx  <- iq_grid + k_grid
  im_idx  <- iq_grid - k_grid
  valid   <- (ip_idx >= 1) & (ip_idx <= nq) & (im_idx >= 1) & (im_idx <= nq)
  safe_ip <- pmin(pmax(ip_idx, 1), nq)
  safe_im <- pmin(pmax(im_idx, 1), nq)
  rho     <- psi_vec[safe_ip] * Conj(psi_vec[safe_im])
  rho[!valid] <- 0
  dim(rho) <- c(nq, nk)

  # FFT along k. Pad to length M (power of 2) for fine p sampling.
  # Standard FFT layout: place rho[k=0] at index 1, rho[k=1] at index 2,
  # ..., rho[k=k_max] at index k_max+1, then zeros, then rho[k=-k_max]
  # at index M-k_max+1, ..., rho[k=-1] at index M.
  M <- 2^ceiling(log2(max(2*nq, 4*np)))
  rho_padded <- matrix(0, nrow=nq, ncol=M) + 0i
  # Positive-k slots (including k=0).
  rho_padded[, 1:(k_max+1)] <- rho[, (k_max+1):nk]   # k=0,1,...,k_max
  # Negative-k slots wrapped to the end.
  rho_padded[, (M-k_max+1):M] <- rho[, 1:k_max]      # k=-k_max,...,-1

  fft_full   <- t(mvfft(t(rho_padded)))
  # The Wigner is W(q, p) = (dq/pi) * sum_k rho(q, k dq) exp(-2 i p k dq).
  # The FFT result at index j (j=0,...,M-1) corresponds to frequency
  # 2*p*dq = -2*pi*j/M (note negative sign; mvfft uses e^{-2 pi i j k / M},
  # but the integrand in W has e^{-2 i p k dq}, so 2*p*dq = 2*pi*j/M
  # gives p = pi*j/(M*dq)). The Wigner is real by Hermitian symmetry of
  # rho, but for complex psi we still take Re() at the end to drop any
  # tiny imaginary residue from finite-grid quadrature.
  half       <- M/2 + 1
  W_native   <- (dq/pi) * Re(fft_full[, 1:half, drop=FALSE])
  p_native   <- pi * (0:(half-1)) / (M * dq)

  abs_p <- abs(p_grid)
  W_mat <- matrix(0, nrow=nq, ncol=np)
  # For real psi (W even in p) abs_p suffices; for complex psi (W not
  # symmetric in p), we need to evaluate W at both +p and -p. The full-
  # bilateral FFT above gives us the second half of frequencies as well
  # (indices half+1, ..., M correspond to negative p), so we can simply
  # interpolate from the full mvfft output rather than just the half-line.
  W_full <- (dq/pi) * Re(fft_full)
  # Frequencies for the full output: standard FFT layout maps index j to
  # p_j = pi*j/(M*dq) for j in 0..M/2, and to p_j = pi*(j-M)/(M*dq) for
  # j > M/2 (i.e., negative p, wrapping).
  p_full_pos <- pi * (0:(half-1)) / (M * dq)
  p_full_neg <- pi * ((half:(M-1)) - M) / (M * dq)
  p_full     <- c(p_full_pos, p_full_neg)
  ord        <- order(p_full)
  p_sorted   <- p_full[ord]

  for (iq in seq_len(nq)) {
    W_row    <- W_full[iq, ord]
    interp   <- approx(p_sorted, W_row, xout=p_grid,
                       rule=2, ties=mean)$y
    interp[p_grid < min(p_sorted) | p_grid > max(p_sorted)] <- 0
    W_mat[iq, ] <- interp
  }
  W_mat
}

# ------------------------------------------------------------------------------
# BUILD WIGNER STATE
#
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
#'         heatmap_dt, norm. Heatmap data uses the colour-ready name w_plot.
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
  # approx() does not natively handle complex y values; split into real
  # and imaginary parts when psi is complex, interpolate each, recombine.
  if (is.complex(psi_vec)) {
    psi_re   <- approx(psi_q_grid, Re(psi_vec), xout=q_int,
                       rule=1, yleft=0, yright=0)$y
    psi_im   <- approx(psi_q_grid, Im(psi_vec), xout=q_int,
                       rule=1, yleft=0, yright=0)$y
    psi_int  <- complex(real=psi_re, imaginary=psi_im)
  } else {
    psi_int  <- approx(psi_q_grid, psi_vec, xout=q_int,
                       rule=1, yleft=0, yright=0)$y
  }

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
#
# Takes a precomputed WignerState and a kernel-builder function, returns
# the convolved P(q, 0) cross-section on the display grid.
#
# The kernel-builder must have signature  kernel_fn(q_grid, p_grid) -> matrix
# evaluated relative to the grid midpoint (so that ifftshift places the
# kernel peak at the FFT origin). Both husimi_kernel_matrix and
# G_delta_q_kernel_matrix follow this convention.
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
