# ==============================================================================
# marginal_tools.R
# Helpers for the marginal-residual comparison: P_{delta q} vs |psi(q)|^2
# and P_{delta p} vs |psi-hat(p)|^2 across non-Gaussian states.
#
# Provides:
#   G_delta_p_kernel_matrix    momentum-squeezed kernel (polar-dual partner
#                              of G_delta_q_kernel_matrix in symplectic_kernel.R)
#   symplectic_q_marginal      integral of P_{delta q}(q,p) over p, on q grid
#                              (renamed/aliased symplectic_marginal_density)
#   symplectic_p_marginal      integral of P_{delta p}(q,p) over q, on p grid
#   exact_position_marginal    |psi(q)|^2 from psi on a position grid
#   exact_momentum_marginal    |psi-hat(p)|^2 via FFT of psi
#   compute_residuals          L^1, L^2, L^inf residuals between two arrays
#                              on a shared grid, normalized to total mass
#
# Convention matches the rest of the pipeline:
#   - hbar = 1
#   - Delta_q = sqrt(2 * <(q - <q>)^2>), Delta_p = sqrt(2 * <(p - <p>)^2>)
#   - delta_q = hbar / Delta_p, delta_p = hbar / Delta_q
#
# Author: Brian S. Mulloy
# ==============================================================================

library(here)
source(here("R", "math_tools.R"))
source(here("R", "symplectic_kernel.R"))
source(here("R", "wigner_density.R"))

# ------------------------------------------------------------------------------
# G_DELTA_P KERNEL MATRIX
#
# Polar-dual partner of G_delta_q. Squeezed in p (width delta_p), spans
# Delta_q in q. Area = pi * Delta_q * delta_p = h/2.
# ------------------------------------------------------------------------------

#' G_delta_p kernel evaluated at offset (q,p) from grid midpoint.
G_delta_p_kernel <- function(q, p, Delta_q, Delta_p, hbar=1.0) {
  w <- symplectic_kernel_widths(Delta_q, Delta_p, hbar=hbar)
  (1/pi) * exp(-(q/Delta_q)^2 - (p/w$delta_p)^2)
}

#' Build the G_delta_p kernel matrix on a (q_grid, p_grid) integration grid.
#'
#' Same calling convention as G_delta_q_kernel_matrix() — kernel centered
#' on grid midpoint for ifftshift compatibility in fft_convolve_2d.
G_delta_p_kernel_matrix <- function(q_grid, p_grid, Delta_q, Delta_p, hbar=1.0) {
  q_mid <- (min(q_grid) + max(q_grid)) / 2
  p_mid <- (min(p_grid) + max(p_grid)) / 2
  outer(q_grid, p_grid,
        FUN = function(q, p) G_delta_p_kernel(q - q_mid, p - p_mid,
                                              Delta_q, Delta_p, hbar=hbar))
}

# ------------------------------------------------------------------------------
# POSITION MARGINAL OF P_{delta q}
#
# rho_{delta q}(q) = integral (W * G_{delta q})(q, p) dp
#
# Re-exposes symplectic_marginal_density() under a name that pairs
# naturally with the momentum partner below.
# ------------------------------------------------------------------------------

#' Position marginal of P_{delta q}: integral over p, returned on q_display.
symplectic_q_marginal <- function(state, Delta_q, Delta_p, q_display, hbar=1.0) {
  kernel_fn <- function(qg, pg) {
    G_delta_q_kernel_matrix(qg, pg, Delta_q, Delta_p, hbar=hbar)
  }
  symplectic_marginal_density(state, kernel_fn, q_display)
}

# ------------------------------------------------------------------------------
# MOMENTUM MARGINAL OF P_{delta p}
#
# rho_{delta p}(p) = integral (W * G_{delta p})(q, p) dq
#
# Mirror of symplectic_q_marginal in the conjugate variable.
# ------------------------------------------------------------------------------

#' Momentum marginal of P_{delta p}: integral over q, returned on p_display.
symplectic_p_marginal <- function(state, Delta_q, Delta_p, p_display, hbar=1.0) {
  K_mat <- G_delta_p_kernel_matrix(state$q_int, state$p_int,
                                   Delta_q, Delta_p, hbar=hbar)
  conv  <- fft_convolve_2d(state$W_matrix, K_mat,
                           state$dq_int, state$dp_int)
  # Marginalize over q: sum across rows (q axis) at each p.
  rho_int <- colSums(conv$P_mat) * state$dq_int
  rho     <- approx(state$p_int, rho_int, xout=p_display, rule=1)$y
  rho[is.na(rho)] <- 0
  rho
}

# ------------------------------------------------------------------------------
# EXACT MARGINALS FROM psi
#
# Position: |psi(q)|^2 directly on q_display (interpolated from psi_q_grid).
# Momentum: |psi-hat(p)|^2 via FFT of psi, with psi-hat normalized so that
#           integral |psi-hat(p)|^2 dp = 1 with hbar = 1.
#
# psi-hat(p) = (1/sqrt(2 pi hbar)) * integral psi(q) exp(-i p q / hbar) dq
# ------------------------------------------------------------------------------

#' |psi(q)|^2 on a display grid, interpolated from psi sampled on psi_q_grid.
exact_position_marginal <- function(psi_vec, psi_q_grid, q_display) {
  dq    <- diff(psi_q_grid)[1]
  norm  <- sqrt(sum(abs(psi_vec)^2)*dq)
  psi_n <- psi_vec / norm
  prob  <- as.numeric(abs(psi_n)^2)
  rho_q <- approx(psi_q_grid, prob, xout=q_display, rule=1)$y
  rho_q[is.na(rho_q)] <- 0
  rho_q
}

#' |psi-hat(p)|^2 on a display grid, via FFT of psi sampled on psi_q_grid.
#'
#' Convention: psi-hat(p) = (1/sqrt(2 pi hbar)) * integral psi(q) e^{-i p q / hbar} dq
#' so integral |psi-hat(p)|^2 dp = 1.
exact_momentum_marginal <- function(psi_vec, psi_q_grid, p_display, hbar=1.0) {
  nq      <- length(psi_q_grid)
  dq      <- diff(psi_q_grid)[1]
  norm    <- sqrt(sum(abs(psi_vec)^2)*dq)
  psi_n   <- psi_vec / norm
  # Pad for fine p sampling.
  M       <- 2^ceiling(log2(max(8*nq, 4*length(p_display))))
  psi_pad <- complex(real=rep(0, M), imaginary=rep(0, M))
  # Center psi in the padded array so the FFT phase corresponds to a
  # symmetric-window transform: place psi_n at indices that put q=0 at
  # the natural FFT origin (index 1 in R / index 0 in standard FFT).
  q0_idx  <- which.min(abs(psi_q_grid))
  # Standard layout: psi(q_k) at FFT index ((k - q0_idx) mod M) + 1
  shift   <- (seq_len(nq) - q0_idx) %% M + 1
  psi_pad[shift] <- psi_n
  # FFT: sum_k psi(q_k) e^{-2 pi i j k / M} for j = 0..M-1.
  fft_out <- fft(psi_pad)
  # Map to physical p: the integral is sum_k psi(q_k) e^{-i p q_k / hbar} dq.
  # With q_k = (k - q0_idx) * dq and j-th frequency in fft, the relation is
  #   p_j = 2 pi hbar j / (M dq) for j = 0..M/2,
  #   p_j = 2 pi hbar (j - M) / (M dq) for j > M/2 (negative p, wrapping).
  # The factor (dq / sqrt(2 pi hbar)) converts the discrete sum to the
  # continuous Fourier convention above.
  psi_hat_full <- (dq / sqrt(2*pi*hbar)) * fft_out
  half     <- M/2 + 1
  p_pos    <- 2*pi*hbar * (0:(half-1)) / (M * dq)
  p_neg    <- 2*pi*hbar * ((half:(M-1)) - M) / (M * dq)
  p_full   <- c(p_pos, p_neg)
  ord      <- order(p_full)
  p_sorted <- p_full[ord]
  psi_hat_sorted <- psi_hat_full[ord]
  prob_p   <- as.numeric(abs(psi_hat_sorted)^2)
  rho_p    <- approx(p_sorted, prob_p, xout=p_display, rule=1)$y
  rho_p[is.na(rho_p)] <- 0
  rho_p
}

# ------------------------------------------------------------------------------
# RESIDUALS
#
# Given two non-negative arrays on a shared grid, compute L^1, L^2, L^inf
# norms of their difference, expressed as percentages of the total mass
# of the reference (exact) array.
#
# The L^1 norm normalized to total mass is the most interpretable single
# number: "the convolved marginal differs from the exact marginal by X%
# of the state's total probability."
# ------------------------------------------------------------------------------

#' Residuals of approx vs. exact, on a shared uniform grid.
#'
#' Sign convention: residual = exact - approx. Positive residual at a
#' grid point means the model under-estimates there (typical at peaks
#' since the model is a smoothed version of the exact); negative residual
#' means the model over-estimates (typical in the valleys/tails into
#' which the smoothing redistributes mass).
#'
#' @param exact   Reference array (e.g. |psi(q)|^2 on q_display).
#' @param approx  Test array on the same grid (e.g. P_{delta q} marginal).
#' @param d       Grid spacing.
#' @return list with l1, l2, linf (each as a fraction; multiply by 100 for %),
#'         plus residual = (exact - approx) for plotting.
compute_residuals <- function(exact, approx, d) {
  total_mass <- sum(exact) * d
  if (total_mass <= 0) total_mass <- 1
  diff_arr   <- exact - approx
  l1         <- sum(abs(diff_arr)) * d / total_mass
  l2         <- sqrt(sum(diff_arr^2) * d) / total_mass
  linf       <- max(abs(diff_arr)) / max(exact)
  list(l1 = l1, l2 = l2, linf = linf,
       residual = diff_arr,
       exact_mass = total_mass)
}
