# ==============================================================================
# wigner_tools.R
# Wigner function computation primitives.
#
# Two main functions:
#   solve_schrodinger() — finite-difference matrix diagonalization
#   wigner_fft()        — Wigner function via FFT method
#                         (Leonhardt 1997 Ch.5; Johansson et al. 2012)
#
# This file computes raw Wigner data on integration grids. Higher-level
# orchestration (extracting cross-sections, building heatmap data, applying
# kernels) lives in wigner_state.R.
#
# Units: positions in q_0, momenta in p_0, hbar = q_0*p_0 = 1
# Author: Brian S. Mulloy
# ==============================================================================

library(here)
library(gsl)
source(here("R", "math_tools.R"))

# Default integration grid sizes.
# For high quantum numbers the wavefunction oscillates fast; pass larger
# values via build_wigner_state() to avoid Nyquist aliasing.
WIGNER_DEFAULT_NQ <- 801
WIGNER_DEFAULT_NP <- 601

# ------------------------------------------------------------------------------
# QHO — exact analytic Wigner function (reference)
# ------------------------------------------------------------------------------

#' Exact QHO Wigner function W_n(q,p) = (-1)^n/pi * exp(-(q^2+p^2)) * L_n(2(q^2+p^2)).
qho_wigner <- function(n, q, p) {
  rho2 <- q^2 + p^2
  (-1)^n / pi * exp(-rho2) * laguerre_n(n, 0, 2*rho2)
}

# ------------------------------------------------------------------------------
# SCHRODINGER SOLVER
# H*psi = E*psi, H = -hbar^2/2 * d^2/dq^2 + V(q)
# ------------------------------------------------------------------------------

solve_schrodinger <- function(V_fn, q_min, q_max, dq=0.01, n_states=6, hbar=1.0) {
  q_grid <- seq(q_min, q_max, by=dq)
  nq     <- length(q_grid)
  V_vec  <- V_fn(q_grid)

  ke_diag <- hbar^2 / dq^2
  ke_off  <- -hbar^2 / (2*dq^2)
  diag_H  <- ke_diag + V_vec

  if (nq <= 3000) {
    H_mat <- diag(diag_H)
    for (j in seq_len(nq-1)) {
      H_mat[j,   j+1] <- ke_off
      H_mat[j+1, j  ] <- ke_off
    }
    eig      <- eigen(H_mat, symmetric=TRUE)
    energies <- rev(eig$values)
    psi_mat  <- eig$vectors[, rev(seq_len(ncol(eig$vectors)))]
  } else {
    if (!requireNamespace("RSpectra", quietly=TRUE)) stop("Install RSpectra")
    if (!requireNamespace("Matrix",   quietly=TRUE)) stop("Install Matrix")
    off_H  <- rep(ke_off, nq-1)
    H_sp   <- Matrix::bandSparse(nq, nq, k=c(-1,0,1),
                                 diagonals=list(off_H, diag_H, off_H))
    eig    <- RSpectra::eigs_sym(H_sp, k=n_states, which="SM")
    ord    <- order(eig$values)
    energies <- eig$values[ord]
    psi_mat  <- eig$vectors[, ord, drop=FALSE]
  }

  psi_mat <- psi_mat[, seq_len(n_states), drop=FALSE]
  for (j in seq_len(n_states)) {
    norm <- sqrt(sum(psi_mat[,j]^2) * dq)
    if (norm > 0) psi_mat[,j] <- psi_mat[,j] / norm
  }

  cat(sprintf("  Schrodinger solved: %d grid points, %d states\n", nq, n_states))
  for (j in seq_len(n_states)) {
    cat(sprintf("    n=%d: E=%.6f norm=%.6f\n",
                j-1, energies[j], sum(psi_mat[,j]^2)*dq))
  }

  list(energies=energies, psi_matrix=psi_mat, q_grid=q_grid, dq=dq)
}

# ------------------------------------------------------------------------------
# WIGNER FFT METHOD
#
# Definition: W(q, p) = (1/pi) * integral psi(q+x) psi(q-x) exp(2i*p*x) dx
#
# Discrete form: for each iq, define
#   rho[iq, k] = psi(q_iq + k*dq) * psi(q_iq - k*dq), k = 0..nq-1, hbar=1
# Then for real psi:
#   W(q_iq, p) = (dq/pi) * [ rho[iq,0] + 2 * sum_{k>=1} rho[iq,k] * cos(2*p*k*dq) ]
#              = (dq/pi) * ( 2 * Re[FFT(rho_padded)] - rho[iq,0] )
#
# The cosine sum over k is the real part of the FFT of rho along k. Padding
# rho with zeros before FFT gives fine sampling of the cosine kernel.
#
# Output frequencies: p_native[m] = pi * m / (M * dq), m = 0..M/2.
# We interpolate from p_native to the requested p_grid (Wigner is even in p
# for real psi, so |p| works).
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
  norm  <- sqrt(sum(psi_vec^2)*dq)
  if (norm > 0) psi_vec <- psi_vec/norm

  # Build rho[iq, k] = psi[iq+k] * psi[iq-k] for valid indices.
  # iq, k are 1-indexed in R; out-of-bounds entries are 0.
  iq_grid <- matrix(seq_len(nq), nrow=nq, ncol=nq)
  k_grid  <- matrix(seq_len(nq) - 1, nrow=nq, ncol=nq, byrow=TRUE)
  ip_idx  <- iq_grid + k_grid
  im_idx  <- iq_grid - k_grid
  valid   <- (ip_idx >= 1) & (ip_idx <= nq) & (im_idx >= 1) & (im_idx <= nq)
  safe_ip <- ifelse(valid, ip_idx, 1)
  safe_im <- ifelse(valid, im_idx, 1)
  rho     <- psi_vec[safe_ip] * psi_vec[safe_im]
  rho     <- ifelse(valid, rho, 0)
  dim(rho) <- c(nq, nq)

  # FFT along k. Pad to length M (power of 2) for fine p sampling.
  # M is chosen so the native p step pi/(M*dq) is at most ~1/4 of the
  # display p step; this controls interpolation error.
  M <- 2^ceiling(log2(max(2*nq, 4*np)))
  rho_padded <- matrix(0, nrow=nq, ncol=M)
  rho_padded[, 1:nq] <- rho

  # mvfft FFTs columns. We want FFT along k (rows of rho_padded) — transpose.
  fft_full   <- t(mvfft(t(rho_padded)))
  half       <- M/2 + 1
  W_native   <- (dq/pi) * (2*Re(fft_full[, 1:half, drop=FALSE]) - rho[, 1])
  p_native   <- pi * (0:(half-1)) / (M * dq)

  # Wigner is even in p for real psi: interpolate at |p|.
  abs_p <- abs(p_grid)
  W_mat <- matrix(0, nrow=nq, ncol=np)
  for (iq in seq_len(nq)) {
    interp <- approx(p_native, W_native[iq, ], xout=abs_p,
                     rule=2, ties=mean)$y
    interp[abs_p > p_native[half]] <- 0
    W_mat[iq, ] <- interp
  }
  W_mat
}
