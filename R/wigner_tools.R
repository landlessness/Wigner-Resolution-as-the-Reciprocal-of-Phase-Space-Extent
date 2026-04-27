# ==============================================================================
# wigner_tools.R
# Wigner function computation — exact (QHO) and numerical (anharmonic).
#
# Schrodinger solver: finite difference matrix diagonalization
# Reference: Numerov 1924, Press et al. Numerical Recipes 3rd ed. Ch.18
#            Landau, Paez & Bordeianu Computational Physics Ch.9
#
# Wigner FFT method: Leonhardt, Measuring the Quantum State of Light
#            (Cambridge 1997) Ch.5
#            Johansson, Nation & Nori, Comp. Phys. Comm. 183, 1760 (2012)
#            (QuTiP _wigner_fourier method)
#
# Units: positions in q_0, momenta in p_0, hbar = q_0*p_0 = 1
# Author: Brian S. Mulloy
# ==============================================================================

library(here)
library(gsl)
source(here("R", "math_tools.R"))

# Default integration grid sizes for Wigner/Husimi pipelines.
# For high quantum numbers the wavefunction oscillates fast; pass larger
# values via n_q_int / n_p_int to avoid Nyquist aliasing.
WIGNER_DEFAULT_NQ <- 801
WIGNER_DEFAULT_NP <- 601

# ------------------------------------------------------------------------------
# QHO — exact analytic Wigner function
# ------------------------------------------------------------------------------

#' Exact QHO Wigner function W_n(q,p).
qho_wigner <- function(n, q, p) {
  rho2 <- q^2 + p^2
  (-1)^n / pi * exp(-rho2) * laguerre_n(n, 0, 2*rho2)
}

# ------------------------------------------------------------------------------
# SCHRODINGER SOLVER
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
# Reference: Leonhardt 1997 Ch.5, Johansson et al. 2012
# ------------------------------------------------------------------------------

#' Wigner function via FFT method on sampled wavefunction.
wigner_fft <- function(psi_vec, q_grid, p_grid) {
  nq    <- length(q_grid)
  np    <- length(p_grid)
  dq    <- diff(q_grid)[1]
  norm  <- sqrt(sum(psi_vec^2)*dq)
  if (norm > 0) psi_vec <- psi_vec/norm

  W_mat <- matrix(0, nrow=nq, ncol=np)
  k_seq <- 0:(nq-1)

  for (iq in seq_len(nq)) {
    rho <- numeric(nq)
    for (k in k_seq) {
      ip <- iq+k; im <- iq-k
      pp <- if (ip >= 1 && ip <= nq) psi_vec[ip] else 0
      pm <- if (im >= 1 && im <= nq) psi_vec[im] else 0
      rho[k+1] <- pp * pm
    }
    for (ip_idx in seq_len(np)) {
      phase   <- exp(2i * p_grid[ip_idx] * k_seq * dq)
      contrib <- rho * phase
      W_mat[iq, ip_idx] <-
        (Re(contrib[1]) + 2*sum(Re(contrib[-1]))) * dq / pi
    }
  }
  W_mat
}

# ------------------------------------------------------------------------------
# WIGNER PIPELINE
# Returns W(q,0) cross-section on display grid.
# Integration grid extends beyond display to avoid boundary artifacts.
# n_q_int / n_p_int control resolution; raise for high quantum numbers.
# ------------------------------------------------------------------------------

#' Compute W(q,0) cross-section from sampled wavefunction.
#' @param n_q_int Number of q points in the integration grid (default 801).
#'   Must satisfy dq < 1 / (2 * highest spatial frequency in psi). For n=100
#'   harmonic states, set to ~2401 to avoid Nyquist aliasing.
#' @param n_p_int Number of p points (default 601). Generally less critical.
compute_wigner_cross_section <- function(psi_vec, psi_q_grid,
                                         q_lo, q_hi, p_lo, p_hi,
                                         q_display,
                                         n_q_int=WIGNER_DEFAULT_NQ,
                                         n_p_int=WIGNER_DEFAULT_NP) {
  disp_width <- q_hi - q_lo
  q_int <- seq(q_lo - disp_width, q_hi + disp_width, length.out=n_q_int)
  p_int <- seq(p_lo - 2,          p_hi + 2,          length.out=n_p_int)

  psi_int <- approx(psi_q_grid, psi_vec, xout=q_int,
                    rule=1, yleft=0, yright=0)$y

  cat(sprintf("    Computing Wigner on %d x %d grid (dq=%.4f)...\n",
              n_q_int, n_p_int, diff(q_int)[1]))
  W_mat  <- wigner_fft(psi_int, q_int, p_int)
  w_norm <- sum(W_mat) * diff(q_int)[1] * diff(p_int)[1]
  cat(sprintf("    Wigner norm: %.6f\n", w_norm))
  if (abs(w_norm-1) > 1e-3) warning(sprintf("Wigner norm=%.6f.", w_norm))

  extract_p0_cross_section(W_mat, q_int, p_int, q_display)
}
