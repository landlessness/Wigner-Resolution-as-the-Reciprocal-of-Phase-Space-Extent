# ==============================================================================
# schroedinger_solver.R
# Time-independent Schroedinger eigenvalue problem on a uniform q grid.
#
# H psi = E psi, with H = -hbar^2/2 * d^2/dq^2 + V(q)
#
# Finite-difference Hamiltonian. Dense eigen() for small grids; sparse
# RSpectra::eigs_sym for large grids. Re-orthonormalizes the returned
# eigenvectors as a numerical convenience.
#
# This solver is shared by all pipelines that need a Schroedinger eigenstate
# psi_n: the Wigner densities and (in the wider trilogy) Schroedinger-truth
# reference plots. The semiclassical pipeline does NOT use this solver — it
# computes the energy shell from V(q) and E directly, and the orbit-derived
# covariance from V(q) and E directly (see classical_action_tools.R).
#
# Reference: Press et al., Numerical Recipes 3rd ed. Ch.18 (eigenvalue
#            problems on uniform grids); Lehoucq, Sorensen & Yang ARPACK
#            Users' Guide (SIAM 1998) for the sparse path.
# Author: Brian S. Mulloy
# ==============================================================================

#' Solve the time-independent Schroedinger equation on a uniform q grid.
#'
#' @param V_fn Function(q) returning potential values.
#' @param q_min Left edge of grid.
#' @param q_max Right edge of grid.
#' @param dq Grid spacing.
#' @param n_states Number of low-lying eigenstates to return.
#' @param hbar Planck constant in chosen units.
#' @return Named list: energies (length n_states), psi_matrix (nq x n_states,
#'         each column an L2-normalized eigenstate), q_grid, dq.
solve_schroedinger <- function(V_fn, q_min, q_max, dq=0.01,
                               n_states=6, hbar=1.0) {
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
    eig    <- RSpectra::eigs_sym(H_sp, k=n_states, which="SM",
                                 opts=list(ncv=max(8*n_states, 100),
                                           maxitr=2000))
    ord    <- order(eig$values)
    energies <- eig$values[ord]
    psi_mat  <- eig$vectors[, ord, drop=FALSE]
  }

  psi_mat <- psi_mat[, seq_len(n_states), drop=FALSE]
  for (j in seq_len(n_states)) {
    norm <- sqrt(sum(psi_mat[,j]^2) * dq)
    if (norm > 0) psi_mat[,j] <- psi_mat[,j] / norm
  }

  cat(sprintf("  Schroedinger solved: %d grid points, %d states\n",
              nq, n_states))
  for (j in seq_len(n_states)) {
    cat(sprintf("    n=%d: E=%.6f norm=%.6f\n",
                j-1, energies[j], sum(psi_mat[,j]^2)*dq))
  }

  list(energies=energies, psi_matrix=psi_mat, q_grid=q_grid, dq=dq)
}

# ------------------------------------------------------------------------------
# EIGENSTATE DENSITY ON A DISPLAY GRID
# ------------------------------------------------------------------------------

#' Probability density |psi_n(q)|^2 of the n-th Schroedinger eigenstate,
#' interpolated to a display grid.
#'
#' Uses the soln object returned by solve_schroedinger() (or by analytic-
#' eigenstate equivalents like harmonic_soln() that mirror the same struct).
#' The wavefunction is already L2-normalized on the solver grid; this
#' function returns its modulus squared on q_display, with linear inter-
#' polation and zero-padding outside the solver grid.
#'
#' @param soln       Output of solve_schroedinger() — list with fields
#'                   energies, psi_matrix, q_grid, dq.
#' @param n          Quantum number (0-indexed).
#' @param q_display  Output grid in q.
#' @return Numeric vector of |psi_n|^2 on q_display.
schroedinger_density <- function(soln, n, q_display) {
  if (n < 0 || n >= length(soln$energies))
    stop(sprintf("n=%d outside available [0, %d]", n, length(soln$energies)-1))
  psi_n  <- soln$psi_matrix[, n + 1]   # column n+1 = state n (0-indexed)
  psi_disp <- approx(soln$q_grid, psi_n, xout=q_display,
                     rule=1, yleft=0, yright=0)$y
  psi_disp[is.na(psi_disp)] <- 0
  psi_disp^2
}

# ------------------------------------------------------------------------------
# EIGENSTATE DENSITY ON A DISPLAY GRID
# ------------------------------------------------------------------------------

#' Probability density |psi_n(q)|^2 of the n-th Schroedinger eigenstate,
#' interpolated to a display grid.
#'
#' Uses the soln object returned by solve_schroedinger() (or by the
#' analytic-eigenstate equivalents like harmonic_soln()). The wavefunction
#' is already L2-normalized on the solver grid; this function returns its
#' modulus squared on q_display, with linear interpolation and zero-padding
#' outside the solver grid.
#'
#' @param soln       Output of solve_schroedinger() — list with fields
#'                   energies, psi_matrix, q_grid, dq.
#' @param n          Quantum number (0-indexed).
#' @param q_display  Output grid in q.
#' @return Numeric vector of |psi_n|^2 on q_display.
schroedinger_density <- function(soln, n, q_display) {
  if (n < 0 || n >= length(soln$energies))
    stop(sprintf("n=%d outside available [0, %d]", n, length(soln$energies)-1))
  psi_n  <- soln$psi_matrix[, n + 1]   # column n+1 = state n (0-indexed)
  psi_disp <- approx(soln$q_grid, psi_n, xout=q_display,
                     rule=1, yleft=0, yright=0)$y
  psi_disp[is.na(psi_disp)] <- 0
  psi_disp^2
}
