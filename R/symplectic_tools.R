# ==============================================================================
# symplectic_tools.R
# Shared symplectic geometry tools for resolving Wigner negativity and
# semiclassical caustics via finite classical action.
#
# Units: all positions in units of q_0, momenta in units of p_0.
# This gives hbar = q_0 * p_0 = 1 and A_0 = pi * q_0 * p_0 = h/2.
# Action ratios A/A_0 = 2n+1 for the QHO, derived from Robertson-Schrodinger.
#
# Reference: de Gosson (2009), Zurek (2001), Robertson (1929), Schrodinger (1930)
# Author: Brian S. Mulloy
# ==============================================================================

library(gsl)
library(data.table)
library(ggforce)

# ------------------------------------------------------------------------------
# UNIT SYSTEM
# Position : q measured in units of q_0  (ground state position uncertainty)
# Momentum : p measured in units of p_0  (ground state momentum uncertainty)
# Action   : A measured in units of A_0 = pi * q_0 * p_0 = h/2
# In these units: hbar = q_0 * p_0 = 1
# ------------------------------------------------------------------------------

#' Robertson-Schrodinger symplectic geometry from a covariance matrix.
#' @param sigma_qq Position variance
#' @param sigma_pp Momentum variance
#' @param sigma_qp Position-momentum covariance (0 for QHO eigenstates)
#' @param hbar Reduced Planck constant (default 1.0 in q_0*p_0 units)
#' @return Named list: Delta_q, Delta_p, delta_q, delta_p, A_over_A0,
#'         rs_satisfied, sp_satisfied
robertson_schrodinger <- function(sigma_qq, sigma_pp, sigma_qp = 0, hbar = 1.0) {
  rs_lhs       <- sigma_qq * sigma_pp - sigma_qp^2
  rs_bound     <- (hbar / 2)^2
  rs_satisfied <- rs_lhs >= rs_bound - .Machine$double.eps * abs(rs_bound)
  if (!rs_satisfied) warning(sprintf("RS inequality violated: %.6e < %.6e", rs_lhs, rs_bound))

  Delta_q <- sqrt(2 * sigma_qq)
  Delta_p <- sqrt(2 * sigma_pp)
  delta_q <- hbar / Delta_p
  delta_p <- hbar / Delta_q

  sp_product   <- delta_q * Delta_p
  sp_satisfied <- abs(sp_product - hbar) < .Machine$double.eps^0.5 * hbar
  if (!sp_satisfied) warning(sprintf("Kernel symplectic positivity not saturated: delta_q * Delta_p = %.6e (expected %.6e)", sp_product, hbar))

  # A/A_0 = (pi * Delta_q * Delta_p) / (pi * hbar) = Delta_q * Delta_p / hbar
  A_over_A0 <- (Delta_q * Delta_p) / hbar

  list(Delta_q=Delta_q, Delta_p=Delta_p, delta_q=delta_q, delta_p=delta_p, A_over_A0=A_over_A0, rs_satisfied=rs_satisfied, sp_satisfied=sp_satisfied)
}

#' QHO Robertson-Schrodinger geometry for eigenstate n.
#' sigma_qq = sigma_pp = (2n+1)*hbar/2, sigma_qp = 0, A/A_0 = 2n+1.
qho_covariance <- function(n, hbar = 1.0) {
  alpha <- 2*n + 1
  robertson_schrodinger(alpha*hbar/2, alpha*hbar/2, 0, hbar)
}

#' True QHO Wigner function W_n(q,p) in q_0,p_0 units.
#' W_n(q,p) = (-1)^n/pi * exp(-(q^2+p^2)) * L_n(2*(q^2+p^2))
#' Spot-check: W_n(0,0) = (-1)^n / pi
qho_wigner <- function(n, q, p) {
  rho2 <- q^2 + p^2
  (-1)^n / pi * exp(-rho2) * laguerre_n(n, 0, 2*rho2)
}

#' Squeezed kernel G_delta_q(q,p): squeezed in q (width delta_q), spans Delta_p.
#' Symplectic area = h/2. Prefactor 1/pi in q_0,p_0 units (hbar=1).
squeezed_kernel_q <- function(q, p, rs) {
  (1/pi) * exp(-q^2/rs$delta_q^2 - p^2/rs$Delta_p^2)
}

#' Conjugate kernel G_delta_p(q,p): squeezed in p (width delta_p), spans Delta_q.
squeezed_kernel_p <- function(q, p, rs) {
  (1/pi) * exp(-q^2/rs$Delta_q^2 - p^2/rs$delta_p^2)
}

#' 2D FFT convolution P = W * K.
fft_convolve_2d <- function(W_mat, K_mat, dq, dp) {
  nq <- nrow(W_mat); np <- ncol(W_mat)
  K_norm <- K_mat / (sum(K_mat) * dq * dp)
  ifftshift2d <- function(m) {
    nr <- nrow(m); nc <- ncol(m); sr <- floor(nr/2); sc <- floor(nc/2)
    rbind(cbind(m[(sr+1):nr,(sc+1):nc], m[(sr+1):nr,1:sc]), cbind(m[1:sr,(sc+1):nc], m[1:sr,1:sc]))
  }
  K_shift      <- ifftshift2d(K_norm)
  P_mat        <- Re(fft(fft(W_mat) * fft(K_shift), inverse=TRUE)) / (nq*np)
  P_mat        <- P_mat * dq * dp
  peak_val     <- max(abs(P_mat))
  tol          <- peak_val * sqrt(.Machine$double.eps) * sqrt(nq*np)
  max_negative <- min(P_mat)
  if (max_negative < -tol) warning(sprintf("fft_convolve_2d: min(P_mat)=%.2e exceeds tolerance %.2e.", max_negative, -tol))
  P_mat[P_mat < 0 & P_mat >= -tol] <- 0
  list(P_mat=P_mat, max_negative=max_negative, tolerance=tol)
}

#' Compute symplectic density on display grid via 2D convolution.
compute_symplectic_density <- function(n, wigner_fn, kernel_fn, rs, q_display, hbar=1.0) {
  Delta_q   <- rs$Delta_q
  integ_lim <- max(Delta_q * 2.0, 4.0 * sqrt(hbar))
  integ_res <- max(601, 20 * ceiling(rs$A_over_A0))
  if (integ_res %% 2 == 0) integ_res <- integ_res + 1
  q_int <- seq(-integ_lim, integ_lim, length.out=integ_res)
  p_int <- seq(-integ_lim, integ_lim, length.out=integ_res)
  dq    <- diff(q_int)[1]; dp <- diff(p_int)[1]
  W_mat  <- outer(q_int, p_int, FUN=function(q,p) wigner_fn(n,q,p))
  w_norm <- sum(W_mat) * dq * dp
  if (abs(w_norm-1) > 1e-3) warning(sprintf("n=%d: Wigner norm=%.6f. Increase integ_lim or integ_res.", n, w_norm))
  w_at_origin  <- wigner_fn(n, 0, 0)
  w_expected   <- (-1)^n / pi
  w_spot_check <- abs(w_at_origin - w_expected)
  if (w_spot_check > 1e-6) warning(sprintf("n=%d: W_n(0,0)=%.6f, expected %.6f (error %.2e)", n, w_at_origin, w_expected, w_spot_check))
  K_mat    <- outer(q_int, p_int, FUN=function(q,p) kernel_fn(q,p,rs))
  conv     <- fft_convolve_2d(W_mat, K_mat, dq, dp)
  P_mat    <- conv$P_mat
  P_q_int  <- rowSums(P_mat) * dp
  norm_int <- sum(P_q_int) * dq
  if (abs(norm_int) < .Machine$double.eps^0.5) stop(sprintf("n=%d: P_delta_q integrates to zero.", n))
  P_q_int <- P_q_int / norm_int
  P_q     <- approx(q_int, P_q_int, xout=q_display, rule=1)$y
  P_q[is.na(P_q)] <- 0
  list(P_q=P_q, w_norm=w_norm, w_spot_check=w_spot_check, max_negative=conv$max_negative, tolerance=conv$tolerance)
}

# ------------------------------------------------------------------------------
# SYMPLECTIC ELLIPSE OVERLAYS
# Three functions:
#   symplectic_ellipse_data()         — build data frames from RS geometry
#   symplectic_ellipse_layers_bottom() — white fill only, drawn BEFORE raster
#   symplectic_ellipse_layers_top()    — outlines only, drawn AFTER raster
#   symplectic_ellipse_layers()        — both combined, for non-raster plots
#
# For the harmonic oscillator r_system = Delta_q so system and Fermi boundaries
# coincide — a visual statement that A = A_ho. For anharmonic systems pass the
# classical turning point radius as r_system and they separate automatically.
# ------------------------------------------------------------------------------

#' Build data frames for symplectic ellipse overlays.
#' @param rs Output of robertson_schrodinger() or qho_covariance()
#' @param r_system Classical system boundary radius. If NULL, defaults to Delta_q.
#' @return Named list: system, fermi, cigars
symplectic_ellipse_data <- function(rs, r_system = NULL) {
  r_fermi <- rs$Delta_q
  if (is.null(r_system)) r_system <- r_fermi
  df_system <- data.frame(x0=0, y0=0, r=r_system)
  df_fermi  <- data.frame(x0=0, y0=0, r=r_fermi)
  df_cigars <- data.frame(x0=0, y0=0, aq_a=rs$delta_q, aq_b=rs$Delta_p, ap_a=rs$Delta_q, ap_b=rs$delta_p)
  list(system=df_system, fermi=df_fermi, cigars=df_cigars)
}

#' Bottom ellipse layers: white fill only. Draw BEFORE geom_raster.
symplectic_ellipse_layers_bottom <- function(ell_data) {
  list(
    geom_circle(data=ell_data$system, aes(x0=x0, y0=y0, r=r), inherit.aes=FALSE, fill="white", color=NA)
  )
}

#' Top ellipse layers: outlines only. Draw AFTER geom_raster.
#' System boundary : solid black  — physical energy shell
#' Fermi blob      : dashed gray  — RS covariance ellipse
#' Quantum blobs   : dotted gray  — Zurek reciprocal scale ellipses
symplectic_ellipse_layers_top <- function(ell_data) {
  list(
    geom_circle(data=ell_data$system, aes(x0=x0, y0=y0, r=r), inherit.aes=FALSE, color="black", linewidth=0.5, linetype="solid"),
    geom_circle(data=ell_data$fermi,  aes(x0=x0, y0=y0, r=r), inherit.aes=FALSE, color="gray20", linewidth=0.4, linetype="dashed"),
    geom_ellipse(data=ell_data$cigars, aes(x0=x0, y0=y0, a=aq_a, b=aq_b, angle=0), inherit.aes=FALSE, color="gray20", linewidth=0.4, linetype="dotted"),
    geom_ellipse(data=ell_data$cigars, aes(x0=x0, y0=y0, a=ap_a, b=ap_b, angle=0), inherit.aes=FALSE, color="gray20", linewidth=0.4, linetype="dotted")
  )
}

#' Combined layers for non-raster plots.
symplectic_ellipse_layers <- function(ell_data) {
  c(symplectic_ellipse_layers_bottom(ell_data), symplectic_ellipse_layers_top(ell_data))
}
