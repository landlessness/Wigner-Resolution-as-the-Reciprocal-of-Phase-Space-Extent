# ==============================================================================
# semiclassical_density.R
# Semiclassical phase-space densities and the per-state bundle that wraps them.
#
# The semiclassical regime takes V(q) and energy E_n as inputs and produces
# phase-space content WITHOUT solving Schroedinger:
#
#   wkb_caustic_density(q, E, V_fn): 1D analytical density
#     P_mc(q) = 1 / (T(E) * |p(q)|), where p(q) = sqrt(2*(E - V(q)))
#     and T(E) is the classical orbit period:
#       T(E) = oint dq / |p(q)|
#     Diverges at classical turning points where E = V(q) — the caustic.
#
#   semiclassical_shell_density(q, p, E, V_fn, epsilon): regularized 2D shell
#     W_cl(q,p) = (1/Z) * exp(-(H(q,p) - E)^2 / (2*epsilon^2))
#     where H(q,p) = p^2/2 + V(q). The regularization is purely numerical
#     (the discrete grid cannot represent a Dirac delta); epsilon is set
#     small relative to E_n and the kernel widths so that downstream
#     convolution is insensitive to it.
#
#   build_semiclassical_state(...): orchestrates the above into a per-state
#     bundle compatible in shape with build_wigner_state(), so kernel
#     application code can treat them uniformly.
#
# Kernel-sizing widths (Delta_q, Delta_p) for the symplectic kernel come
# from orbit_covariance() in classical_action_tools.R — also Schroedinger-free.
# The whole semiclassical pipeline is independent of the wavefunction.
#
# Reference: Berry 1977 J.Phys.A 10, 2083 (semiclassical phase space)
#            Berry & Mount 1972 Rep.Prog.Phys. 35, 315 (WKB caustics)
# Author: Brian S. Mulloy
# ==============================================================================

library(here)
library(data.table)
source(here("R", "math_tools.R"))

# Default regularization for the 2D shell:
# epsilon = SEMICLASSICAL_REG_FRACTION * |E_n|. Small relative to the energy
# scale, but large enough to span several grid cells of the integration grid.
SEMICLASSICAL_REG_FRACTION <- 0.02

# ------------------------------------------------------------------------------
# 1D WKB CAUSTIC DENSITY (TIME-AVERAGED, smooth)
# ------------------------------------------------------------------------------

#' Analytical WKB density at energy E for potential V_fn on q_grid.
#'
#' Computes P_mc(q) = 1 / (T(E) * |p(q)|) where p(q) = sqrt(2*(E - V(q))).
#' This is the time-averaged classical density — fraction of period spent
#' near q. Diverges at classical turning points where E = V(q); these
#' points are rendered as Inf in the returned vector. Caller is expected
#' to handle the infinities (e.g., clip in plot, render with arrows).
#'
#' Note: the OSCILLATING WKB |psi_WKB|^2 — which is what the canonical
#' middle column of the semiclassical figures actually shows — is in
#' oscillating_wkb_density() below. This function is retained for
#' analytic comparison and for any callers that want the smooth envelope.
#'
#' @param q_grid Position grid (numeric vector).
#' @param E Energy.
#' @param V_fn Function(q) returning potential values.
#' @return Numeric vector of P_mc values, with Inf at turning points and 0
#'         in classically forbidden region E < V(q).
wkb_caustic_density <- function(q_grid, E, V_fn) {
  V_vals <- V_fn(q_grid)
  arg    <- 2 * (E - V_vals)
  p_q    <- sqrt(pmax(arg, 0))   # |p(q)| where defined; 0 in forbidden region

  # Compute period via fine-grid integration over the classically-allowed region.
  q_fine  <- seq(min(q_grid), max(q_grid), length.out=10001)
  V_fine  <- V_fn(q_fine)
  arg_f   <- 2 * (E - V_fine)
  allowed <- arg_f > 0
  # Tiny floor avoids 0/0; integrand falls off as 1/sqrt away from the
  # turning point, integrable on a fine grid with trapezoidal quadrature.
  p_fine    <- sqrt(pmax(arg_f, 1e-15))
  dq_fine   <- diff(q_fine)[1]
  integrand <- ifelse(allowed, 1/p_fine, 0)
  T_E       <- sum(integrand) * dq_fine

  # Caustic density: 1 / (T * |p|), with Inf at turning points (p=0).
  caustic <- ifelse(p_q > 0, 1/(T_E * p_q), Inf)
  caustic[!is.finite(caustic) & p_q == 0 & arg == 0] <- Inf  # turning points
  caustic[arg < 0] <- 0  # forbidden region: zero, not infinity

  caustic
}

# ------------------------------------------------------------------------------
# 1D OSCILLATING WKB DENSITY |psi_WKB(q)|^2
#
# psi_WKB(q) = (A / sqrt(p(q))) * cos(S(q)/hbar - pi/4)
# |psi_WKB(q)|^2 = (A^2 / p(q)) * cos^2(S(q)/hbar - pi/4)
#
# where S(q) = int_{q_-}^q p(q') dq' is the WKB action measured from
# the inner turning point q_-. Oscillates n+1 times between turning
# points (Bohr-Sommerfeld phase ramps up by pi(n+1/2)). Diverges at
# turning points like the time-averaged caustic — the cos^2 modulation
# preserves the 1/p envelope's singularity.
#
# Normalization: A^2 chosen so int |psi_WKB|^2 dq = 1 over the orbit.
# Since <cos^2> = 1/2 and int 1/p dq over half-orbit = T/2, we have
# A^2 = 2/T(E).
#
# This is the canonical Griffiths-Ch.8 WKB density displayed in the
# middle column of the semiclassical figures.
# ------------------------------------------------------------------------------

#' Oscillating WKB density |psi_WKB(q)|^2 on q_grid.
#'
#' @param q_grid    Position grid (numeric vector).
#' @param E         Energy.
#' @param V_fn      Function(q) returning potential values.
#' @param q_minus   Inner turning point (action reference).
#' @param q_plus    Outer turning point.
#' @param hbar      Planck constant in chosen units.
#' @param n_action  Number of integration points for the action integral
#'                  per cell (default 2001).
#' @return Numeric vector of |psi_WKB|^2 values on q_grid; 0 in the
#'         forbidden region, Inf at turning points.
oscillating_wkb_density <- function(q_grid, E, V_fn, q_minus, q_plus,
                                    hbar=1.0, n_action=2001) {
  V_vals  <- V_fn(q_grid)
  arg     <- 2 * (E - V_vals)
  allowed <- arg > 0
  p_q     <- sqrt(pmax(arg, 0))

  # Period for normalization
  q_fine  <- seq(q_minus, q_plus, length.out=10001)
  arg_f   <- 2 * (E - V_fn(q_fine))
  p_fine  <- sqrt(pmax(arg_f, 1e-15))
  dq_fine <- diff(q_fine)[1]
  T_E     <- 2 * sum(1/p_fine) * dq_fine
  A_sq    <- 2 / T_E

  # Action integral S(q) from q_minus to each q in q_grid (clipped to orbit)
  S <- numeric(length(q_grid))
  for (k in seq_along(q_grid)) {
    q <- q_grid[k]
    if (q <= q_minus) { S[k] <- 0; next }
    q_eff <- min(q, q_plus)
    qq <- seq(q_minus, q_eff, length.out=n_action)
    S[k] <- sum(sqrt(2*pmax(E - V_fn(qq), 0))) * diff(qq)[1]
  }

  cos2 <- cos(S/hbar - pi/4)^2

  # |psi_WKB|^2 = (A^2/p) * cos^2; 0 in forbidden, Inf at turning points
  rho <- ifelse(allowed & p_q > 0, A_sq * cos2 / p_q, 0)
  # Mark turning points as Inf so the plot helper can draw infinity arrows
  rho[allowed & p_q == 0] <- Inf
  rho
}

# ------------------------------------------------------------------------------
# 2D REGULARIZED ENERGY SHELL (smooth, time-averaged)
# ------------------------------------------------------------------------------

#' Regularized classical energy shell on a (q,p) grid.
#'
#' W_cl(q,p) = (1/Z) * exp(-(H(q,p) - E)^2 / (2*epsilon^2))
#' where H(q,p) = p^2/2 + V(q). Z normalizes so that integral W_cl dq dp = 1.
#'
#' Smooth shell — its q-marginal is the time-averaged WKB caustic
#' 1/(T*|p|). Used historically; now superseded by wkb_phase_space_lift()
#' below for the oscillating-input pipeline. Retained for backward
#' compatibility and for any callers that want the smooth shell.
#'
#' @param q_grid Position grid.
#' @param p_grid Momentum grid.
#' @param E Energy.
#' @param V_fn Function(q) returning potential values.
#' @param epsilon Regularization width (energy units).
#' @return Matrix W_cl[iq, ip].
semiclassical_shell_density <- function(q_grid, p_grid, E, V_fn, epsilon) {
  V_vals <- V_fn(q_grid)
  H_mat  <- outer(V_vals, p_grid^2/2, FUN="+")
  W_mat  <- exp(-(H_mat - E)^2 / (2*epsilon^2))
  dq     <- diff(q_grid)[1]
  dp     <- diff(p_grid)[1]
  Z      <- sum(W_mat) * dq * dp
  if (Z > 0) W_mat <- W_mat / Z
  W_mat
}

# ------------------------------------------------------------------------------
# 2D WKB PHASE-SPACE LIFT (oscillating)
#
# Lifts |psi_WKB(q)|^2 onto phase space:
#   W_WKB(q,p) = (1/Z) * exp(-(H(q,p) - E)^2 / (2*eps^2)) * cos^2(S(q)/hbar - pi/4)
#
# The Gaussian shell concentrates probability on the energy contour
# H(q,p) = E. The cos^2 factor modulates along the contour with the
# WKB phase. NO 1/p factor — that emerges from the projection onto q
# (integrating across the contour band, which is steeper-in-p where
# |p| is small, gathering more length per dq there). Including a 1/p
# in W would double-count the divergence.
#
# Marginalizing W_WKB over p reproduces the oscillating WKB density
# |psi_WKB|^2 = (A^2/p)*cos^2(S/hbar - pi/4), up to the regularization
# softening at turning points.
# ------------------------------------------------------------------------------

#' Phase-space lift of oscillating |psi_WKB|^2 onto a (q,p) grid.
#'
#' @param q_grid    Position grid.
#' @param p_grid    Momentum grid.
#' @param E         Energy.
#' @param V_fn      Function(q) returning potential values.
#' @param q_minus   Inner turning point (action reference).
#' @param q_plus    Outer turning point.
#' @param epsilon   Regularization width (energy units).
#' @param hbar      Planck constant in chosen units.
#' @return Matrix W_WKB[iq, ip], normalized so integral = 1.
wkb_phase_space_lift <- function(q_grid, p_grid, E, V_fn,
                                 q_minus, q_plus, epsilon, hbar=1.0,
                                 n_action=2001) {
  V_vals <- V_fn(q_grid)
  H_mat  <- outer(V_vals, p_grid^2/2, FUN="+")
  W_shell <- exp(-(H_mat - E)^2 / (2*epsilon^2))

  # Action S(q) at each q_grid point, from q_minus to q (clipped to orbit)
  S_q <- numeric(length(q_grid))
  for (k in seq_along(q_grid)) {
    q <- q_grid[k]
    if (q <= q_minus) { S_q[k] <- 0; next }
    q_eff <- min(q, q_plus)
    qq <- seq(q_minus, q_eff, length.out=n_action)
    S_q[k] <- sum(sqrt(2*pmax(E - V_fn(qq), 0))) * diff(qq)[1]
  }
  cos2_q <- cos(S_q/hbar - pi/4)^2

  # Modulate shell by cos^2(S/hbar - pi/4) along q
  W_mat <- W_shell * matrix(cos2_q, nrow=length(q_grid), ncol=length(p_grid))

  dq <- diff(q_grid)[1]
  dp <- diff(p_grid)[1]
  Z  <- sum(W_mat) * dq * dp
  if (Z > 0) W_mat <- W_mat / Z
  W_mat
}

# ------------------------------------------------------------------------------
# BUILD SEMICLASSICAL STATE
# ------------------------------------------------------------------------------

#' Build a per-state semiclassical bundle.
#'
#' Drop-in replacement for build_wigner_state(): returns a bundle with the
#' same field names (q_int, p_int, dq_int, dp_int, W_matrix, W_cross,
#' heatmap_dt, norm) so apply_kernel_cross_section() works without
#' modification.
#'
#' Now uses the OSCILLATING WKB phase-space lift as the W_matrix (the
#' input to the symplectic kernel convolution), and the OSCILLATING WKB
#' density |psi_WKB(q)|^2 as wkb_density (the analytic 1D content for
#' the middle column of the figures). This is the canonical Griffiths-
#' Ch.8 starting point: |psi_WKB|^2 oscillates with n+1 lobes between
#' turning points and diverges at the turning points themselves.
#'
#' Also stores wkb_density_smooth (the time-averaged 1/(T*|p|)) for
#' reference / backward compatibility.
#'
#' @param E_n Energy of the eigenstate.
#' @param V_fn Function(q) returning potential values.
#' @param q_lo,q_hi Display window in q.
#' @param p_lo,p_hi Display window in p.
#' @param q_display Position grid for cross-section / projection rendering.
#' @param q_minus,q_plus Classical turning points at energy E_n. Required
#'                       for the oscillating WKB construction; the action
#'                       integral is measured from q_minus and the WKB
#'                       density is supported on [q_minus, q_plus].
#' @param epsilon Regularization width for 2D shell. If NULL, defaults to
#'                SEMICLASSICAL_REG_FRACTION * |E_n|.
#' @param nq_int,np_int Integration grid resolution (default 801 x 601).
#' @return List with fields matching build_wigner_state() bundle plus
#'         wkb_density (oscillating |psi_WKB|^2) and wkb_density_smooth
#'         (time-averaged 1/(T|p|)).
build_semiclassical_state <- function(E_n, V_fn,
                                      q_lo, q_hi, p_lo, p_hi, q_display,
                                      q_minus, q_plus,
                                      epsilon=NULL,
                                      nq_int=801, np_int=601) {
  if (is.null(epsilon)) epsilon <- SEMICLASSICAL_REG_FRACTION * abs(E_n)

  q_int  <- seq(q_lo, q_hi, length.out=nq_int)
  p_int  <- seq(p_lo, p_hi, length.out=np_int)
  dq_int <- diff(q_int)[1]
  dp_int <- diff(p_int)[1]

  cat(sprintf("    Building oscillating WKB phase-space lift on %d x %d grid (epsilon=%.4f)...\n",
              nq_int, np_int, epsilon))

  W_mat <- wkb_phase_space_lift(q_int, p_int, E_n, V_fn,
                                q_minus, q_plus, epsilon)
  norm  <- sum(W_mat) * dq_int * dp_int
  cat(sprintf("    Shell norm: %.6f\n", norm))

  # Extract p=0 cross-section (interpolated to q_display)
  W_cross <- extract_p0_cross_section(W_mat, q_int, p_int, q_display)

  # 1D oscillating WKB density on q_display (canonical for middle column)
  wkb_osc    <- oscillating_wkb_density(q_display, E_n, V_fn,
                                        q_minus, q_plus)
  # Time-averaged smooth caustic on q_display (kept for reference)
  wkb_smooth <- wkb_caustic_density(q_display, E_n, V_fn)

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
    q_int              = q_int,
    p_int              = p_int,
    dq_int             = dq_int,
    dp_int             = dp_int,
    W_matrix           = W_mat,
    W_cross            = W_cross,
    wkb_density        = wkb_osc,      # oscillating, canonical
    wkb_density_smooth = wkb_smooth,   # time-averaged, kept for reference
    heatmap_dt         = heatmap_dt,
    norm               = norm,
    epsilon            = epsilon
  )
}
