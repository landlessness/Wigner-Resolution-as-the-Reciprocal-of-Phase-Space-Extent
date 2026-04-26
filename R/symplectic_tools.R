# ==============================================================================
# symplectic_tools.R
# Shared symplectic geometry tools for resolving Wigner negativity and
# semiclassical caustics via finite classical action.
#
# Units: all positions in units of q_0, momenta in units of p_0.
# This gives hbar = q_0 * p_0 = 1 and A_0 = pi * q_0 * p_0 = h/2.
#
# Reference: de Gosson (2009), Zurek (2001), Robertson (1929), Schrodinger (1930)
# Author: Brian S. Mulloy
# ==============================================================================

library(gsl)
library(data.table)
library(ggforce)

# ------------------------------------------------------------------------------
# UNIT SYSTEM
# Position : q measured in units of q_0
# Momentum : p measured in units of p_0
# Action   : A measured in units of A_0 = pi * q_0 * p_0 = h/2
# hbar = q_0 * p_0 = 1
# ------------------------------------------------------------------------------

robertson_schrodinger <- function(sigma_qq, sigma_pp, sigma_qp = 0, hbar = 1.0) {
  rs_lhs       <- sigma_qq * sigma_pp - sigma_qp^2
  rs_bound     <- (hbar / 2)^2
  rs_satisfied <- rs_lhs >= rs_bound - .Machine$double.eps * abs(rs_bound)
  if (!rs_satisfied) warning(sprintf("RS inequality violated: %.6e < %.6e", rs_lhs, rs_bound))
  Delta_q      <- sqrt(2 * sigma_qq)
  Delta_p      <- sqrt(2 * sigma_pp)
  delta_q      <- hbar / Delta_p
  delta_p      <- hbar / Delta_q
  sp_product   <- delta_q * Delta_p
  sp_satisfied <- abs(sp_product - hbar) < .Machine$double.eps^0.5 * hbar
  if (!sp_satisfied) warning(sprintf("Kernel symplectic positivity not saturated: %.6e (expected %.6e)", sp_product, hbar))
  A_over_A0    <- (Delta_q * Delta_p) / hbar
  list(Delta_q=Delta_q, Delta_p=Delta_p, delta_q=delta_q, delta_p=delta_p, A_over_A0=A_over_A0, rs_satisfied=rs_satisfied, sp_satisfied=sp_satisfied)
}

qho_covariance <- function(n, hbar = 1.0) {
  alpha <- 2*n + 1
  robertson_schrodinger(alpha*hbar/2, alpha*hbar/2, 0, hbar)
}

qho_wigner <- function(n, q, p) {
  rho2 <- q^2 + p^2
  (-1)^n / pi * exp(-rho2) * laguerre_n(n, 0, 2*rho2)
}

squeezed_kernel_q <- function(q, p, rs) {
  (1/pi) * exp(-q^2/rs$delta_q^2 - p^2/rs$Delta_p^2)
}

squeezed_kernel_p <- function(q, p, rs) {
  (1/pi) * exp(-q^2/rs$Delta_q^2 - p^2/rs$delta_p^2)
}

fft_convolve_2d <- function(W_mat, K_mat, dq, dp) {
  nq <- nrow(W_mat); np <- ncol(W_mat)
  K_norm <- K_mat / (sum(K_mat) * dq * dp)
  ifftshift2d <- function(m) {
    nr <- nrow(m); nc <- ncol(m); sr <- floor(nr/2); sc <- floor(nc/2)
    rbind(cbind(m[(sr+1):nr,(sc+1):nc], m[(sr+1):nr,1:sc]), cbind(m[1:sr,(sc+1):nc], m[1:sr,1:sc]))
  }
  K_shift      <- ifftshift2d(K_norm)
  P_mat        <- Re(fft(fft(W_mat) * fft(K_shift), inverse=TRUE)) / (nq*np) * dq*dp
  peak_val     <- max(abs(P_mat))
  tol          <- peak_val * sqrt(.Machine$double.eps) * sqrt(nq*np)
  max_negative <- min(P_mat)
  if (max_negative < -tol) warning(sprintf("fft_convolve_2d: min=%.2e exceeds tol %.2e.", max_negative, -tol))
  P_mat[P_mat < 0 & P_mat >= -tol] <- 0
  list(P_mat=P_mat, max_negative=max_negative, tolerance=tol)
}

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
  if (abs(w_norm-1) > 1e-3) warning(sprintf("n=%d: Wigner norm=%.6f.", n, w_norm))
  w_at_origin  <- wigner_fn(n, 0, 0)
  w_expected   <- (-1)^n / pi
  w_spot_check <- abs(w_at_origin - w_expected)
  if (w_spot_check > 1e-6) warning(sprintf("n=%d: W_n(0,0)=%.6f, expected %.6f", n, w_at_origin, w_expected))
  K_mat    <- outer(q_int, p_int, FUN=function(q,p) kernel_fn(q,p,rs))
  conv     <- fft_convolve_2d(W_mat, K_mat, dq, dp)
  P_mat    <- conv$P_mat
  P_q_int  <- rowSums(P_mat) * dp
  norm_int <- sum(P_q_int) * dq
  if (abs(norm_int) < .Machine$double.eps^0.5) stop(sprintf("n=%d: integrates to zero.", n))
  P_q_int <- P_q_int / norm_int
  P_q     <- approx(q_int, P_q_int, xout=q_display, rule=1)$y
  P_q[is.na(P_q)] <- 0
  list(P_q=P_q, w_norm=w_norm, w_spot_check=w_spot_check, max_negative=conv$max_negative, tolerance=conv$tolerance)
}

# ------------------------------------------------------------------------------
# NUMERICAL WIGNER FUNCTION (for anharmonic potentials)
# Computes Wigner function numerically from wavefunction psi(q) on a grid.
# Used for Morse and double-well where no closed-form Wigner exists.
# ------------------------------------------------------------------------------

#' Compute Wigner function numerically from wavefunction on a grid.
#' W(q,p) = (1/pi) * integral psi*(q+x) psi(q-x) exp(2ipx) dx
#' @param psi_fn Function(q) returning complex wavefunction values
#' @param q_grid Position grid for Wigner evaluation
#' @param p_grid Momentum grid for Wigner evaluation
#' @param x_half Half-width of integration in x (should exceed wavefunction support)
#' @param nx Number of integration points in x
#' @return Matrix W[iq, ip] of Wigner function values
wigner_numerical <- function(psi_fn, q_grid, p_grid, x_half=10, nx=501) {
  x_int <- seq(-x_half, x_half, length.out=nx)
  dx    <- diff(x_int)[1]
  nq    <- length(q_grid); np <- length(p_grid)
  W_mat <- matrix(0, nrow=nq, ncol=np)
  for (iq in seq_len(nq)) {
    q0      <- q_grid[iq]
    psi_plus  <- psi_fn(q0 + x_int)
    psi_minus <- Conj(psi_fn(q0 - x_int))
    integrand <- psi_minus * psi_plus
    for (ip in seq_len(np)) {
      p0         <- p_grid[ip]
      kernel     <- exp(2i * p0 * x_int)
      W_mat[iq,ip] <- Re(sum(integrand * kernel) * dx) / pi
    }
  }
  W_mat
}

#' Compute Robertson-Schrodinger covariance from numerical wavefunction.
#' @param psi_fn Function(q) returning wavefunction values
#' @param q_grid Position grid
#' @param hbar Reduced Planck constant (default 1.0)
#' @return Output of robertson_schrodinger()
numerical_covariance <- function(psi_fn, q_grid, hbar=1.0) {
  dq      <- diff(q_grid)[1]
  psi_vec <- psi_fn(q_grid)
  prob    <- Re(Conj(psi_vec) * psi_vec)
  prob    <- prob / (sum(prob) * dq)
  q_mean  <- sum(q_grid * prob) * dq
  sigma_qq <- sum((q_grid - q_mean)^2 * prob) * dq
  # Momentum space via FFT
  n_pts   <- length(q_grid)
  psi_k   <- fft(psi_vec) * dq
  dk      <- 2*pi / (n_pts * dq)
  k_grid  <- c(seq(0, (n_pts/2-1)*dk, by=dk), seq(-(n_pts/2)*dk, -dk, by=dk))
  prob_k  <- Re(Conj(psi_k) * psi_k)
  prob_k  <- prob_k / (sum(prob_k) * dk)
  p_mean  <- sum(k_grid * prob_k) * dk
  sigma_pp <- sum((k_grid - p_mean)^2 * prob_k) * dk
  # qp covariance: <(q-<q>)(p-<p>) + (p-<p>)(q-<q>)>/2
  # For real wavefunctions sigma_qp = 0 by symmetry
  sigma_qp <- 0
  robertson_schrodinger(sigma_qq, sigma_pp, sigma_qp, hbar)
}

# ------------------------------------------------------------------------------
# DISPLAY GEOMETRY
# Computes display limits and axis breaks — identical across all figure files.
# ------------------------------------------------------------------------------

#' Compute display geometry from classical boundary.
#' @param Delta_q Classical position semi-axis
#' @param ell_scale Fractional padding beyond Delta_q for heatmap (default 1.25)
#' @param plot_extra Absolute padding for cross-section display (default 2.5)
display_geometry <- function(Delta_q, ell_scale=1.25, plot_extra=2.5) {
  ell_lim       <- Delta_q * ell_scale
  plot_lim      <- max(Delta_q * ell_scale, Delta_q + plot_extra)
  break_val     <- min(round(Delta_q,1), floor(ell_lim*10)/10)
  custom_breaks <- c(-break_val, 0, break_val)
  label_format  <- function(x) sprintf("%.1f", x)
  list(ell_lim=ell_lim, plot_lim=plot_lim, custom_breaks=custom_breaks, label_format=label_format)
}

# ------------------------------------------------------------------------------
# HEATMAP DATA
# Builds the 2D Wigner heatmap data table with contrast boost.
# Identical across all figure files.
# ------------------------------------------------------------------------------

#' Build 2D Wigner heatmap data with contrast boost.
#' @param n_val Quantum number (passed to wigner_fn)
#' @param wigner_fn Function(n, q, p) returning Wigner values
#' @param ell_lim Display half-width
#' @param grid_pts Resolution of heatmap grid (default 400)
wigner_heatmap_data <- function(n_val, wigner_fn, ell_lim, grid_pts=400) {
  q_ell <- seq(-ell_lim, ell_lim, length.out=grid_pts)
  p_ell <- seq(-ell_lim, ell_lim, length.out=grid_pts)
  dt    <- as.data.table(expand.grid(q=q_ell, p=p_ell))
  dt[, w := wigner_fn(n_val, q, p)]
  dt[, w_plot := sign(w) * abs(w)^0.4]
  max_w <- max(abs(dt$w_plot), na.rm=TRUE)
  if (max_w > 0) dt[, w_plot := w_plot/max_w]
  dt
}

#' Build 2D Wigner heatmap from precomputed matrix (for numerical wavefunctions).
#' @param W_mat Matrix of Wigner values (nq x np)
#' @param q_grid Position grid
#' @param p_grid Momentum grid
wigner_heatmap_data_from_matrix <- function(W_mat, q_grid, p_grid) {
  dt <- as.data.table(expand.grid(q=q_grid, p=p_grid))
  dt[, w := as.vector(W_mat)]
  dt[, w_plot := sign(w) * abs(w)^0.4]
  max_w <- max(abs(dt$w_plot), na.rm=TRUE)
  if (max_w > 0) dt[, w_plot := w_plot/max_w]
  dt
}

# ------------------------------------------------------------------------------
# ROW LABEL
# ------------------------------------------------------------------------------

#' Build row label panel.
plot_row_label <- function(label_str, parse=TRUE, base_font="") {
  ggplot() + theme_void() +
    coord_cartesian(xlim=c(-1,1), ylim=c(0,1), clip="off") +
    annotate("text", x=0, y=0.5, label=label_str, parse=parse, family=base_font, size=4.5, hjust=0.5)
}

# ------------------------------------------------------------------------------
# PHASE SPACE HEATMAP PLOT
# ------------------------------------------------------------------------------

#' Build phase space heatmap ggplot panel.
plot_phase_space_heatmap <- function(dt_w2d, ell_data, ell_lim, custom_breaks, label_format, base_font="") {
  ax_x <- expression(italic(q)/italic(q)[0])
  ax_y <- expression(italic(p)/italic(p)[0])
  ggplot(dt_w2d, aes(x=q, y=p)) +
    symplectic_ellipse_layers_bottom(ell_data) +
    geom_raster(aes(fill=w_plot), interpolate=TRUE) +
    scale_fill_gradient2(low="gray10", mid="white", high="gray40", midpoint=0, limits=c(-1,1), guide="none") +
    symplectic_ellipse_layers_top(ell_data) +
    coord_fixed(xlim=c(-ell_lim,ell_lim), ylim=c(-ell_lim,ell_lim), expand=FALSE) +
    scale_x_continuous(breaks=custom_breaks, labels=label_format) +
    scale_y_continuous(breaks=custom_breaks, labels=label_format) +
    theme_bw(base_family=base_font) +
    theme(panel.grid.minor=element_blank(), panel.background=element_rect(fill="white"), axis.text=element_text(size=8), plot.margin=margin(2,4,2,4)) +
    labs(x=ax_x, y=ax_y)
}

# ------------------------------------------------------------------------------
# SAVE FIGURE
# ------------------------------------------------------------------------------

#' Save figure to PDF with standard dimensions.
save_figure <- function(p, filepath, n_rows, fig_width=7.0) {
  ggsave(filename=filepath, plot=p, device=cairo_pdf, width=fig_width, height=n_rows*1.8+0.5, limitsize=FALSE)
  cat("Done.", filepath, "\n")
}

# ------------------------------------------------------------------------------
# SYMPLECTIC ELLIPSE OVERLAYS
# ------------------------------------------------------------------------------

symplectic_ellipse_data <- function(rs, r_system=NULL) {
  r_fermi <- rs$Delta_q
  if (is.null(r_system)) r_system <- r_fermi
  df_system <- data.frame(x0=0, y0=0, r=r_system)
  df_fermi  <- data.frame(x0=0, y0=0, r=r_fermi)
  df_cigars <- data.frame(x0=0, y0=0, aq_a=rs$delta_q, aq_b=rs$Delta_p, ap_a=rs$Delta_q, ap_b=rs$delta_p)
  list(system=df_system, fermi=df_fermi, cigars=df_cigars)
}

symplectic_ellipse_layers_bottom <- function(ell_data) {
  list(geom_circle(data=ell_data$system, aes(x0=x0, y0=y0, r=r), inherit.aes=FALSE, fill="white", color=NA))
}

symplectic_ellipse_layers_top <- function(ell_data) {
  list(
    geom_circle(data=ell_data$system, aes(x0=x0, y0=y0, r=r), inherit.aes=FALSE, color="black", linewidth=0.5, linetype="solid"),
    geom_circle(data=ell_data$fermi,  aes(x0=x0, y0=y0, r=r), inherit.aes=FALSE, color="gray20", linewidth=0.4, linetype="dashed"),
    geom_ellipse(data=ell_data$cigars, aes(x0=x0, y0=y0, a=aq_a, b=aq_b, angle=0), inherit.aes=FALSE, color="gray20", linewidth=0.4, linetype="dotted"),
    geom_ellipse(data=ell_data$cigars, aes(x0=x0, y0=y0, a=ap_a, b=ap_b, angle=0), inherit.aes=FALSE, color="gray20", linewidth=0.4, linetype="dotted")
  )
}

symplectic_ellipse_layers <- function(ell_data) {
  c(symplectic_ellipse_layers_bottom(ell_data), symplectic_ellipse_layers_top(ell_data))
}
