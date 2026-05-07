# ==============================================================================
# tomography_compute.R
# Squeezed-uncertainty tomography pipeline. Builds the state bundle for
# each of the eight states using the same shared state_builder.R logic
# the existing plot_wigner.R and plot_cats.R rely on, then adds:
#   - angular sweep theta in [0, pi)
#   - rotated kernel G_{delta x_theta} convolution at each theta
#   - 1D projection at each theta -> sinogram
#   - filtered back-projection of the sinogram -> tilde_W_delta
#   - 1D marginals of tilde_W_delta over p (rho_q) and over q (rho_p)
#   - exact |psi(q)|^2 and |psi-hat(p)|^2 for comparison
#
# Output: data/tomography_data.rds, consumed by the five
# plot_tomography_*.R figure renderers.
#
# Author: Brian S. Mulloy
# ==============================================================================

library(here)
library(data.table)

source(here("R", "math_tools.R"))
source(here("R", "quantum_tools.R"))
source(here("R", "wigner_density.R"))
source(here("R", "symplectic_kernel.R"))
source(here("R", "marginal_tools.R"))
source(here("R", "schroedinger_solver.R"))

source(here("R", "harmonic_system.R"))
source(here("R", "morse_system.R"))
source(here("R", "double_well_system.R"))
source(here("R", "cat_system.R"))

source(here("R", "state_builder.R"))

OUTPUT_RDS <- here("data", "tomography_data.rds")

# ------------------------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------------------------

N_THETA         <- 91
N_X             <- 401
GRID_OVERSAMPLE <- 3
N_GRID_MIN      <- 201
N_GRID_MAX      <- 801

# ------------------------------------------------------------------------------
# DESCRIPTORS
# ------------------------------------------------------------------------------

# Cached Schroedinger solutions.
.morse_soln <- NULL
ensure_morse <- function() {
  if (is.null(.morse_soln)) {
    cat("Solving Schroedinger for Morse...\n")
    .morse_soln <<- suppressMessages(solve_schroedinger(
      morse_V, MORSE_Q_MIN, MORSE_Q_MAX, MORSE_DQ,
      n_states=MORSE_N_STATES))
  }
  .morse_soln
}
.dw_soln <- NULL
ensure_dw <- function() {
  if (is.null(.dw_soln)) {
    cat("Solving Schroedinger for double well...\n")
    .dw_soln <<- suppressMessages(solve_schroedinger(
      double_well_V, DOUBLE_WELL_Q_MIN, DOUBLE_WELL_Q_MAX,
      DOUBLE_WELL_DQ, n_states=DOUBLE_WELL_N_STATES))
  }
  .dw_soln
}

SQUEEZED_R <- 0.5

squeezed_vacuum_descriptor <- list(
  kind = "eigenstate", name = "squeezed_vacuum", V = harmonic_V,
  n_target = 0,
  E_fn = function(n) 0.5 * cosh(2*SQUEEZED_R),
  psi_fn = function(n, q) {
    sigma_q <- exp(-SQUEEZED_R) / sqrt(2)
    (1 / (pi * sigma_q^2))^(1/4) * exp(-q^2 / (2 * sigma_q^2))
  },
  q_window = function(E) { Delta_q <- exp(-SQUEEZED_R)
                           list(q_lo=-3*Delta_q, q_hi=3*Delta_q) },
  p_window = function(E) { Delta_p <- exp(+SQUEEZED_R)
                           list(p_lo=-1.5*Delta_p, p_hi=1.5*Delta_p) },
  q_breaks_fn = function(E) { Delta_q <- exp(-SQUEEZED_R)
                              round(c(-Delta_q, Delta_q), 2) },
  p_breaks_fn = function(E) { Delta_p <- exp(+SQUEEZED_R)
                              round(c(-Delta_p, 0, Delta_p), 1) },
  psi_q_grid = seq(-15, 15, by=0.005)
)

harmonic_n1_descriptor <- list(
  kind = "eigenstate", name = "harmonic_n1", V = harmonic_V,
  n_target = 1,
  E_fn = function(n) n + 0.5,
  psi_fn = function(n, q) harmonic_psi(n, q),
  q_window = function(E) { qt <- sqrt(2*E); span <- 2*qt
                           list(q_lo=-qt-0.3*span/2, q_hi=qt+0.3*span/2) },
  p_window = function(E) { p_max <- sqrt(2*E)
                           list(p_lo=-1.3*p_max, p_hi=1.3*p_max) },
  q_breaks_fn = function(E) { qt <- sqrt(2*E); round(c(-qt, qt), 1) },
  p_breaks_fn = function(E) { p_max <- sqrt(2*E); round(c(-p_max, 0, p_max), 1) },
  psi_q_grid = seq(-25, 25, by=0.02)
)

morse_n8_descriptor <- list(
  kind = "eigenstate", name = "morse_n8", V = morse_V,
  n_target = 8,
  E_fn = function(n) ensure_morse()$energies[n + 1],
  psi_fn = function(n, q) {
    soln <- ensure_morse()
    psi_solver <- soln$psi_matrix[, n + 1]
    psi_q <- approx(soln$q_grid, psi_solver, xout=q,
                    rule=2, yleft=0, yright=0)$y
    psi_q[is.na(psi_q)] <- 0
    psi_q
  },
  q_window = function(E) {
    qm <- -log(1 + sqrt(E/De)) / alpha
    qp <- -log(1 - sqrt(E/De)) / alpha
    span <- qp - qm
    list(q_lo = qm - 0.15*span, q_hi = qp + 0.15*span)
  },
  p_window = function(E) { p_max <- sqrt(2*E)
                           list(p_lo=-1.3*p_max, p_hi=1.3*p_max) },
  q_breaks_fn = function(E) {
    qm <- -log(1 + sqrt(E/De)) / alpha
    qp <- -log(1 - sqrt(E/De)) / alpha
    round(c(qm, qp), 1)
  },
  p_breaks_fn = function(E) { p_max <- sqrt(2*E); round(c(-p_max, 0, p_max), 1) },
  psi_q_grid = NULL
)

double_well_n5_descriptor <- list(
  kind = "eigenstate", name = "double_well_n5", V = double_well_V,
  n_target = 5,
  E_fn = function(n) ensure_dw()$energies[n + 1],
  psi_fn = function(n, q) {
    soln <- ensure_dw()
    psi_solver <- soln$psi_matrix[, n + 1]
    psi_q <- approx(soln$q_grid, psi_solver, xout=q,
                    rule=2, yleft=0, yright=0)$y
    psi_q[is.na(psi_q)] <- 0
    psi_q
  },
  q_window = function(E) {
    roots <- polyroot(c(-E, 0, -mu2/2, 0, lambda/4))
    real_roots <- sort(Re(roots[abs(Im(roots)) < 1e-8]))
    q_lo <- min(real_roots); q_hi <- max(real_roots)
    span <- q_hi - q_lo
    list(q_lo = q_lo - 0.2*span, q_hi = q_hi + 0.2*span)
  },
  p_window = function(E) {
    V_min <- -double_well_barrier
    p_max <- sqrt(2 * (E - V_min))
    list(p_lo = -1.3*p_max, p_hi = 1.3*p_max)
  },
  q_breaks_fn = function(E) {
    roots <- polyroot(c(-E, 0, -mu2/2, 0, lambda/4))
    real_roots <- sort(Re(roots[abs(Im(roots)) < 1e-8]))
    round(c(min(real_roots), max(real_roots)), 1)
  },
  p_breaks_fn = function(E) {
    V_min <- -double_well_barrier
    p_max <- sqrt(2 * (E - V_min))
    round(c(-p_max, 0, p_max), 1)
  },
  psi_q_grid = NULL
)

fill_lazy_psi_grid <- function(d) {
  if (d$name == "morse_n8")        d$psi_q_grid <- ensure_morse()$q_grid
  if (d$name == "double_well_n5")  d$psi_q_grid <- ensure_dw()$q_grid
  d
}

make_cat_descriptor <- function(name, n_cats, variant) {
  list(
    kind = "cat", name = name,
    n_cats = n_cats, variant = variant, hbar = CAT_HBAR,
    q_lo = -CAT_Q_DISPLAY, q_hi = CAT_Q_DISPLAY,
    p_lo = -CAT_Q_DISPLAY, p_hi = CAT_Q_DISPLAY,
    custom_breaks_q = c(-CAT_P_MAX, 0, CAT_P_MAX),
    custom_breaks_p = c(-CAT_P_MAX, 0, CAT_P_MAX),
    label_format    = function(x) sprintf("%.0f", x)
  )
}

cat2_descriptor    <- make_cat_descriptor("cat_2",        2, "diag")
cat3_descriptor    <- make_cat_descriptor("cat_3",        3, "diag")
cat4sq_descriptor  <- make_cat_descriptor("cat_4_square", 4, "diag")
compass_descriptor <- make_cat_descriptor("cat_compass",  4, "axis")

ALL_DESCRIPTORS <- list(
  squeezed_vacuum_descriptor,
  harmonic_n1_descriptor,
  morse_n8_descriptor,
  double_well_n5_descriptor,
  cat2_descriptor,
  cat3_descriptor,
  cat4sq_descriptor,
  compass_descriptor
)

# ------------------------------------------------------------------------------
# DISPATCH
# ------------------------------------------------------------------------------

build_state_for_descriptor <- function(d) {
  if (d$kind == "eigenstate") {
    d <- fill_lazy_psi_grid(d)
    build_eigenstate_state(d)
  } else if (d$kind == "cat") {
    build_cat_state(d)
  } else {
    stop(sprintf("Unknown descriptor kind: %s", d$kind))
  }
}

# ------------------------------------------------------------------------------
# ROTATED KERNEL, RADON, FBP
# ------------------------------------------------------------------------------

rotated_kernel_matrix <- function(q_grid, p_grid, theta,
                                   Delta_q, Delta_p, hbar=1.0) {
  Delta_x <- sqrt(Delta_q^2 * cos(theta)^2 + Delta_p^2 * sin(theta)^2)
  Delta_y <- sqrt(Delta_q^2 * sin(theta)^2 + Delta_p^2 * cos(theta)^2)
  delta_x <- hbar / Delta_y

  q_mid <- (min(q_grid) + max(q_grid)) / 2
  p_mid <- (min(p_grid) + max(p_grid)) / 2

  outer(q_grid, p_grid, FUN = function(q, p) {
    qq <- q - q_mid; pp <- p - p_mid
    x  <-  qq * cos(theta) + pp * sin(theta)
    y  <- -qq * sin(theta) + pp * cos(theta)
    (1/pi) * exp(-(x/delta_x)^2 - (y/Delta_y)^2)
  })
}

bilinear_interp <- function(f_mat, q_grid, p_grid, q_query, p_query) {
  nq <- length(q_grid); np <- length(p_grid)
  dq <- diff(q_grid)[1]; dp <- diff(p_grid)[1]
  iq <- (q_query - q_grid[1]) / dq + 1
  ip <- (p_query - p_grid[1]) / dp + 1
  iq_lo <- floor(iq); iq_hi <- iq_lo + 1
  ip_lo <- floor(ip); ip_hi <- ip_lo + 1
  valid <- iq_lo >= 1 & iq_hi <= nq & ip_lo >= 1 & ip_hi <= np
  out   <- rep(0, length(q_query))
  if (!any(valid)) return(out)
  iq_lo_v <- iq_lo[valid]; iq_hi_v <- iq_hi[valid]
  ip_lo_v <- ip_lo[valid]; ip_hi_v <- ip_hi[valid]
  fq <- iq[valid] - iq_lo_v; fp <- ip[valid] - ip_lo_v
  v00 <- f_mat[cbind(iq_lo_v, ip_lo_v)]
  v10 <- f_mat[cbind(iq_hi_v, ip_lo_v)]
  v01 <- f_mat[cbind(iq_lo_v, ip_hi_v)]
  v11 <- f_mat[cbind(iq_hi_v, ip_hi_v)]
  out[valid] <- (1-fq)*(1-fp)*v00 + fq*(1-fp)*v10 +
                (1-fq)*fp*v01 + fq*fp*v11
  out
}

project_at_angle <- function(f_mat, q_grid, p_grid, theta, x_grid, y_grid) {
  ct <- cos(theta); st <- sin(theta)
  dy <- diff(y_grid)[1]
  result <- numeric(length(x_grid))
  for (i in seq_along(x_grid)) {
    x_i    <- x_grid[i]
    q_line <- x_i * ct - y_grid * st
    p_line <- x_i * st + y_grid * ct
    vals   <- bilinear_interp(f_mat, q_grid, p_grid, q_line, p_line)
    result[i] <- sum(vals) * dy
  }
  result
}

ramp_filter <- function(proj, dx) {
  n <- length(proj)
  m <- 2^ceiling(log2(2*n))
  pad <- c(proj, rep(0, m - n))
  fft_p <- fft(pad)
  freqs <- c(0:(m/2), -(m/2 - 1):-1) / (m * dx)
  filt  <- fft_p * abs(freqs)
  filt_t <- Re(fft(filt, inverse=TRUE)) / m
  filt_t[1:n]
}

inverse_radon <- function(sinogram, x_grid, theta_grid, q_grid, p_grid) {
  n_theta <- length(theta_grid)
  dx      <- diff(x_grid)[1]
  recon   <- matrix(0, nrow=length(q_grid), ncol=length(p_grid))
  for (i in seq_len(n_theta)) {
    theta <- theta_grid[i]
    proj  <- sinogram[, i]
    filt  <- ramp_filter(proj, dx)
    ct <- cos(theta); st <- sin(theta)
    q_mat <- matrix(q_grid, nrow=length(q_grid), ncol=length(p_grid))
    p_mat <- matrix(p_grid, nrow=length(q_grid), ncol=length(p_grid), byrow=TRUE)
    x_query <- q_mat * ct + p_mat * st
    iq      <- (x_query - x_grid[1]) / dx + 1
    iq_lo   <- floor(iq); iq_hi <- iq_lo + 1
    valid   <- iq_lo >= 1 & iq_hi <= length(x_grid)
    fq      <- iq - iq_lo
    contrib <- matrix(0, nrow=length(q_grid), ncol=length(p_grid))
    contrib[valid] <- (1 - fq[valid]) * filt[iq_lo[valid]] +
                            fq[valid]  * filt[iq_hi[valid]]
    recon <- recon + contrib
  }
  recon * pi / n_theta
}

bilinear_resample <- function(f_mat, q_src, p_src, q_dst, p_dst) {
  nq <- length(q_dst); np <- length(p_dst)
  out <- matrix(0, nrow=nq, ncol=np)
  for (j in seq_len(np)) {
    p_j <- p_dst[j]
    for (i in seq_len(nq)) {
      q_i <- q_dst[i]
      out[i, j] <- bilinear_interp(f_mat, q_src, p_src, q_i, p_j)
    }
  }
  out
}

# ------------------------------------------------------------------------------
# PER-STATE TOMOGRAPHY
# ------------------------------------------------------------------------------

run_tomography <- function(ps) {
  cat(sprintf("  Running tomography for %s...\n", ps$name))

  # Adaptive grid sizing.
  n_q <- min(N_GRID_MAX,
             max(N_GRID_MIN,
                 ceiling((ps$q_hi - ps$q_lo) * GRID_OVERSAMPLE / ps$rs$delta_q)))
  n_p <- min(N_GRID_MAX,
             max(N_GRID_MIN,
                 ceiling((ps$p_hi - ps$p_lo) * GRID_OVERSAMPLE / ps$rs$delta_p)))
  q_grid <- seq(ps$q_lo, ps$q_hi, length.out=n_q)
  p_grid <- seq(ps$p_lo, ps$p_hi, length.out=n_p)
  dq     <- diff(q_grid)[1]
  dp     <- diff(p_grid)[1]
  cat(sprintf("  Tomography grid: %d x %d  dq=%.4f  dp=%.4f\n",
              n_q, n_p, dq, dp))

  # Resample W from the build_*_state grid (q_int, p_int) onto the
  # tomography grid.
  cat("  Resampling W onto tomography grid...\n")
  W_on_grid <- bilinear_resample(ps$state$W_matrix,
                                  ps$state$q_int, ps$state$p_int,
                                  q_grid, p_grid)

  half_diag <- sqrt(((ps$q_hi - ps$q_lo)/2)^2 +
                      ((ps$p_hi - ps$p_lo)/2)^2)
  x_grid    <- seq(-half_diag, half_diag, length.out=N_X)
  y_grid    <- seq(-half_diag, half_diag, length.out=max(n_q, n_p))

  theta_grid <- seq(0, pi, length.out=N_THETA + 1)[-(N_THETA + 1)]
  sinogram   <- matrix(0, nrow=N_X, ncol=N_THETA)

  cat(sprintf("  Sweeping %d angles...\n", N_THETA))
  for (i in seq_len(N_THETA)) {
    if (i %% 10 == 0) cat(sprintf("    theta = %d/%d\n", i, N_THETA))
    theta   <- theta_grid[i]
    K_mat   <- rotated_kernel_matrix(q_grid, p_grid, theta,
                                      ps$rs$Delta_q, ps$rs$Delta_p,
                                      hbar=ps$hbar)
    P_theta <- fft_convolve_2d(W_on_grid, K_mat, dq, dp)$P_mat
    sinogram[, i] <- project_at_angle(P_theta, q_grid, p_grid,
                                       theta, x_grid, y_grid)
  }

  cat("  Inverse Radon...\n")
  tilde_W <- inverse_radon(sinogram, x_grid, theta_grid, q_grid, p_grid)

  rho_q <- rowSums(tilde_W) * dp
  rho_p <- colSums(tilde_W) * dq

  exact_q <- exact_position_marginal(ps$psi_vec, ps$psi_q_grid, q_grid)
  exact_p <- exact_momentum_marginal(ps$psi_vec, ps$psi_q_grid,
                                      p_grid, hbar=ps$hbar)

  list(
    q_grid_tomo = q_grid,
    p_grid_tomo = p_grid,
    x_grid      = x_grid,
    theta_grid  = theta_grid,
    W_tomo      = W_on_grid,
    sinogram    = sinogram,
    tilde_W     = tilde_W,
    rho_q       = rho_q,
    rho_p       = rho_p,
    exact_q     = exact_q,
    exact_p     = exact_p
  )
}

# ------------------------------------------------------------------------------
# DRIVE
# ------------------------------------------------------------------------------

cat("Computing tomographic reconstruction across eight states...\n")
cat(sprintf("Angular sweep: %d angles. Sinogram radial samples: %d.\n",
            N_THETA, N_X))

results <- list()
for (d in ALL_DESCRIPTORS) {
  ps     <- build_state_for_descriptor(d)
  tomo   <- run_tomography(ps)
  bundle <- c(ps, tomo)
  results[[ps$name]] <- bundle
}

dir.create(here("data"), showWarnings=FALSE)
attr(results, "version") <- TOMOGRAPHY_CACHE_VERSION
saveRDS(results, OUTPUT_RDS)
cat(sprintf("\nWrote %s\n", OUTPUT_RDS))
