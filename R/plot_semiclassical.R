# ==============================================================================
# plot_semiclassical.R
# Semiclassical resolution of WKB caustics across three systems.
# Each row shows one system at one chosen action capacity:
#   row 1: harmonic   n=0   (ground state, isotropic Fermi blob)
#   row 2: Morse      n=8   (anharmonic mid-state, asymmetric horseshoe orbit)
#   row 3: double-well n=1  (sub-barrier lower-doublet partner, four turning
#                            points, two disconnected orbital loops; aligned
#                            with the Wigner figure's Schroedinger n=1 row at
#                            A_BS/A_0 = 1.97593)
#
# Three columns:
#   left:    classical orbit at E_n with the symplectic quantum of action
#            (Fermi blob + conjugate quantum blobs) overlaid.
#   center:  oscillating WKB density |psi_WKB(q)|^2 with infinity-arrows
#            at every classical turning point.
#   right:   symplectic resolution rho_{delta q}(q) (the kernel-convolved
#            classical microcanonical density, marginalized over p) with
#            the Airy uniform density overlaid as the prior-art comparator.
#
# SEMICLASSICAL UNIVERSE -- WHAT THIS PIPELINE TOUCHES (and does NOT touch):
#
#   This pipeline lives entirely in the semiclassical universe:
#     * potential V(q), classical orbit at E, classical turning points
#     * orbit time-averaged second moments -> orbit covariance
#     * action capacity A = pi * Delta_q * Delta_p (orbit-derived)
#     * kernel widths delta_q = hbar/Delta_p, delta_p = hbar/Delta_q
#     * oscillating WKB density (built from classical action S(q), not psi)
#     * symplectic kernel G_{delta q} convolved with classical
#       microcanonical density W_cl propto delta(H - E)
#     * Airy uniform construction (built from classical action and turning
#       points; does NOT use psi)
#     * 1D position marginal rho_{delta q}(q)
#
#   This pipeline does NOT touch:
#     * Schroedinger wavefunctions psi_n(q)
#     * Robertson-Schroedinger covariance of any psi
#     * Wigner functions of any kind
#     * Husimi kernel or Husimi Q-function
#
#   The single shared object across the two universes (semiclassical and
#   quantum) is the symplectic kernel G_{delta q} itself: a Gaussian of
#   phase-space area h/2. Its widths in this file come from the classical
#   orbit's covariance -- NOT from any quantum state.
#
# Author: Brian S. Mulloy
# ==============================================================================

library(here)
library(data.table)
library(patchwork)

# Plot helpers (shared across both universes; no quantum-state machinery here).
source(here("R", "plot_tools.R"))

# Semiclassical-universe modules ONLY. No schroedinger_solver, no
# wigner_density, no husimi_kernel, no quantum_tools.
source(here("R", "harmonic_system.R"))           # harmonic_V, harmonic_turning_points,
# harmonic_E_at_action_capacity
source(here("R", "morse_system.R"))              # morse_V, morse_turning_points,
# morse_E_BS, morse_E_at_action_capacity
source(here("R", "double_well_system.R"))        # double_well_V,
# double_well_turning_points,
# double_well_E_at_action_capacity
source(here("R", "classical_action_tools.R"))    # orbit_covariance,
# classical_trajectory
source(here("R", "semiclassical_density.R"))     # build_semiclassical_state,
# wkb_phase_space_lift,
# oscillating_wkb_density,
# extract_p0_cross_section
source(here("R", "symplectic_kernel.R"))         # G_delta_q_kernel_matrix,
# symplectic_marginal_density,
# symplectic_overlay_layers
source(here("R", "airy_uniform.R"))              # airy_uniform_density

latex_font  <- "CMU Serif"
dir_figures <- here("figures")
if (!dir.exists(dir_figures)) dir.create(dir_figures, recursive=TRUE)
file_output_pdf <- file.path(dir_figures, "semiclassical.pdf")

# Vertical-fill fraction for ribbon plots.
GOLDEN_FILL <- 0.6

# ------------------------------------------------------------------------------
# MULTI-SEGMENT SEMICLASSICAL STATE BUILDER
#
# build_semiclassical_state() in semiclassical_density.R takes a single
# (q_minus, q_plus) pair. For sub-barrier double-well states with four
# turning points, the WKB density is supported on two disconnected orbital
# loops; we build it on each loop separately and sum, then run the same
# downstream symplectic convolution.
# ------------------------------------------------------------------------------

build_semiclassical_state_segments <- function(E_n, V_fn,
                                               q_lo, q_hi, p_lo, p_hi,
                                               q_display,
                                               turning_points,
                                               epsilon=NULL,
                                               nq_int=801, np_int=601) {
  if (is.null(epsilon)) epsilon <- 0.05 * abs(E_n)

  q_int  <- seq(q_lo, q_hi, length.out=nq_int)
  p_int  <- seq(p_lo, p_hi, length.out=np_int)
  dq_int <- diff(q_int)[1]
  dp_int <- diff(p_int)[1]

  cat(sprintf("    Building multi-segment WKB lift on %d x %d grid (epsilon=%.4f)...\n",
              nq_int, np_int, epsilon))

  W_mat         <- matrix(0, nrow=nq_int, ncol=np_int)
  wkb_osc_total <- numeric(length(q_display))
  n_segments    <- length(turning_points) / 2

  for (k in seq(1, length(turning_points), by=2)) {
    q_minus_k <- turning_points[k]
    q_plus_k  <- turning_points[k+1]
    cat(sprintf("      segment %d: q in [%.3f, %.3f]\n",
                (k+1)/2, q_minus_k, q_plus_k))

    W_seg         <- wkb_phase_space_lift(q_int, p_int, E_n, V_fn,
                                          q_minus_k, q_plus_k, epsilon)
    W_mat         <- W_mat + W_seg
    wkb_osc_seg   <- oscillating_wkb_density(q_display, E_n, V_fn,
                                             q_minus_k, q_plus_k)
    wkb_osc_total <- wkb_osc_total + wkb_osc_seg
  }

  Z <- sum(W_mat) * dq_int * dp_int
  if (Z > 0) W_mat <- W_mat / Z
  cat(sprintf("    Multi-segment shell norm: %.6f -> renormalized to 1.\n", Z))

  if (n_segments > 1) wkb_osc_total <- wkb_osc_total / n_segments

  W_cross <- extract_p0_cross_section(W_mat, q_int, p_int, q_display)

  heatmap_dt <- as.data.table(expand.grid(q=q_int, p=p_int))
  heatmap_dt[, w := as.vector(W_mat)]
  max_w <- max(heatmap_dt$w, na.rm=TRUE)
  if (max_w > 0) {
    heatmap_dt[, w_plot := sqrt(pmax(w/max_w, 0))]
  } else {
    heatmap_dt[, w_plot := 0]
  }

  list(
    q_int       = q_int,
    p_int       = p_int,
    dq_int      = dq_int,
    dp_int      = dp_int,
    W_matrix    = W_mat,
    W_cross     = W_cross,
    wkb_density = wkb_osc_total,
    heatmap_dt  = heatmap_dt,
    norm        = 1.0,
    epsilon     = epsilon
  )
}

# ------------------------------------------------------------------------------
# SYSTEM DESCRIPTORS
#
# Each descriptor packages the system-specific functions and parameters
# needed by the semiclassical row builder. The builder is otherwise
# generic: same pipeline, same overlay logic, same column structure.
#
# Each descriptor exposes:
#   name           -- short identifier for console output
#   V              -- potential V(q)
#   turning_points -- function(E) -> numeric vector of turning q-values
#                     (length 2 for connected orbits, 4 for sub-barrier
#                     double-well)
#   E_BS           -- function(n) -> classical Bohr-Sommerfeld energy at
#                     quantum number n. Pure-classical solver only.
#   n_target       -- the chosen quantum number for this row
#   q_window       -- function(tp) -> list(q_lo, q_hi) for the display
#   p_window       -- function(E)  -> list(p_lo, p_hi) for the display
#   custom_breaks_q, custom_breaks_p, label_format -- axis breaks/format
# ------------------------------------------------------------------------------

# ---- Harmonic ----------------------------------------------------------------
harmonic_descriptor <- list(
  name = "harmonic",
  V    = harmonic_V,
  turning_points = function(E) {
    tp <- harmonic_turning_points(E)
    c(tp$q_minus, tp$q_plus)
  },

  # E_n = (n + 1/2) hbar omega; in dimensionless harmonic units (omega=hbar=1),
  # E_n = n + 1/2. This is BS quantization, not Schroedinger eigenvalues
  # (the two coincide for the harmonic oscillator -- a special property of
  # the harmonic system that does not hold elsewhere).
  E_BS = function(n) n + 0.5,

  n_target = 0,

  q_window = function(tp) {
    q_lo <- min(tp); q_hi <- max(tp)
    span <- q_hi - q_lo
    list(q_lo = q_lo - 0.3*span, q_hi = q_hi + 0.3*span)
  },
  p_window = function(E) {
    p_max <- sqrt(2 * E)
    list(p_lo = -1.3*p_max, p_hi = 1.3*p_max)
  },
  q_breaks_fn = function(tp) round(c(min(tp), max(tp)), 1),
  p_breaks_fn = function(E)  round(c(-sqrt(2*E), 0, sqrt(2*E)), 1)
)

# ---- Morse -------------------------------------------------------------------
morse_descriptor <- list(
  name = "morse",
  V    = morse_V,
  turning_points = function(E) {
    tp <- morse_turning_points(E)
    c(tp$q_minus, tp$q_plus)
  },
  E_BS = morse_E_BS,
  n_target = 8,

  q_window = function(tp) {
    q_lo <- min(tp); q_hi <- max(tp)
    span <- q_hi - q_lo
    list(q_lo = q_lo - 0.15*span, q_hi = q_hi + 0.15*span)
  },
  p_window = function(E) {
    p_max <- sqrt(2 * E)
    list(p_lo = -1.3*p_max, p_hi = 1.3*p_max)
  },
  q_breaks_fn = function(tp) round(c(min(tp), max(tp)), 1),
  p_breaks_fn = function(E)  round(c(-sqrt(2*E), 0, sqrt(2*E)), 1)
)

# ---- Double well -------------------------------------------------------------
#
# Classical Bohr-Sommerfeld quantization on the double well, computed via
# the textbook closed-orbit integral A_BS = oint p dq returned by
# classical_action() in classical_action_tools.R. NOT the Fermi-blob
# moment-product returned by orbit_covariance(); those are different
# objects (they coincide only for harmonic orbits).
#
# Sub-barrier states (4 turning points) form tunneling doublets: per-well
# BS quantization gives single-well action (n_well + 1/2) h on each well,
# combined orbit area 2 * (n_well + 1/2) h = (2 n_well + 1) * 2 A_0 in
# units of A_0 = h/2. Above-barrier states have a single connected
# (peanut) orbit; their BS action is the full peanut integral.
#
# For the lowest sub-barrier doublet, classical action A_BS/A_0 ~ 2
# gives E ~ -2.64 with 4 turning points at +/- 1.29, +/- 2.52 -- the
# canonical doublet visual with two well-localized classical loops.
#
# This row is aligned to Schroedinger n=1 in the Wigner figure: the
# antisymmetric lower-doublet partner has a node at q=0 (psi_1(0)=0),
# making the title's "Schroedinger nodes" pathology visible there. The
# classical action at E_1 = -2.65175 is A_BS/A_0 = 1.97593, and the
# semiclassical row uses that exact value so both figures sit at the
# same physical energy.
double_well_E_at_classical_action <- function(target_A_over_A0, tol=1e-4) {
  if (target_A_over_A0 <= 0) stop("target_A_over_A0 must be positive")

  V_min <- -mu2^2 / (4 * lambda)   # = -4

  A_at <- function(E) {
    tp_obj <- double_well_turning_points(E)
    rt     <- tp_obj$roots
    if (length(rt) %% 2 != 0) return(NA)
    classical_action(double_well_V, E, rt)
  }

  # Search range: just above the well bottom up to a value well above the
  # barrier. classical_action is monotonically increasing in E (across both
  # 4-tp and 2-tp regimes), so bisection converges.
  E_lo <- V_min + 1e-3
  E_hi <- 50.0
  A_lo <- A_at(E_lo)
  A_hi <- A_at(E_hi)
  if (is.na(A_lo) || is.na(A_hi))
    stop("double_well_E_at_classical_action: classical_action undefined at search bounds")
  if (target_A_over_A0 < A_lo || target_A_over_A0 > A_hi)
    stop(sprintf("target A/A_0=%.3f outside reachable range [%.3f, %.3f]",
                 target_A_over_A0, A_lo, A_hi))

  for (iter in 1:100) {
    E_mid <- 0.5 * (E_lo + E_hi)
    A_mid <- A_at(E_mid)
    if (is.na(A_mid)) break
    if (abs(A_mid - target_A_over_A0) < tol) return(E_mid)
    if (A_mid < target_A_over_A0) E_lo <- E_mid else E_hi <- E_mid
  }
  warning(sprintf("double_well_E_at_classical_action: bisection did not converge to %.1e",
                  tol))
  0.5 * (E_lo + E_hi)
}

# Target action capacity aligned with Schroedinger n=1 in the Wigner
# figure. Computed by evaluating classical_action at E_1 = -2.65175
# (the antisymmetric lower-doublet partner's eigenvalue). This makes
# the semiclassical and Wigner figures look at the same physical state.
DOUBLE_WELL_TARGET_A_OVER_A0 <- 1.97593

double_well_descriptor <- list(
  name = "double_well",
  V    = double_well_V,
  turning_points = function(E) double_well_turning_points(E)$roots,
  # Energy is set by a target action capacity (continuous), not an integer
  # quantum number. The chosen target aligns with Schroedinger n=1 in the
  # Wigner figure so the two figures sit at the same physical energy.
  E_BS = function(n) double_well_E_at_classical_action(DOUBLE_WELL_TARGET_A_OVER_A0),
  # n_target is informational; the row's content is set by the action
  # capacity above, not by n. Kept as 1 to record alignment with the
  # Wigner figure's Schroedinger n=1 row.
  n_target = 1,

  q_window = function(tp) {
    q_lo <- min(tp); q_hi <- max(tp)
    span <- q_hi - q_lo
    list(q_lo = q_lo - 0.2*span, q_hi = q_hi + 0.2*span)
  },
  p_window = function(E) {
    V_min <- -double_well_barrier
    p_max <- sqrt(2 * (E - V_min))
    list(p_lo = -1.3*p_max, p_hi = 1.3*p_max)
  },
  q_breaks_fn = function(tp) round(c(min(tp), max(tp)), 1),
  p_breaks_fn = function(E)  {
    V_min <- -double_well_barrier
    p_max <- sqrt(2 * (E - V_min))
    round(c(-p_max, 0, p_max), 1)
  }
)

# ------------------------------------------------------------------------------
# GENERIC SEMICLASSICAL ROW BUILDER
#
# Takes one system descriptor, builds the three-panel semiclassical row.
# Dispatch on number of turning points (2 or 4) selects single-segment
# vs multi-segment semiclassical-state construction.
# ------------------------------------------------------------------------------

build_semiclassical_row <- function(descriptor, base_font="") {
  n   <- descriptor$n_target
  E_n <- descriptor$E_BS(n)
  tp  <- descriptor$turning_points(E_n)

  cat(sprintf("\n== %s | n=%d | E_n=%.4f | tps=[%s] ==\n",
              descriptor$name, n, E_n,
              paste(sprintf("%.3f", tp), collapse=", ")))

  is_multi_segment <- length(tp) == 4

  # Orbit covariance gives Delta_q, Delta_p (and via the reciprocal scales,
  # delta_q, delta_p). All from the classical orbit. Pure semiclassical.
  cov <- orbit_covariance(descriptor$V, E_n, tp)
  cat(sprintf("  A_orbit/A0=%.4f | <q>=%.4f | Delta_q=%.3f Delta_p=%.3f\n",
              cov$A_over_A0, cov$q_mean, cov$Delta_q, cov$Delta_p))
  cat(sprintf("  delta_q=%.3f delta_p=%.3f\n", cov$delta_q, cov$delta_p))

  # Display windows.
  qw <- descriptor$q_window(tp)
  pw <- descriptor$p_window(E_n)
  q_lo <- qw$q_lo; q_hi <- qw$q_hi
  p_lo <- pw$p_lo; p_hi <- pw$p_hi

  custom_breaks_q <- descriptor$q_breaks_fn(tp)
  custom_breaks_p <- descriptor$p_breaks_fn(E_n)
  label_format    <- function(x) sprintf("%.1f", x)
  q_display       <- seq(q_lo, q_hi, length.out=500)

  # Build the semiclassical state. Single-orbit (2 tps) uses the standard
  # builder; sub-barrier (4 tps) uses the multi-segment builder above.
  if (is_multi_segment) {
    state <- build_semiclassical_state_segments(
      E_n=E_n, V_fn=descriptor$V,
      q_lo=q_lo, q_hi=q_hi, p_lo=p_lo, p_hi=p_hi,
      q_display=q_display, turning_points=tp)
  } else {
    state <- build_semiclassical_state(
      E_n=E_n, V_fn=descriptor$V,
      q_lo=q_lo, q_hi=q_hi, p_lo=p_lo, p_hi=p_hi,
      q_display=q_display,
      q_minus=tp[1], q_plus=tp[2])
  }

  # Convolve with the symplectic kernel (orbit-derived widths) and
  # marginalize over p to get the 1D position density rho_{delta q}(q).
  symplectic_kernel_for_state <- function(qg, pg) {
    G_delta_q_kernel_matrix(qg, pg, cov$Delta_q, cov$Delta_p)
  }
  rho_sympl <- symplectic_marginal_density(state, symplectic_kernel_for_state,
                                           q_display)

  # Airy uniform comparator. airy_uniform_density() handles both 2-tp and
  # 4-tp cases via its existing per-segment-Langer-summed loop.
  rho_airy <- airy_uniform_density(q_display, E_n, descriptor$V, tp)

  # Y-scaling.
  rho_peak  <- max(rho_sympl, na.rm=TRUE)
  if (!is.finite(rho_peak) || rho_peak == 0) rho_peak <- 1
  y_lim_rho <- rho_peak / (GOLDEN_FILL + 0.2)

  # Center-column y-limit: drop near-tp Airy spikes from the peak detection
  # (sub-barrier states would otherwise have the inner-tp WKB peaks dominate).
  span <- q_hi - q_lo
  inside_pad <- 0.02 * span
  near_tp <- rep(FALSE, length(q_display))
  for (qt in tp) near_tp <- near_tp | (abs(q_display - qt) < inside_pad)
  finite_inside <- state$wkb_density[!near_tp &
                                       is.finite(state$wkb_density) &
                                       state$wkb_density > 0]
  if (length(finite_inside) > 0) {
    osc_peak     <- max(finite_inside)
    y_lim_caustic <- osc_peak / GOLDEN_FILL
  } else {
    y_lim_caustic <- 1.0
  }

  dt_caustic   <- data.table(q=q_display, wkb_density=state$wkb_density)
  dt_rho_sympl <- data.table(q=q_display, rho_sympl=rho_sympl)

  # QoA overlay (orbit-derived widths).
  overlay_layers <- symplectic_overlay_layers(cov$Delta_q, cov$Delta_p,
                                              q_center=cov$q_mean)
  df_traj <- classical_trajectory(descriptor$V, E_n, tp)

  # Airy overlay (semiclassical comparator, prior art). Drawn on all
  # rows including the sub-barrier doublet (multi-orbit) row. Per-well
  # Langer/Airy uniform approximation is the textbook treatment for
  # WKB caustics: matched to harmonic-oscillator wavefunctions near
  # each well minimum, with exponential decay through the barrier. At
  # the doublet's inner turning points the per-well construction
  # spikes -- this is a genuine failure of the textbook-level method,
  # not an artifact of the implementation. More sophisticated uniform
  # approximations exist (parabolic-cylinder matching at the barrier
  # region, Mathieu-function uniform approximations) but require
  # state-specific machinery beyond Langer's universal turning-point
  # patch. We show Langer as-is because it is the prior art whose
  # limits define this paper's contribution; the manuscript's
  # discussion of the figure names the research-level alternatives
  # and explains the visual clipping.
  airy_overlay <- list(
    list(
      data       = data.frame(q = q_display, rho = rho_airy),
      color      = "gray30",
      linewidth  = 0.35
    )
  )

  list(
    plot_classical_orbit_phase_space(
      df_traj,
      q_lim=c(q_lo,q_hi), p_lim=c(p_lo,p_hi),
      custom_breaks_q=custom_breaks_q,
      custom_breaks_p=custom_breaks_p,
      label_format=label_format, base_font=base_font,
      overlay_layers=overlay_layers,
      orbit_color="gray45",
      orbit_linewidth=0.4),
    plot_wkb_caustic_cross_section(
      dt_caustic, q_lim=c(q_lo,q_hi), y_lim=y_lim_caustic,
      custom_breaks=custom_breaks_q,
      label_format=label_format,
      q_minus=NULL, q_plus=NULL,  # auto-detected from data
      base_font=base_font),
    plot_semiclassical_resolution(
      dt_rho_sympl, q_lim=c(q_lo,q_hi), y_lim=y_lim_rho,
      custom_breaks=custom_breaks_q,
      label_format=label_format, base_font=base_font,
      overlays=airy_overlay)
  )
}

# ------------------------------------------------------------------------------
# DRIVE: build all three rows, assemble, save.
# ------------------------------------------------------------------------------

cat("Computing semiclassical figure (3 systems x 3 panels)...\n")

descriptors <- list(harmonic_descriptor,
                    morse_descriptor,
                    double_well_descriptor)

rows <- lapply(descriptors,
               function(d) build_semiclassical_row(d, base_font=latex_font))

p_final <- assemble_grid_unlabeled(rows,
                                   COLUMN_TITLE_CENTER_SEMICLASSICAL,
                                   COLUMN_TITLE_RIGHT_SYMPLECTIC,
                                   base_font=latex_font)

save_figure(p_final, file_output_pdf, length(descriptors))
