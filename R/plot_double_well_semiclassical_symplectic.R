# ==============================================================================
# plot_double_well_semiclassical_symplectic.R
# Symplectic resolution of WKB caustics — symmetric double well.
#
# Three rows indexed by quantum number:
#   n=0  lowest sub-barrier tunneling-doublet partner (4 turning points,
#        2 disconnected classical orbits)
#   n=2  upper sub-barrier doublet partner (4 turning points, 2 orbits)
#   n=7  above-barrier connected orbit (2 turning points, 1 peanut-shaped
#        orbit threading through the barrier region)
#
# Three columns:
#   Left:    classical orbit at energy E_n (1D curve in 2D phase space)
#            with QoA overlay. classical_trajectory() handles multi-
#            segment orbits via NA-row separators.
#   Center:  oscillating WKB density |psi_WKB(q)|^2. For sub-barrier
#            states the density is supported on two disconnected wells
#            (4 turning points -> 4 infinity arrows); for above-barrier
#            it's one connected support region with 2 turning points.
#   Right:   1D position density rho_{delta q}(q), the marginal of
#            (W_cl * G_{delta q}). Airy uniform density rho_Airy(q)
#            overlaid in dashed gray as the prior-art semiclassical
#            comparator. airy_uniform_density() handles both sub-barrier
#            (4 turning points -> two per-well Langer densities summed)
#            and above-barrier (2 turning points -> single Langer
#            density) regimes via its existing per-orbit-segment loop.
#            Same overlay style as the Morse semiclassical figure so
#            the comparator reads consistently across both figures.
#
# Multi-segment pipeline (sub-barrier only): we cannot use
# build_semiclassical_state() directly because it takes a single
# (q_minus, q_plus) pair. Instead we call its primitives —
# wkb_phase_space_lift() and oscillating_wkb_density() — on each
# orbital segment separately and sum the results. The 2D phase-space
# lift W_mat is then a sum of two single-orbit lifts; the 1D
# oscillating WKB density is the sum of two single-orbit densities;
# both go through the symplectic convolution as in the Morse case.
# ==============================================================================

library(here)
library(patchwork)

source(here("R", "plot_tools.R"))
source(here("R", "double_well_system.R"))
source(here("R", "classical_action_tools.R"))   # orbit_covariance,
# classical_trajectory
source(here("R", "semiclassical_density.R"))    # build_semiclassical_state,
# wkb_phase_space_lift,
# oscillating_wkb_density,
# extract_p0_cross_section
source(here("R", "symplectic_kernel.R"))        # G_delta_q_kernel_matrix,
# symplectic_marginal_density,
# symplectic_overlay_layers
source(here("R", "airy_uniform.R"))             # airy_uniform_density
# (prior-art semiclassical
# comparator for the
# right-column overlay)
source(here("R", "schroedinger_solver.R"))      # solve_schroedinger
# (used only for E_n, to keep
# the heatmap and orbit
# visually exactly co-energetic
# with the Wigner figure;
# the WKB and convolution
# pipelines are otherwise
# Schroedinger-free)

latex_font      <- "CMU Serif"
dir_figures     <- here("figures")
if (!dir.exists(dir_figures)) dir.create(dir_figures, recursive=TRUE)
file_output_pdf <- file.path(dir_figures, "double_well_semiclassical_symplectic.pdf")

target_quantum_numbers <- c(0, 2, 7)

# Vertical-fill fraction for ribbons (matches Morse semiclassical figure).
GOLDEN_FILL <- 0.6

# Solve Schroedinger once just to get eigenvalues E_n. We do NOT use the
# wavefunctions here — only the energies, which serve as the action-
# quantization condition for the classical orbits. The Wigner figure
# uses the same energies, so the two figures' rows match exactly.
cat("Solving Schroedinger for double well (energies only)...\n")
dw_soln <- solve_schroedinger(double_well_V,
                              DOUBLE_WELL_Q_MIN, DOUBLE_WELL_Q_MAX,
                              DOUBLE_WELL_DQ,
                              n_states=DOUBLE_WELL_N_STATES)

# ------------------------------------------------------------------------------
# MULTI-SEGMENT WKB BUILDER (sub-barrier only)
#
# For 4-turning-point states we cannot use build_semiclassical_state()
# directly — it takes a single (q_minus, q_plus) pair. Instead we call
# wkb_phase_space_lift() and oscillating_wkb_density() on each orbital
# segment separately and sum. The result is a state bundle with the
# same field names that build_semiclassical_state() produces, so the
# downstream symplectic convolution and rendering work without changes.
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

  # Per-segment 2D phase-space lift, summed across segments. Each segment's
  # lift is normalized to integrate to 1; we sum and then renormalize so
  # the total integrates to 1 (so the symplectic convolution preserves
  # normalization downstream).
  W_mat <- matrix(0, nrow=nq_int, ncol=np_int)
  wkb_osc_total <- numeric(length(q_display))

  n_segments <- length(turning_points) / 2
  for (k in seq(1, length(turning_points), by=2)) {
    q_minus_k <- turning_points[k]
    q_plus_k  <- turning_points[k+1]
    cat(sprintf("      segment %d: q in [%.3f, %.3f]\n",
                (k+1)/2, q_minus_k, q_plus_k))

    W_seg <- wkb_phase_space_lift(q_int, p_int, E_n, V_fn,
                                  q_minus_k, q_plus_k, epsilon)
    W_mat <- W_mat + W_seg

    wkb_osc_seg <- oscillating_wkb_density(q_display, E_n, V_fn,
                                           q_minus_k, q_plus_k)
    wkb_osc_total <- wkb_osc_total + wkb_osc_seg
  }

  # Renormalize the 2D lift so the total integrates to 1
  Z <- sum(W_mat) * dq_int * dp_int
  if (Z > 0) W_mat <- W_mat / Z
  cat(sprintf("    Multi-segment shell norm before renormalization: %.6f (%d segments combined)\n",
              Z, n_segments))
  cat(sprintf("    Renormalized to 1.000000\n"))

  # Renormalize the 1D oscillating density similarly. The per-segment
  # oscillating_wkb_density is normalized to integrate to 1 over its own
  # support; summing 2 of them gives total integral 2, so we divide by
  # the number of segments to keep total integral = 1.
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

build_double_well_row <- function(n, base_font="") {

  E_n    <- dw_soln$energies[n + 1]
  tp_obj <- double_well_turning_points(E_n)
  tp     <- tp_obj$roots
  is_sub_barrier <- length(tp) == 4

  if (is_sub_barrier) {
    cat(sprintf("\n== n=%d | E_n=%.4f | sub-barrier | tps=[%.3f, %.3f, %.3f, %.3f] ==\n",
                n, E_n, tp[1], tp[2], tp[3], tp[4]))
  } else {
    cat(sprintf("\n== n=%d | E_n=%.4f | above-barrier | tps=[%.3f, %.3f] ==\n",
                n, E_n, tp[1], tp[2]))
  }

  cov <- orbit_covariance(double_well_V, E_n, tp)
  cat(sprintf("  A_orbit/A0=%.4f | <q>=%.4f | Delta_q=%.3f Delta_p=%.3f\n",
              cov$A_over_A0, cov$q_mean, cov$Delta_q, cov$Delta_p))

  # Display window — same convention as the double-well Wigner figure so
  # rows align across the two figures.
  q_outer_lo <- min(tp)
  q_outer_hi <- max(tp)
  q_span     <- q_outer_hi - q_outer_lo
  q_pad      <- q_span * 0.2
  q_lo       <- q_outer_lo - q_pad
  q_hi       <- q_outer_hi + q_pad

  V_min  <- -double_well_barrier
  p_max  <- sqrt(2 * (E_n - V_min))
  p_lo   <- -1.3 * p_max
  p_hi   <-  1.3 * p_max

  custom_breaks_q <- round(c(q_outer_lo, q_outer_hi), 1)
  custom_breaks_p <- round(c(-p_max, 0, p_max), 1)
  label_format    <- function(x) sprintf("%.1f", x)
  q_display       <- seq(q_lo, q_hi, length.out=500)

  # Build the semiclassical state. Above-barrier (2 tps): use the standard
  # build_semiclassical_state() with q_minus, q_plus pair. Sub-barrier
  # (4 tps): use the multi-segment builder defined above.
  if (is_sub_barrier) {
    state <- build_semiclassical_state_segments(
      E_n=E_n, V_fn=double_well_V,
      q_lo=q_lo, q_hi=q_hi, p_lo=p_lo, p_hi=p_hi,
      q_display=q_display,
      turning_points=tp)
  } else {
    state <- build_semiclassical_state(
      E_n=E_n, V_fn=double_well_V,
      q_lo=q_lo, q_hi=q_hi, p_lo=p_lo, p_hi=p_hi,
      q_display=q_display,
      q_minus=tp[1], q_plus=tp[2])
  }

  # Symplectic kernel and marginalized density. Same kernel-agnostic
  # helper as the Morse figure.
  symplectic_kernel_for_state <- function(q_grid, p_grid) {
    G_delta_q_kernel_matrix(q_grid, p_grid, cov$Delta_q, cov$Delta_p)
  }
  rho_sympl <- symplectic_marginal_density(state, symplectic_kernel_for_state,
                                           q_display)

  # Airy uniform density on the display grid. airy_uniform_density()
  # already handles multi-orbit cases by looping over orbital segments
  # and summing per-segment Langer densities — so it works directly for
  # both sub-barrier (4 turning points, 2 wells, summed) and above-
  # barrier (2 turning points, single orbit) regimes.
  #
  # For sub-barrier states the per-well Langer construction has tall
  # spikes at the inner turning points: the textbook Langer formula
  # places an Airy lobe peak exactly where the well's classical orbit
  # ends, which for the doublet is inside the barrier rather than at
  # the well center. These peaks are real Langer features (not
  # numerical artifacts), and they exemplify why the Langer/Miller
  # uniform-Airy construction does not handle two-well coherent states
  # cleanly. We display Airy as-is and let the inner-tp spikes clip
  # out the top of the panel; the y-axis is set by the symplectic
  # peak so the symplectic result stays readable. The clipping is
  # itself the visual statement.
  rho_airy <- airy_uniform_density(q_display, E_n, double_well_V, tp)

  # Y-scaling: driven by the symplectic peak, NOT the Airy peak. Airy
  # spikes at inner turning points (sub-barrier doublets) will clip out
  # the top of the panel — visually honest about Langer's behavior at
  # interior caustics.
  rho_peak  <- max(rho_sympl, na.rm=TRUE)
  if (!is.finite(rho_peak) || rho_peak == 0) rho_peak <- 1
  y_lim_rho <- rho_peak / (GOLDEN_FILL + 0.2)

  # Middle column y_lim — peak of finite WKB amplitude inside any orbit
  # segment, excluding small bands at every turning point so a near-
  # turning-point spike doesn't dominate the scale. We use a global
  # exclusion band rather than per-segment because the y_lim is shared
  # across all segments in the panel.
  inside_pad <- 0.02 * q_span
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

  dt_caustic <- data.table(q=q_display, wkb_density=state$wkb_density)
  dt_rho     <- data.table(q=q_display, rho_sympl=rho_sympl)

  overlay_layers <- symplectic_overlay_layers(cov$Delta_q, cov$Delta_p,
                                              q_center=cov$q_mean)
  df_traj <- classical_trajectory(double_well_V, E_n, tp)

  # Airy overlay for the right column. Same style as the Morse
  # semiclassical figure — line-only (no fill) at gray30/0.35 so the
  # comparator reads the same across both manuscript figures.
  airy_overlay <- list(
    list(
      data       = data.frame(q = q_display, rho = rho_airy),
      color      = "gray30",
      linewidth  = 0.35
    )
  )

  row_label <- sprintf("italic(n)==%d", n)

  list(
    row_label,
    plot_classical_orbit_phase_space(
      df_traj,
      q_lim=c(q_lo,q_hi), p_lim=c(p_lo,p_hi),
      custom_breaks_q=custom_breaks_q,
      custom_breaks_p=custom_breaks_p,
      label_format=label_format, base_font=base_font,
      overlay_layers=overlay_layers,
      orbit_color=HEATMAP_COLOR_HIGH,
      orbit_linewidth=0.3),
    plot_wkb_caustic_cross_section(
      dt_caustic, q_lim=c(q_lo,q_hi), y_lim=y_lim_caustic,
      custom_breaks=custom_breaks_q,
      label_format=label_format,
      q_minus=NULL, q_plus=NULL,  # auto-detected from data
      base_font=base_font),
    plot_semiclassical_resolution(
      dt_rho, q_lim=c(q_lo,q_hi), y_lim=y_lim_rho,
      custom_breaks=custom_breaks_q,
      label_format=label_format, base_font=base_font,
      overlays=airy_overlay)
  )
}

cat("Computing double-well semiclassical-symplectic grid...\n")
rows    <- lapply(target_quantum_numbers,
                  function(n) build_double_well_row(n, base_font=latex_font))
p_final <- assemble_grid(rows,
                         COLUMN_TITLE_CENTER_SEMICLASSICAL,
                         COLUMN_TITLE_RIGHT_SYMPLECTIC,
                         base_font=latex_font)

save_figure(p_final, file_output_pdf, length(target_quantum_numbers))
