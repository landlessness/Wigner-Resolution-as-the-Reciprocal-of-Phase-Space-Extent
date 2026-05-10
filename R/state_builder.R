# ==============================================================================
# state_builder.R
# Shared per-state pipeline: from descriptor to fully-built phase-space state.
#
# This file extracts the state-building logic that was duplicated between
# plot_wigner.R (build_wigner_row) and plot_cats.R (build_cat_row), so that
# every consumer -- the existing plot files, the new tomography pipeline,
# and the marginal-residual scripts -- starts from the same data.
#
# A state bundle is a single list containing everything downstream code
# needs:
#   $name         short identifier
#   $psi_vec, $psi_q_grid    wavefunction and its sampling grid
#   $rs           Robertson-Schroedinger covariance:
#                   Delta_q, Delta_p, delta_q, delta_p, q_mean, A_over_A0
#   $q_lo, $q_hi, $p_lo, $p_hi     display windows
#   $q_display    display q-grid (length 500 by default)
#   $custom_breaks_q, $custom_breaks_p    axis breaks
#   $label_format function
#   $hbar         Planck constant in chosen units
#   $state        output of build_wigner_state(): q_int, p_int, dq_int,
#                 dp_int, W_matrix, W_cross, heatmap_dt, norm
#   $overlay_layers     symplectic overlay layers (subset of Fermi blob A
#                       and conjugate quantum blobs a_q, a_p) ready to add
#                       to a heatmap plot. Subset selected by the `cells`
#                       argument to the build function (default: all
#                       three).
#
# Two descriptor flavors are supported:
#
#   Eigenstate descriptor (matches plot_wigner.R style):
#     $name, $V, $n_target, $E_fn(n), $psi_fn(n, q),
#     $q_window(E), $p_window(E),
#     $q_breaks_fn(E), $p_breaks_fn(E),
#     $psi_q_grid
#
#   Cat descriptor (matches plot_cats.R style):
#     $name, $n_cats, $variant, $hbar (defaults to CAT_HBAR)
#     $q_window, $p_window, $custom_breaks_q, $custom_breaks_p
#       (these are static for cats; the cat figure shares one window
#        across all four rows so they're passed in directly rather than
#        as functions of E)
#
# Use build_eigenstate_state() for harmonic / Morse / double-well /
# squeezed-vacuum and build_cat_state() for cat configurations. Both
# produce the same bundle structure.
#
# Author: Brian S. Mulloy
# ==============================================================================

library(here)
library(data.table)

source(here("R", "quantum_tools.R"))         # numerical_covariance
source(here("R", "wigner_density.R"))         # build_wigner_state
source(here("R", "symplectic_kernel.R"))      # symplectic_overlay_layers
source(here("R", "cat_system.R"))             # cat_psi (only used for cat
                                              # path; harmless to source
                                              # in either path)

# ------------------------------------------------------------------------------
# CACHE SCHEMA
#
# Every entry in data/tomography_data.rds is a bundle returned by
# c(state_builder_output, run_tomography_output). Bumping this version
# any time the bundle's expected fields change forces stale caches to
# be rebuilt rather than producing cryptic ggplot errors downstream.
# ------------------------------------------------------------------------------

TOMOGRAPHY_CACHE_VERSION <- 2

# Field list used to validate a bundle. Bumping the version is the
# preferred mechanism but this list provides a structural fallback when
# (e.g.) two pipeline branches share a version number but differ in
# field set.
TOMOGRAPHY_REQUIRED_FIELDS <- c(
  "name", "rs", "psi_vec", "psi_q_grid",
  "q_lo", "q_hi", "p_lo", "p_hi",
  "custom_breaks_q", "custom_breaks_p", "label_format",
  "state", "overlay_layers",
  "q_grid_tomo", "p_grid_tomo", "x_grid", "theta_grid",
  "W_tomo", "sinogram", "tilde_W",
  "rho_q", "rho_p", "exact_q", "exact_p"
)

#' Load the tomography cache, rebuilding it if stale or missing.
#'
#' @param cache_path  Path to the .rds file.
#' @param compute_path Path to the tomography_compute.R script.
#' @return A named list of state bundles.
load_tomography_cache <- function(cache_path, compute_path) {
  rebuild <- FALSE
  reason  <- ""

  if (!file.exists(cache_path)) {
    rebuild <- TRUE
    reason  <- "cache file missing"
  } else {
    blob <- readRDS(cache_path)
    if (is.null(attr(blob, "version")) ||
        attr(blob, "version") != TOMOGRAPHY_CACHE_VERSION) {
      rebuild <- TRUE
      reason  <- sprintf("cache version mismatch (have %s, need %d)",
                         as.character(attr(blob, "version")),
                         TOMOGRAPHY_CACHE_VERSION)
    } else if (length(blob) < 1) {
      rebuild <- TRUE
      reason  <- "cache empty"
    } else {
      sample_bundle <- blob[[1]]
      missing       <- setdiff(TOMOGRAPHY_REQUIRED_FIELDS, names(sample_bundle))
      if (length(missing) > 0) {
        rebuild <- TRUE
        reason  <- sprintf("cache missing fields: %s",
                           paste(missing, collapse=", "))
      }
    }
  }

  if (rebuild) {
    cat(sprintf("Cache rebuild required: %s\n", reason))
    cat(sprintf("Sourcing %s ...\n", compute_path))
    source(compute_path)
    blob <- readRDS(cache_path)
  }
  blob
}

# ------------------------------------------------------------------------------
# EIGENSTATE PATH
#
# Mirrors steps 1-9 of build_wigner_row() in plot_wigner.R. Output is a
# state bundle ready for either the existing 3-panel plot or the new
# tomography pipeline.
#
# The `cells` argument controls which symplectic overlay cells are
# rendered on the heatmap. Defaults to all three; pass e.g.
# c("A", "a_q") to show only the extended cell and the q-squeezed cell.
# ------------------------------------------------------------------------------

build_eigenstate_state <- function(descriptor, base_font="",
                                   cells = c("A", "a_q", "a_p")) {
  n   <- descriptor$n_target
  E_n <- descriptor$E_fn(n)

  cat(sprintf("\n== %s | n=%d | E_n=%.4f ==\n",
              descriptor$name, n, E_n))

  # 1. Sample psi.
  psi_q_grid <- descriptor$psi_q_grid
  psi_vec    <- descriptor$psi_fn(n, psi_q_grid)

  # 2. RS covariance.
  rs <- numerical_covariance(psi_vec, psi_q_grid, hbar=1.0)
  cat(sprintf("  RS: A_RS/A0=%.4f Delta_q=%.3f Delta_p=%.3f delta_q=%.3f delta_p=%.3f\n",
              rs$A_over_A0, rs$Delta_q, rs$Delta_p, rs$delta_q, rs$delta_p))

  # 3. Display windows and axis breaks.
  qw <- descriptor$q_window(E_n)
  pw <- descriptor$p_window(E_n)
  q_lo <- qw$q_lo; q_hi <- qw$q_hi
  p_lo <- pw$p_lo; p_hi <- pw$p_hi

  custom_breaks_q <- descriptor$q_breaks_fn(E_n)
  custom_breaks_p <- descriptor$p_breaks_fn(E_n)
  label_format    <- function(x) sprintf("%.1f", x)
  q_display       <- seq(q_lo, q_hi, length.out=500)

  # 4. Wigner state.
  cat("  Building Wigner...\n")
  state <- build_wigner_state(
    psi_vec=psi_vec, psi_q_grid=psi_q_grid,
    q_lo=q_lo, q_hi=q_hi, p_lo=p_lo, p_hi=p_hi,
    q_display=q_display)

  # 5. Symplectic overlay layers, centered at the visual peak of W(q,0).
  q_center_visual <- q_display[which.max(abs(state$W_cross))]
  overlay_layers  <- symplectic_overlay_layers(rs$Delta_q, rs$Delta_p,
                                                q_center=q_center_visual,
                                                cells=cells)

  list(
    name            = descriptor$name,
    psi_vec         = psi_vec,
    psi_q_grid      = psi_q_grid,
    rs              = rs,
    q_lo            = q_lo, q_hi = q_hi,
    p_lo            = p_lo, p_hi = p_hi,
    q_display       = q_display,
    custom_breaks_q = custom_breaks_q,
    custom_breaks_p = custom_breaks_p,
    label_format    = label_format,
    hbar            = 1.0,
    state           = state,
    overlay_layers  = overlay_layers,
    n               = n,
    E_n             = E_n
  )
}

# ------------------------------------------------------------------------------
# CAT PATH
#
# Mirrors steps 1-8 of build_cat_row() in plot_cats.R. The cat figure
# uses a single shared display window across all four rows; we accept
# that window as descriptor fields rather than re-computing per-row.
#
# The `cells` argument controls which symplectic overlay cells are
# rendered on the heatmap. Defaults to all three.
# ------------------------------------------------------------------------------

build_cat_state <- function(descriptor, base_font="",
                            cells = c("A", "a_q", "a_p")) {
  n_cats  <- descriptor$n_cats
  variant <- descriptor$variant
  hbar    <- if (is.null(descriptor$hbar)) CAT_HBAR else descriptor$hbar

  cat(sprintf("\n== n_cats=%d variant=%s ==\n", n_cats, variant))

  # 1. Sample psi on the cat wavefunction grid.
  q_psi   <- seq(CAT_Q_MIN, CAT_Q_MAX, by=CAT_DQ)
  psi_vec <- cat_psi(q_psi, n_cats, variant=variant,
                     p_max=CAT_P_MAX, xi=CAT_XI, hbar=hbar)

  # 2. RS covariance.
  rs <- numerical_covariance(psi_vec, q_psi, hbar=hbar)
  cat(sprintf("  RS: A_RS/A0=%.4f Delta_q=%.3f Delta_p=%.3f delta_q=%.3f delta_p=%.3f q_mean=%.3f\n",
              rs$A_over_A0, rs$Delta_q, rs$Delta_p,
              rs$delta_q, rs$delta_p, rs$q_mean))

  # 3. Display window comes from the descriptor (cats share one window).
  q_lo <- descriptor$q_lo; q_hi <- descriptor$q_hi
  p_lo <- descriptor$p_lo; p_hi <- descriptor$p_hi
  custom_breaks_q <- descriptor$custom_breaks_q
  custom_breaks_p <- descriptor$custom_breaks_p
  label_format    <- if (is.null(descriptor$label_format))
                       function(x) sprintf("%.0f", x)
                     else descriptor$label_format
  q_display       <- seq(q_lo, q_hi, length.out=500)

  # 4. Wigner state.
  cat("  Building Wigner...\n")
  state <- build_wigner_state(
    psi_vec=psi_vec, psi_q_grid=q_psi,
    q_lo=q_lo, q_hi=q_hi, p_lo=p_lo, p_hi=p_hi,
    q_display=q_display)

  # 5. Symplectic overlay layers.
  overlay_layers <- symplectic_overlay_layers(rs$Delta_q, rs$Delta_p,
                                               q_center=rs$q_mean,
                                               hbar=hbar,
                                               cells=cells)

  list(
    name            = if (is.null(descriptor$name))
                         sprintf("cat_%d_%s", n_cats, variant)
                       else descriptor$name,
    psi_vec         = psi_vec,
    psi_q_grid      = q_psi,
    rs              = rs,
    q_lo            = q_lo, q_hi = q_hi,
    p_lo            = p_lo, p_hi = p_hi,
    q_display       = q_display,
    custom_breaks_q = custom_breaks_q,
    custom_breaks_p = custom_breaks_p,
    label_format    = label_format,
    hbar            = hbar,
    state           = state,
    overlay_layers  = overlay_layers,
    n_cats          = n_cats,
    variant         = variant
  )
}
