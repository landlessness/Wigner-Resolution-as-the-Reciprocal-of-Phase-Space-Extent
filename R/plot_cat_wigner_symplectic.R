# ==============================================================================
# plot_cat_wigner_symplectic.R
# Symplectic resolution of the n-cat Wigner negativity using both conjugate
# quantum blobs G_{a_q} and G_{a_p}, displayed as separate measurements
# plus a joint-min ensemble heatmap.
#
# Four rows:
#   row 1: n=2          (lobes at N, S)
#   row 2: n=3          (apex at N, base lobes at p=-5; triangle pointing up)
#   row 3: n=4 diagonal (lobes at NE, NW, SW, SE — Zurek rotated 45 deg)
#   row 4: n=4 cardinal (lobes at N, E, S, W — Zurek original compass)
#
# Six columns:
#   col 1: bare Wigner heatmap with three-ellipse symplectic overlay
#          (outer Fermi blob A, inner conjugate quantum blobs a_q, a_p)
#   col 2: cross-section P_aq(q, 0) at p=0 (along q)
#   col 3: 2D heatmap of P_aq = W * G_{a_q}
#   col 4: cross-section P_ap(0, p) at q=0 (along p) — complementary to col 2
#   col 5: 2D heatmap of P_ap = W * G_{a_p}
#   col 6: 2D heatmap of P_ens(q, p) = min(P_aq(q,p), P_ap(q,p))
#
# G_{a_q} is the position-squeezed conjugate quantum blob (narrow in q,
# broad in p); G_{a_p} is its complementary momentum-squeezed partner
# (broad in q, narrow in p). Each is the Wigner of a coherent (squeezed)
# state, so by Hudson's theorem each individual P_aq, P_ap is non-negative.
#
# Per the de Gosson / Zurek framing in the manuscript, P_aq and P_ap are
# accessible only through ENSEMBLE TOMOGRAPHY ON SEPARATE COPIES of the
# state — they are two independent measurements on different ensembles,
# not averaged together. Column 6 visualizes their joint structure via
# the pointwise minimum: a phase-space point is "well-resolved by both
# measurements" only if both P_aq and P_ap assign it high density. The
# minimum is a visualization device, NOT a literal joint distribution
# (no such thing exists for complementary observables).
#
# When run standalone, produces figures/cat_wigner_symplectic.pdf at
# 13.5"x9.6" (4 rows x 6 columns at 2.25"x2.4" per panel).
# When sourced under ASSEMBLY_MODE, builds named panels but skips ggsave.
#
# Reference: de Gosson 2009, "Symplectic Methods in Harmonic Analysis";
#            Hudson 1974 Rep. Math. Phys. 6, 249;
#            Zurek 2001 Nature 412, 712 (reciprocal scales delta_q, delta_p).
# ==============================================================================

library(here)
library(patchwork)

source(here("R", "plot_tools.R"))
source(here("R", "cat_system.R"))               # cat_psi, CAT_*
source(here("R", "wigner_density.R"))           # build_wigner_state,
                                                # apply_kernel_density
source(here("R", "symplectic_kernel.R"))        # G_delta_q_kernel_matrix,
                                                # G_delta_p_kernel_matrix,
                                                # symplectic_overlay_layers
source(here("R", "quantum_tools.R"))            # numerical_covariance

latex_font      <- "CMU Serif"
dir_figures     <- here("figures")
if (!dir.exists(dir_figures)) dir.create(dir_figures, recursive=TRUE)
file_output_pdf <- file.path(dir_figures, "cat_wigner_symplectic.pdf")

# ------------------------------------------------------------------------------
# DISPLAY WINDOW (must match all cat_wigner_*.R scripts)
# ------------------------------------------------------------------------------

Q_DISPLAY <- 8

# ------------------------------------------------------------------------------
# DISPLAY GRIDS AND BREAKS (shared across rows)
# ------------------------------------------------------------------------------

q_lo <- -Q_DISPLAY
q_hi <-  Q_DISPLAY
p_lo <- -Q_DISPLAY
p_hi <-  Q_DISPLAY

custom_breaks_q <- c(-CAT_P_MAX, 0, CAT_P_MAX)
custom_breaks_p <- c(-CAT_P_MAX, 0, CAT_P_MAX)
label_format    <- function(x) sprintf("%.0f", x)
q_display       <- seq(q_lo, q_hi, length.out=500)
p_display       <- seq(p_lo, p_hi, length.out=500)
custom_breaks_s <- c(-CAT_P_MAX, 0, CAT_P_MAX)

# ------------------------------------------------------------------------------
# CONFIGURATIONS (same as plot_cat_wigner.R)
# ------------------------------------------------------------------------------

cat_configurations <- list(
  list(label="n2",  desc="n=2 (N, S)",                        n_cats=2, orientation="cardinal"),
  list(label="n3",  desc="n=3 (apex at N, base at p=-5)",     n_cats=3, orientation="cardinal"),
  list(label="n4d", desc="n=4 diagonal (NE, NW, SW, SE)",     n_cats=4, orientation="diagonal"),
  list(label="n4",  desc="n=4 cardinal (N, E, S, W = Zurek)", n_cats=4, orientation="cardinal")
)

# ------------------------------------------------------------------------------
# PER-ROW BUILDER
# Convolves with both G_{a_q} and G_{a_p}, builds cross-sections and 2D
# heatmaps for each, plus the joint-min ensemble heatmap.
# ------------------------------------------------------------------------------

build_cat_symplectic_row <- function(cfg) {
  cat(sprintf("\n=== %s ===\n", cfg$desc))

  cat("Building cat wavefunction...\n")
  q_psi   <- seq(CAT_Q_MIN, CAT_Q_MAX, by=CAT_DQ)
  psi_vec <- cat_psi(q_psi, n_cats=cfg$n_cats, p_max=CAT_P_MAX,
                     xi=CAT_XI, orientation=cfg$orientation)

  rs <- numerical_covariance(psi_vec, q_psi, hbar=1.0)
  cat(sprintf("  p_max=%.1f  xi=%.1f\n", CAT_P_MAX, CAT_XI))
  cat(sprintf("  A_RS/A0=%.2f | Delta_q=%.3f Delta_p=%.3f | delta_q=%.3f delta_p=%.3f\n",
              rs$A_over_A0, rs$Delta_q, rs$Delta_p, rs$delta_q, rs$delta_p))

  cat("\nBuilding Wigner on extended integration grid...\n")
  state <- build_wigner_state(psi_vec, q_psi,
                              q_lo, q_hi, p_lo, p_hi, q_display,
                              n_q_int=1601, n_p_int=1601)

  # ------------------------------------------------------------------------
  # CONVOLVE WITH BOTH CONJUGATE BLOBS (independent measurements).
  # ------------------------------------------------------------------------

  cat("\nApplying G_{a_q} kernel (position-squeezed conjugate blob)...\n")
  gaq_for_state <- function(qg, pg) G_delta_q_kernel_matrix(qg, pg,
                                                            rs$Delta_q, rs$Delta_p)
  P_aq <- apply_kernel_density(state, gaq_for_state,
                               q_lo, q_hi, p_lo, p_hi, q_display)

  cat("\nApplying G_{a_p} kernel (momentum-squeezed conjugate blob)...\n")
  gap_for_state <- function(qg, pg) G_delta_p_kernel_matrix(qg, pg,
                                                            rs$Delta_q, rs$Delta_p)
  P_ap <- apply_kernel_density(state, gap_for_state,
                               q_lo, q_hi, p_lo, p_hi, q_display)

  cat(sprintf("  P_aq cross-section at p=0 range: [%.4e, %.4e]\n",
              min(P_aq$P_cross), max(P_aq$P_cross)))
  cat(sprintf("  P_ap min over heatmap: %.4e\n", min(P_ap$heatmap_dt$w)))

  # ------------------------------------------------------------------------
  # CROSS-SECTIONS
  # P_aq sliced along p=0 (vs q); P_ap sliced along q=0 (vs p).
  # ------------------------------------------------------------------------

  P_ap_cross_q0 <- extract_q0_cross_section(P_ap$P_matrix,
                                            state$q_int, state$p_int, p_display)
  cat(sprintf("  P_ap cross-section at q=0 range: [%.4e, %.4e]\n",
              min(P_ap_cross_q0), max(P_ap_cross_q0)))

  dt_cross_aq <- data.table(q=q_display, P_sympl=P_aq$P_cross)
  dt_cross_ap <- data.table(p=p_display, P_sympl=P_ap_cross_q0)

  y_lim_aq <- max(dt_cross_aq$P_sympl, na.rm=TRUE) * 1.1
  y_lim_ap <- max(dt_cross_ap$P_sympl, na.rm=TRUE) * 1.1

  # ------------------------------------------------------------------------
  # JOINT-MIN ENSEMBLE HEATMAP
  # P_ens(q, p) = min(P_aq(q, p), P_ap(q, p)) -- pointwise minimum of the
  # two complementary measurements. Visualizes phase-space points well-
  # resolved by BOTH ensemble measurements.
  # ------------------------------------------------------------------------

  cat("\nBuilding joint-min ensemble heatmap (pointwise min of P_aq, P_ap)...\n")
  P_ens_matrix <- pmin(P_aq$P_matrix, P_ap$P_matrix)

  q_mask <- state$q_int >= q_lo & state$q_int <= q_hi
  p_mask <- state$p_int >= p_lo & state$p_int <= p_hi
  ens_heatmap_dt <- as.data.table(expand.grid(
    q = state$q_int[q_mask],
    p = state$p_int[p_mask]
  ))
  ens_heatmap_dt[, w := as.vector(P_ens_matrix[q_mask, p_mask])]
  max_w_ens <- max(ens_heatmap_dt$w, na.rm=TRUE)
  ens_heatmap_dt[, w_plot := if (max_w_ens > 0) w/max_w_ens else w]
  cat(sprintf("  P_ens range: [%.4e, %.4e]\n",
              min(ens_heatmap_dt$w), max(ens_heatmap_dt$w)))

  # ------------------------------------------------------------------------
  # PANEL CONSTRUCTION
  # ------------------------------------------------------------------------

  overlay_layers <- symplectic_overlay_layers(rs$Delta_q, rs$Delta_p, q_center=0)
  empty_overlay  <- list()

  # col 1: bare Wigner with three-ellipse overlay
  panel_overlay_heatmap <- plot_wigner_heatmap(
    state$heatmap_dt, overlay_layers,
    q_lim=c(q_lo, q_hi), p_lim=c(p_lo, p_hi),
    custom_breaks_q=custom_breaks_q,
    custom_breaks_p=custom_breaks_p,
    label_format=label_format, base_font=latex_font)

  # col 2: P_aq cross-section at p=0
  panel_aq_cross <- plot_symplectic_cross_section(
    dt_cross_aq,
    q_lim=c(q_lo, q_hi), y_lim=y_lim_aq,
    custom_breaks=custom_breaks_s,
    label_format=label_format,
    base_font=latex_font)

  # col 3: P_aq 2D heatmap
  panel_aq_heatmap <- plot_wigner_heatmap(
    P_aq$heatmap_dt, empty_overlay,
    q_lim=c(q_lo, q_hi), p_lim=c(p_lo, p_hi),
    custom_breaks_q=custom_breaks_q,
    custom_breaks_p=custom_breaks_p,
    label_format=label_format, base_font=latex_font)

  # col 4: P_ap cross-section at q=0
  panel_ap_cross <- plot_symplectic_cross_section_p(
    dt_cross_ap,
    p_lim=c(p_lo, p_hi), y_lim=y_lim_ap,
    custom_breaks=custom_breaks_s,
    label_format=label_format,
    base_font=latex_font)

  # col 5: P_ap 2D heatmap
  panel_ap_heatmap <- plot_wigner_heatmap(
    P_ap$heatmap_dt, empty_overlay,
    q_lim=c(q_lo, q_hi), p_lim=c(p_lo, p_hi),
    custom_breaks_q=custom_breaks_q,
    custom_breaks_p=custom_breaks_p,
    label_format=label_format, base_font=latex_font)

  # col 6: P_ens joint-min 2D heatmap
  panel_ens_heatmap <- plot_wigner_heatmap(
    ens_heatmap_dt, empty_overlay,
    q_lim=c(q_lo, q_hi), p_lim=c(p_lo, p_hi),
    custom_breaks_q=custom_breaks_q,
    custom_breaks_p=custom_breaks_p,
    label_format=label_format, base_font=latex_font)

  list(
    overlay_heatmap = panel_overlay_heatmap,
    aq_cross        = panel_aq_cross,
    aq_heatmap      = panel_aq_heatmap,
    ap_cross        = panel_ap_cross,
    ap_heatmap      = panel_ap_heatmap,
    ens_heatmap     = panel_ens_heatmap
  )
}

# ------------------------------------------------------------------------------
# BUILD ALL FOUR ROWS
# ------------------------------------------------------------------------------

rows <- lapply(cat_configurations, build_cat_symplectic_row)
names(rows) <- vapply(cat_configurations, function(cfg) cfg$label, character(1))

# Named panels for assemble_cat.R.
cat_wigner_symplectic_n2_overlay_heatmap  <- rows$n2$overlay_heatmap
cat_symplectic_n2_aq_cross                <- rows$n2$aq_cross
cat_symplectic_n2_aq_heatmap              <- rows$n2$aq_heatmap
cat_symplectic_n2_ap_cross                <- rows$n2$ap_cross
cat_symplectic_n2_ap_heatmap              <- rows$n2$ap_heatmap
cat_symplectic_n2_ens_heatmap             <- rows$n2$ens_heatmap

cat_wigner_symplectic_n3_overlay_heatmap  <- rows$n3$overlay_heatmap
cat_symplectic_n3_aq_cross                <- rows$n3$aq_cross
cat_symplectic_n3_aq_heatmap              <- rows$n3$aq_heatmap
cat_symplectic_n3_ap_cross                <- rows$n3$ap_cross
cat_symplectic_n3_ap_heatmap              <- rows$n3$ap_heatmap
cat_symplectic_n3_ens_heatmap             <- rows$n3$ens_heatmap

cat_wigner_symplectic_n4d_overlay_heatmap <- rows$n4d$overlay_heatmap
cat_symplectic_n4d_aq_cross               <- rows$n4d$aq_cross
cat_symplectic_n4d_aq_heatmap             <- rows$n4d$aq_heatmap
cat_symplectic_n4d_ap_cross               <- rows$n4d$ap_cross
cat_symplectic_n4d_ap_heatmap             <- rows$n4d$ap_heatmap
cat_symplectic_n4d_ens_heatmap            <- rows$n4d$ens_heatmap

cat_wigner_symplectic_n4_overlay_heatmap  <- rows$n4$overlay_heatmap
cat_symplectic_n4_aq_cross                <- rows$n4$aq_cross
cat_symplectic_n4_aq_heatmap              <- rows$n4$aq_heatmap
cat_symplectic_n4_ap_cross                <- rows$n4$ap_cross
cat_symplectic_n4_ap_heatmap              <- rows$n4$ap_heatmap
cat_symplectic_n4_ens_heatmap             <- rows$n4$ens_heatmap

# ------------------------------------------------------------------------------
# SAVE STANDALONE PDF (skipped when sourced by assemble_cat.R)
# Four rows of six panels each at 2.25"x2.4" per panel -> 13.5"x9.6" total.
# ------------------------------------------------------------------------------

if (!exists("ASSEMBLY_MODE")) {
  p_full <-
    (cat_wigner_symplectic_n2_overlay_heatmap |
     cat_symplectic_n2_aq_cross |
     cat_symplectic_n2_aq_heatmap |
     cat_symplectic_n2_ap_cross |
     cat_symplectic_n2_ap_heatmap |
     cat_symplectic_n2_ens_heatmap) /
    (cat_wigner_symplectic_n3_overlay_heatmap |
     cat_symplectic_n3_aq_cross |
     cat_symplectic_n3_aq_heatmap |
     cat_symplectic_n3_ap_cross |
     cat_symplectic_n3_ap_heatmap |
     cat_symplectic_n3_ens_heatmap) /
    (cat_wigner_symplectic_n4d_overlay_heatmap |
     cat_symplectic_n4d_aq_cross |
     cat_symplectic_n4d_aq_heatmap |
     cat_symplectic_n4d_ap_cross |
     cat_symplectic_n4d_ap_heatmap |
     cat_symplectic_n4d_ens_heatmap) /
    (cat_wigner_symplectic_n4_overlay_heatmap |
     cat_symplectic_n4_aq_cross |
     cat_symplectic_n4_aq_heatmap |
     cat_symplectic_n4_ap_cross |
     cat_symplectic_n4_ap_heatmap |
     cat_symplectic_n4_ens_heatmap)
  ggsave(filename=file_output_pdf, plot=p_full, device=cairo_pdf,
         width=13.5, height=9.6, limitsize=FALSE)
  cat("\nDone.", file_output_pdf, "\n")
}
