# ==============================================================================
# plot_cat_wigner.R
# Cat-state Wigner ground truth for four configurations:
#
#   row 1: n=2          (lobes at N, S)
#   row 2: n=3          (apex at N, base lobes at p=-5; triangle pointing up)
#   row 3: n=4 diagonal (lobes at NE, NW, SW, SE — Zurek rotated 45 deg)
#   row 4: n=4 cardinal (lobes at N, E, S, W — Zurek original compass)
#
# Each row has two columns:
#   Left:    bare Wigner heatmap with display window +/-Q_DISPLAY
#   Right:   horizontal cross-section W(q, p=0) along q
#
# All rows are anchored so that the extreme p coordinates of the lobes
# are at +/- 5, guaranteeing the p=0 horizontal cross-section sits halfway
# between extreme-p lobes in every row. The cross-section therefore cuts
# through canonical interference structure with visible Wigner negativity
# for each configuration.
#
# When run standalone, produces figures/cat_wigner.pdf at 4.5"x9.6".
# When sourced under ASSEMBLY_MODE, builds named panels but skips ggsave.
#
# Reference: Wikipedia "Cat state" for the canonical n=2,3,4 visualization;
#            Zurek 2001 Nature 412, 712 for the n=4 compass case.
# ==============================================================================

library(here)
library(patchwork)

source(here("R", "plot_tools.R"))
source(here("R", "cat_system.R"))               # cat_psi, cat_lobe_positions, CAT_*
source(here("R", "wigner_density.R"))           # build_wigner_state

latex_font      <- "CMU Serif"
dir_figures     <- here("figures")
if (!dir.exists(dir_figures)) dir.create(dir_figures, recursive=TRUE)
file_output_pdf <- file.path(dir_figures, "cat_wigner.pdf")

# ------------------------------------------------------------------------------
# DISPLAY WINDOW
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
custom_breaks_s <- c(-CAT_P_MAX, 0, CAT_P_MAX)

# ------------------------------------------------------------------------------
# CONFIGURATIONS
# Each entry: short label (used for named-panel suffixes), human-readable
# description, n_cats, orientation. Row order matches the figure order.
# ------------------------------------------------------------------------------

cat_configurations <- list(
  list(label="n2",  desc="n=2 (N, S)",                       n_cats=2, orientation="cardinal"),
  list(label="n3",  desc="n=3 (apex at N, base at p=-5)",    n_cats=3, orientation="cardinal"),
  list(label="n4d", desc="n=4 diagonal (NE, NW, SW, SE)",    n_cats=4, orientation="diagonal"),
  list(label="n4",  desc="n=4 cardinal (N, E, S, W = Zurek)", n_cats=4, orientation="cardinal")
)

# ------------------------------------------------------------------------------
# PER-ROW BUILDER
# ------------------------------------------------------------------------------

build_cat_row <- function(cfg) {
  cat(sprintf("\n=== %s ===\n", cfg$desc))

  cat("Building cat wavefunction...\n")
  q_psi   <- seq(CAT_Q_MIN, CAT_Q_MAX, by=CAT_DQ)
  psi_vec <- cat_psi(q_psi, n_cats=cfg$n_cats, p_max=CAT_P_MAX,
                     xi=CAT_XI, orientation=cfg$orientation)
  lobes   <- cat_lobe_positions(cfg$n_cats, CAT_P_MAX, cfg$orientation)
  cat(sprintf("  p_max=%.1f  xi=%.1f  nq=%d  dq=%.4f\n",
              CAT_P_MAX, CAT_XI, length(q_psi), CAT_DQ))
  cat(sprintf("  Lobe positions:\n"))
  for (k in seq_along(lobes$q)) {
    cat(sprintf("    %d: (%.3f, %.3f)\n", k, lobes$q[k], lobes$p[k]))
  }

  cat("\nBuilding Wigner on extended integration grid...\n")
  state <- build_wigner_state(psi_vec, q_psi,
                              q_lo, q_hi, p_lo, p_hi, q_display,
                              n_q_int=1601, n_p_int=1601)

  cat("\nExtracting horizontal cross-section at p=0...\n")
  W_cross <- extract_p0_cross_section(state$W_matrix,
                                      state$q_int, state$p_int, q_display)
  dt_cross <- data.table(q=q_display, W_raw=W_cross)
  cat(sprintf("  W cross-section range: [%.4f, %.4f]\n",
              min(W_cross), max(W_cross)))

  w_max <- max(abs(dt_cross$W_raw), na.rm=TRUE)
  y_lim <- w_max * 1.1

  empty_overlay <- list()

  heatmap_panel <- plot_wigner_heatmap(
    state$heatmap_dt, empty_overlay,
    q_lim=c(q_lo, q_hi), p_lim=c(p_lo, p_hi),
    custom_breaks_q=custom_breaks_q,
    custom_breaks_p=custom_breaks_p,
    label_format=label_format, base_font=latex_font)

  cross_panel <- plot_wigner_cross_section(
    dt_cross,
    q_lim=c(q_lo, q_hi), y_lim=y_lim,
    custom_breaks=custom_breaks_s,
    label_format=label_format,
    base_font=latex_font)

  list(heatmap=heatmap_panel, cross=cross_panel)
}

# ------------------------------------------------------------------------------
# BUILD ALL FOUR ROWS
# ------------------------------------------------------------------------------

rows <- lapply(cat_configurations, build_cat_row)
names(rows) <- vapply(cat_configurations, function(cfg) cfg$label, character(1))

# Named panels available to assemble_cat.R (when sourced under ASSEMBLY_MODE).
cat_wigner_n2_heatmap  <- rows$n2$heatmap
cat_wigner_n2_cross    <- rows$n2$cross
cat_wigner_n3_heatmap  <- rows$n3$heatmap
cat_wigner_n3_cross    <- rows$n3$cross
cat_wigner_n4d_heatmap <- rows$n4d$heatmap
cat_wigner_n4d_cross   <- rows$n4d$cross
cat_wigner_n4_heatmap  <- rows$n4$heatmap
cat_wigner_n4_cross    <- rows$n4$cross

# ------------------------------------------------------------------------------
# SAVE STANDALONE PDF (skipped under ASSEMBLY_MODE)
# Four rows of two panels each at 4.5"x2.4" per row -> 4.5"x9.6" total.
# ------------------------------------------------------------------------------

if (!exists("ASSEMBLY_MODE")) {
  p_full <-
    (cat_wigner_n2_heatmap  | cat_wigner_n2_cross)  /
    (cat_wigner_n3_heatmap  | cat_wigner_n3_cross)  /
    (cat_wigner_n4d_heatmap | cat_wigner_n4d_cross) /
    (cat_wigner_n4_heatmap  | cat_wigner_n4_cross)
  ggsave(filename=file_output_pdf, plot=p_full, device=cairo_pdf,
         width=4.5, height=9.6, limitsize=FALSE)
  cat("\nDone.", file_output_pdf, "\n")
}
