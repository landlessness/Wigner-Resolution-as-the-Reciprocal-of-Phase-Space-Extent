# ==============================================================================
# plot_cats.R
# Symplectic resolution of Wigner negativity for n-cat states.
# Each row shows one n-cat configuration:
#   row 1: n_cats=2          (lobes at (0, +/-5))
#   row 2: n_cats=3          (equilateral triangle pointing up)
#   row 3: n_cats=4 diag     (compass with lobes on the diagonals)
#   row 4: n_cats=4 axis     (diag rotated 45 deg)
#
# Three columns:
#   left:    bare Wigner heatmap with symplectic overlay.
#   center:  raw Wigner cross-section W(q, 0).
#   right:   symplectic resolution P_{delta q}(q, 0).
#
# State-building (psi, RS covariance, Wigner via FFT, symplectic overlay)
# is delegated to build_cat_state() in state_builder.R.
#
# Author: Brian S. Mulloy
# ==============================================================================

library(here)
library(data.table)
library(patchwork)

source(here("R", "plot_tools.R"))

source(here("R", "cat_system.R"))                # cat_psi, CAT_*
source(here("R", "wigner_density.R"))            # apply_kernel_cross_section
source(here("R", "symplectic_kernel.R"))         # G_delta_q_kernel_matrix
source(here("R", "husimi_kernel.R"))             # husimi_kernel_matrix
source(here("R", "state_builder.R"))             # build_cat_state

latex_font  <- "CMU Serif"
dir_figures <- here("figures")
if (!dir.exists(dir_figures)) dir.create(dir_figures, recursive=TRUE)
file_output_pdf <- file.path(dir_figures, "cats.pdf")

# ------------------------------------------------------------------------------
# OVERLAY CELLS
#
# Which symplectic cells to draw on the heatmap. Subset of:
#   "A"   - extended (Fermi blob)
#   "a_q" - q-squeezed (tall ellipse: small in q, full extent in p)
#   "a_p" - p-squeezed (wide ellipse: full extent in q, small in p)
# ------------------------------------------------------------------------------
OVERLAY_CELLS <- c("A", "a_q")

# ------------------------------------------------------------------------------
# DISPLAY WINDOW (shared across all four rows)
# ------------------------------------------------------------------------------

q_lo <- -CAT_Q_DISPLAY
q_hi <-  CAT_Q_DISPLAY
p_lo <- -CAT_Q_DISPLAY
p_hi <-  CAT_Q_DISPLAY

custom_breaks_q <- c(-CAT_P_MAX, 0, CAT_P_MAX)
custom_breaks_p <- c(-CAT_P_MAX, 0, CAT_P_MAX)
label_format    <- function(x) sprintf("%.0f", x)

# ------------------------------------------------------------------------------
# ROW DESCRIPTORS
#
# Each descriptor names an n-cat configuration plus the shared display
# window so build_cat_state() can produce a uniform state bundle.
# ------------------------------------------------------------------------------

make_cat_descriptor <- function(n_cats, variant) {
  list(
    n_cats          = n_cats,
    variant         = variant,
    hbar            = CAT_HBAR,
    q_lo            = q_lo, q_hi = q_hi,
    p_lo            = p_lo, p_hi = p_hi,
    custom_breaks_q = custom_breaks_q,
    custom_breaks_p = custom_breaks_p,
    label_format    = label_format
  )
}

row_descriptors <- list(
  make_cat_descriptor(2, "diag"),
  make_cat_descriptor(3, "diag"),
  make_cat_descriptor(4, "diag"),
  make_cat_descriptor(4, "axis")
)

# ------------------------------------------------------------------------------
# ROW BUILDER
#
# Calls build_cat_state() for the shared per-state pipeline, then adds
# the symplectic cross-section and Husimi cross-section that this figure
# needs (and the tomography pipeline does not).
# ------------------------------------------------------------------------------

build_cat_row <- function(descriptor, base_font="",
                          cells=c("A", "a_q", "a_p")) {
  ps <- build_cat_state(descriptor, base_font=base_font, cells=cells)

  # Symplectic kernel cross-section.
  symplectic_kernel_for_state <- function(qg, pg) {
    G_delta_q_kernel_matrix(qg, pg, ps$rs$Delta_q, ps$rs$Delta_p,
                            hbar=ps$hbar)
  }
  cat("  Convolving with symplectic kernel...\n")
  P_sympl_cross <- apply_kernel_cross_section(ps$state,
                                              symplectic_kernel_for_state,
                                              ps$q_display)

  # Husimi cross-section.
  husimi_kernel_for_state <- function(qg, pg) {
    husimi_kernel_matrix(qg, pg)
  }
  cat("  Computing Husimi cross-section at p=0...\n")
  Q_husimi_cross <- apply_kernel_cross_section(ps$state,
                                               husimi_kernel_for_state,
                                               ps$q_display)

  # Y-scaling.
  W_cross_peak <- max(abs(ps$state$W_cross), na.rm=TRUE)
  if (!is.finite(W_cross_peak) || W_cross_peak == 0) W_cross_peak <- 1
  y_lim_W <- W_cross_peak * 1.1

  P_peak_data <- max(P_sympl_cross, na.rm=TRUE)
  if (!is.finite(P_peak_data) || P_peak_data == 0) P_peak_data <- 1
  y_lim_P <- P_peak_data * 1.1

  dt_W       <- data.table(q=ps$q_display, W_raw=ps$state$W_cross)
  dt_P_sympl <- data.table(q=ps$q_display, rho_sympl=P_sympl_cross)

  list(
    plot_eigen_heatmap(
      ps$state$heatmap_dt, ps$overlay_layers, df_traj=NULL,
      q_lim=c(ps$q_lo, ps$q_hi), p_lim=c(ps$p_lo, ps$p_hi),
      custom_breaks_q=ps$custom_breaks_q,
      custom_breaks_p=ps$custom_breaks_p,
      label_format=ps$label_format, base_font=base_font),
    plot_eigen_cross_section(
      dt_W, q_lim=c(ps$q_lo, ps$q_hi), y_lim=y_lim_W,
      custom_breaks=ps$custom_breaks_q,
      label_format=ps$label_format, base_font=base_font),
    plot_semiclassical_resolution(
      dt_P_sympl, q_lim=c(ps$q_lo, ps$q_hi), y_lim=y_lim_P,
      custom_breaks=ps$custom_breaks_q,
      label_format=ps$label_format, base_font=base_font,
      overlays=NULL,
      y_label=expression(italic(P)[italic(delta*q)](italic(q)*","*0)))
  )
}

# ------------------------------------------------------------------------------
# DRIVE
# ------------------------------------------------------------------------------

cat("Computing cats figure (4 cat configurations x 3 panels)...\n")

rows <- lapply(row_descriptors,
               function(d) build_cat_row(d, base_font=latex_font,
                                         cells=OVERLAY_CELLS))

p_final <- assemble_grid_unlabeled(rows,
                                   COLUMN_TITLE_CENTER_WIGNER,
                                   COLUMN_TITLE_RIGHT_SYMPLECTIC,
                                   base_font=latex_font)

save_figure(p_final, file_output_pdf, length(row_descriptors))
