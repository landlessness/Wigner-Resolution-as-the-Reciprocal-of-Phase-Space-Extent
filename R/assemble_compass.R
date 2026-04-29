# ==============================================================================
# assemble_compass.R
# Compose the three compass-state plot rows into a single 3-row figure
# for the manuscript: figures/compass.pdf.
#
# Sources the three plot scripts under ASSEMBLY_MODE so they build their
# named panels but do NOT save their per-script standalone PDFs. Then
# adds row labels (Wigner / Husimi / Symplectic) and column titles
# (Phase-Space Cells / Cross Sections / Resolved), places the colormap
# legend in the upper-right empty cell, and writes only compass.pdf.
#
# Layout:
#   Row 1 (Wigner):     bare W heatmap   | W cross-section | colormap legend
#   Row 2 (Husimi):     W + H overlay    | Q cross-section | Q heatmap
#   Row 3 (Symplectic): W + S overlay    | P cross-section | P heatmap
#
# Shared amplitude scale: each plot script computes wigner_max_abs from
# its own Wigner state; all three converge to the same value, so heatmap
# colormaps and cross-section y-axes are consistent across rows. The
# legend tick labels read off this physical scale.
#
# ASSEMBLY_MODE is removed from the global environment after sourcing,
# so subsequent runs of the standalone plot scripts produce their PDFs.
# ==============================================================================

library(here)
library(patchwork)

# Set the gating flag BEFORE sourcing the three plot scripts so they
# build their named panels but skip their per-script ggsave block.
ASSEMBLY_MODE <- TRUE

source(here("R", "plot_compass_wigner.R"))
source(here("R", "plot_compass_wigner_husimi.R"))
source(here("R", "plot_compass_wigner_symplectic.R"))

# Clean up the gating flag so subsequent standalone plot-script runs
# (in the same R session) save their PDFs as expected.
on.exit(rm("ASSEMBLY_MODE", envir=globalenv()), add=TRUE)

dir_figures <- here("figures")
if (!dir.exists(dir_figures)) dir.create(dir_figures, recursive=TRUE)
file_output_pdf <- file.path(dir_figures, "compass.pdf")

# ------------------------------------------------------------------------------
# UPPER-RIGHT EMPTY CELL
# Each cross-section panel now displays its own y-axis tick labels, so the
# numerical scales are visible in-place rather than being summarized in a
# shared legend. The upper-right cell is intentionally left blank.
# ------------------------------------------------------------------------------

upper_right_empty <- patchwork::plot_spacer()

# ------------------------------------------------------------------------------
# ATTACH ROW LABELS to the leftmost panel of each row.
# Plain text, rotated 90 deg (reads bottom-to-top).
# ------------------------------------------------------------------------------

row_top_left <- attach_compass_row_tag(compass_wigner_heatmap,
                                       COMPASS_ROW_LABEL_TOP,
                                       base_font=latex_font)
row_mid_left <- attach_compass_row_tag(compass_wigner_husimi_overlay_heatmap,
                                       COMPASS_ROW_LABEL_MIDDLE,
                                       base_font=latex_font)
row_bot_left <- attach_compass_row_tag(compass_wigner_symplectic_overlay_heatmap,
                                       COMPASS_ROW_LABEL_BOTTOM,
                                       base_font=latex_font)

# ------------------------------------------------------------------------------
# ATTACH COLUMN TITLES to the top row only.
# add_column_titles returns a list of the three top-row panels with titles.
# ------------------------------------------------------------------------------

top_row_titled <- add_column_titles(
  row_top_left, compass_wigner_cross, upper_right_empty,
  COMPASS_COLUMN_TITLE_LEFT,
  COMPASS_COLUMN_TITLE_CENTER,
  COMPASS_COLUMN_TITLE_RIGHT)

# ------------------------------------------------------------------------------
# ASSEMBLE THE 3x3 GRID
# ------------------------------------------------------------------------------

p_full <- (top_row_titled[[1]] | top_row_titled[[2]] | top_row_titled[[3]]) /
          (row_mid_left        | compass_husimi_cross     | compass_husimi_heatmap)   /
          (row_bot_left        | compass_symplectic_cross | compass_symplectic_heatmap)

# ------------------------------------------------------------------------------
# SAVE THE ASSEMBLED PDF
# Three 2.25"-wide columns x three 2.4"-tall rows, plus margin for column
# titles and row labels.
# ------------------------------------------------------------------------------

ggsave(filename=file_output_pdf, plot=p_full, device=cairo_pdf,
       width=7.0, height=7.6, limitsize=FALSE)
cat("Done.", file_output_pdf, "\n")
