# ==============================================================================
# plot_tools.R
# Display functions only — no computation, no physics.
#
# Layout architecture:
#   Three columns: heatmap, Wigner cross-section, Husimi cross-section.
#   Row labels are attached as patchwork plot tags on each row's heatmap
#   panel. tag_position is set in panel npc coordinates, so the label
#   tracks the panel regardless of how coord_fixed() squeezes the panel's
#   render width.
#
#   The heatmap panel reserves extra left margin (HEATMAP_LEFT_MARGIN_PT)
#   to give the tag horizontal room. Without this margin, the tag's
#   negative-x npc position renders outside the panel's drawn area and
#   gets clipped by cairo_pdf at the figure's left edge.
#
# Heatmap colormap:
#   Asymmetric diverging — positive W renders darkest, negative W renders
#   medium gray, zero is white. The heatmap itself shows where W is
#   negative.
#
# Author: Brian S. Mulloy
# ==============================================================================

library(data.table)
library(ggplot2)
library(ggforce)
library(patchwork)

# ------------------------------------------------------------------------------
# FIGURE-LEVEL CONSTANTS
# ------------------------------------------------------------------------------

COLUMN_TITLE_LEFT   <- "Phase-Space Cells"
COLUMN_TITLE_CENTER <- "Wigner Negativity"
COLUMN_TITLE_RIGHT  <- "Husimi Resolution"

# Three-column layout (heatmap, Wigner, Husimi). Row labels attach as tags.
PANEL_WIDTHS    <- c(1, 1, 1)
FIGURE_WIDTH_IN <- 7.5
ROW_HEIGHT_IN   <- 1.8
FIGURE_PAD_IN   <- 0.5

# Left margin reserved on the heatmap panel for the row tag. The tag sits
# at npc x = TAG_X_NPC (negative => left of panel's left edge); the panel's
# left margin must be wide enough that this position lies within the panel's
# rendered area, otherwise cairo_pdf clips it at the figure boundary.
# 80 pt accommodates labels up to "n = 100" with a small safety margin.
HEATMAP_LEFT_MARGIN_PT <- 80

# Heatmap colormap: asymmetric diverging.
HEATMAP_COLOR_NEG  <- "gray45"   # negative W: medium gray
HEATMAP_COLOR_ZERO <- "white"
HEATMAP_COLOR_POS  <- "gray10"   # positive W: darkest

# Husimi-circle linetype: short dashes so small circles don't look like an
# orphaned dash.
HUSIMI_LINETYPE <- "22"

# Row label tag.
ROW_LABEL_SIZE_PT <- 11
TAG_X_NPC <- -0.18   # left of panel's left edge
TAG_Y_NPC <-  0.5    # vertical center of panel

# ------------------------------------------------------------------------------
# PHASE SPACE HEATMAP — Husimi variant
# ------------------------------------------------------------------------------

plot_phase_space_heatmap_husimi <- function(dt_w2d, husimi_ell, df_traj=NULL,
                                            q_lim, p_lim,
                                            custom_breaks_q, custom_breaks_p,
                                            label_format, base_font="") {
  ax_x <- expression(italic(q)/italic(q)[0])
  ax_y <- expression(italic(p)/italic(p)[0])

  p <- ggplot(dt_w2d, aes(x=q, y=p)) +
    geom_raster(aes(fill=w_plot), interpolate=TRUE) +
    scale_fill_gradient2(low=HEATMAP_COLOR_NEG,
                         mid=HEATMAP_COLOR_ZERO,
                         high=HEATMAP_COLOR_POS,
                         midpoint=0, limits=c(-1,1), guide="none")

  if (!is.null(df_traj))
    p <- p + geom_path(data=df_traj, aes(x=q, y=p), inherit.aes=FALSE,
                       color="black", linewidth=0.5, linetype="solid")

  p +
    geom_circle(data=husimi_ell$circle, aes(x0=x0, y0=y0, r=r),
                inherit.aes=FALSE,
                color="gray20", linewidth=0.4, linetype=HUSIMI_LINETYPE) +
    coord_fixed(xlim=q_lim, ylim=p_lim, expand=FALSE) +
    scale_x_continuous(breaks=custom_breaks_q, labels=label_format) +
    scale_y_continuous(breaks=custom_breaks_p, labels=label_format) +
    theme_bw(base_family=base_font) +
    theme(panel.grid.minor=element_blank(),
          panel.background=element_rect(fill="white"),
          axis.text=element_text(size=8),
          plot.margin=margin(2, 4, 2, 4)) +
    labs(x=ax_x, y=ax_y)
}

# ------------------------------------------------------------------------------
# CROSS-SECTION PANELS
# ------------------------------------------------------------------------------

plot_wigner_cross_section <- function(dt, q_lim, y_lim, custom_breaks, label_format, base_font="") {
  ax_x <- expression(italic(q)/italic(q)[0])
  ax_y <- expression(italic(W)[italic(n)](italic(q)*","*0))
  ggplot(dt, aes(x=q, y=W_raw)) +
    geom_hline(yintercept=0, color="black", linewidth=0.3) +
    geom_ribbon(aes(ymin=pmin(W_raw,0), ymax=0), fill="gray60", alpha=0.6, color=NA) +
    geom_ribbon(aes(ymin=0, ymax=pmax(W_raw,0)), fill="gray85", color=NA) +
    geom_path(color="black", linewidth=0.4) +
    coord_cartesian(xlim=q_lim, ylim=c(-y_lim,y_lim), expand=FALSE) +
    scale_x_continuous(breaks=custom_breaks, labels=label_format) +
    theme_bw(base_family=base_font) +
    theme(panel.grid.minor=element_blank(),
          axis.text=element_text(size=8),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          aspect.ratio=1, plot.margin=margin(2,4,2,4)) +
    labs(x=ax_x, y=ax_y)
}

plot_husimi_cross_section <- function(dt, q_lim, y_lim, custom_breaks, label_format, base_font="") {
  ax_x <- expression(italic(q)/italic(q)[0])
  ax_y <- expression(italic(Q)(italic(q)*","*0))
  ggplot(dt, aes(x=q, y=Q_husimi)) +
    geom_hline(yintercept=0, color="black", linewidth=0.3) +
    geom_ribbon(aes(ymin=0, ymax=pmax(Q_husimi,0)), fill="gray85", color=NA) +
    geom_path(color="black", linewidth=0.4) +
    coord_cartesian(xlim=q_lim, ylim=c(-y_lim,y_lim), expand=FALSE) +
    scale_x_continuous(breaks=custom_breaks, labels=label_format) +
    theme_bw(base_family=base_font) +
    theme(panel.grid.minor=element_blank(),
          axis.text=element_text(size=8),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          aspect.ratio=1, plot.margin=margin(2,4,2,4)) +
    labs(x=ax_x, y=ax_y)
}

# ------------------------------------------------------------------------------
# ROW LABEL via PATCHWORK TAG
# ------------------------------------------------------------------------------

attach_row_tag <- function(p, label_str, base_font="") {
  p +
    labs(tag = parse(text=label_str)) +
    theme(
      plot.tag = element_text(size=ROW_LABEL_SIZE_PT, family=base_font,
                              hjust=1, vjust=0.5),
      plot.tag.position = c(TAG_X_NPC, TAG_Y_NPC),
      plot.margin = margin(2, 4, 2, HEATMAP_LEFT_MARGIN_PT)
    )
}

# ------------------------------------------------------------------------------
# COLUMN TITLES (applied to first row only)
# ------------------------------------------------------------------------------

add_column_titles <- function(p_left, p_center, p_right,
                              title_left, title_center, title_right) {
  list(
    p_left   + labs(title=title_left)  + theme(plot.title=element_text(size=11, hjust=0.5)),
    p_center + labs(title=title_center)+ theme(plot.title=element_text(size=11, hjust=0.5)),
    p_right  + labs(title=title_right) + theme(plot.title=element_text(size=11, hjust=0.5))
  )
}

# ------------------------------------------------------------------------------
# AXIS TITLE MANAGEMENT
# ------------------------------------------------------------------------------

suppress_x_titles <- function(...) {
  lapply(list(...), function(p) p + theme(axis.title.x=element_blank()))
}

suppress_y_titles <- function(...) {
  lapply(list(...), function(p) p + theme(axis.title.y=element_blank()))
}

# ------------------------------------------------------------------------------
# GRID ASSEMBLY
# ------------------------------------------------------------------------------

assemble_wigner_husimi_grid <- function(rows, base_font="") {
  num_rows <- length(rows)
  plot_list <- list()

  for (i in seq_along(rows)) {
    row      <- rows[[i]]
    label    <- row[[1]]
    p_left   <- row[[2]]
    p_center <- row[[3]]
    p_right  <- row[[4]]

    p_left <- attach_row_tag(p_left, label, base_font=base_font)

    if (i == 1) {
      panels <- add_column_titles(p_left, p_center, p_right,
                                  COLUMN_TITLE_LEFT,
                                  COLUMN_TITLE_CENTER,
                                  COLUMN_TITLE_RIGHT)
      p_left <- panels[[1]]; p_center <- panels[[2]]; p_right <- panels[[3]]
    }
    if (i != num_rows) {
      lst <- suppress_x_titles(p_left, p_center, p_right)
      p_left <- lst[[1]]; p_center <- lst[[2]]; p_right <- lst[[3]]
    }
    if (i != ceiling(num_rows/2)) {
      lst <- suppress_y_titles(p_left, p_center, p_right)
      p_left <- lst[[1]]; p_center <- lst[[2]]; p_right <- lst[[3]]
    }

    plot_list <- c(plot_list, list(p_left, p_center, p_right))
  }

  wrap_plots(plot_list, ncol=3, widths=PANEL_WIDTHS) +
    theme(plot.margin=margin(10,10,10,10))
}

# ------------------------------------------------------------------------------
# SAVE FIGURE
# ------------------------------------------------------------------------------

save_figure <- function(p, filepath, n_rows) {
  ggsave(filename=filepath, plot=p, device=cairo_pdf,
         width=FIGURE_WIDTH_IN,
         height=n_rows*ROW_HEIGHT_IN + FIGURE_PAD_IN,
         limitsize=FALSE)
  cat("Done.", filepath, "\n")
}
