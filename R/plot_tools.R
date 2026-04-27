# ==============================================================================
# plot_tools.R
# Display functions only — no computation, no physics.
# Used by all plot files regardless of system or convolution type.
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

# Panel widths: [row label, heatmap, Wigner, Husimi]
PANEL_WIDTHS    <- c(0.13, 1, 1, 1)
FIGURE_WIDTH_IN <- 7.0
ROW_HEIGHT_IN   <- 1.8
FIGURE_PAD_IN   <- 0.5

ROW_LABEL_SIZE  <- 4.5

# ------------------------------------------------------------------------------
# ROW LABEL
# ------------------------------------------------------------------------------

plot_row_label <- function(label_str, parse=TRUE, base_font="") {
  ggplot() +
    theme_void(base_family=base_font) +
    coord_cartesian(xlim=c(0,1), ylim=c(0,1), clip="off") +
    annotate("text", x=1.0, y=0.5, label=label_str, parse=parse,
             family=base_font, size=ROW_LABEL_SIZE, hjust=1, vjust=0.5) +
    theme(plot.margin=margin(0, 0, 0, 0))
}

# ------------------------------------------------------------------------------
# PHASE SPACE HEATMAP — Husimi variant
# Symmetric diverging colormap centered at 0. No power compression: positive
# and negative regions show their amplitude faithfully so oscillation
# structure is visible without blurring high-amplitude lobes together.
# ------------------------------------------------------------------------------

plot_phase_space_heatmap_husimi <- function(dt_w2d, husimi_ell, df_traj=NULL,
                                            q_lim, p_lim,
                                            custom_breaks_q, custom_breaks_p,
                                            label_format, base_font="") {
  ax_x <- expression(italic(q)/italic(q)[0])
  ax_y <- expression(italic(p)/italic(p)[0])

  p <- ggplot(dt_w2d, aes(x=q, y=p)) +
    geom_raster(aes(fill=w_plot), interpolate=TRUE) +
    scale_fill_gradient2(low="gray20", mid="white", high="gray20",
                         midpoint=0, limits=c(-1,1), guide="none")

  if (!is.null(df_traj))
    p <- p + geom_path(data=df_traj, aes(x=q, y=p), inherit.aes=FALSE,
                       color="black", linewidth=0.5, linetype="solid")

  p +
    geom_circle(data=husimi_ell$circle, aes(x0=x0, y0=y0, r=r),
                inherit.aes=FALSE,
                color="gray20", linewidth=0.4, linetype="dashed") +
    coord_fixed(xlim=q_lim, ylim=p_lim, expand=FALSE) +
    scale_x_continuous(breaks=custom_breaks_q, labels=label_format) +
    scale_y_continuous(breaks=custom_breaks_p, labels=label_format) +
    theme_bw(base_family=base_font) +
    theme(panel.grid.minor=element_blank(),
          panel.background=element_rect(fill="white"),
          axis.text=element_text(size=8),
          plot.margin=margin(2,4,2,4)) +
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
# COLUMN TITLES
# ------------------------------------------------------------------------------

add_column_titles <- function(p_label, p_left, p_center, p_right,
                              title_left, title_center, title_right) {
  list(
    p_label  + labs(title=" ")         + theme(plot.title=element_text(size=11, hjust=0.5)),
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

    p_label <- plot_row_label(label, base_font=base_font)

    if (i == 1) {
      panels <- add_column_titles(p_label, p_left, p_center, p_right,
                                  COLUMN_TITLE_LEFT,
                                  COLUMN_TITLE_CENTER,
                                  COLUMN_TITLE_RIGHT)
      p_label <- panels[[1]]; p_left <- panels[[2]]
      p_center <- panels[[3]]; p_right <- panels[[4]]
    }
    if (i != num_rows) {
      lst <- suppress_x_titles(p_left, p_center, p_right)
      p_left <- lst[[1]]; p_center <- lst[[2]]; p_right <- lst[[3]]
    }
    if (i != ceiling(num_rows/2)) {
      lst <- suppress_y_titles(p_left, p_center, p_right)
      p_left <- lst[[1]]; p_center <- lst[[2]]; p_right <- lst[[3]]
    }

    plot_list <- c(plot_list, list(p_label, p_left, p_center, p_right))
  }

  wrap_plots(plot_list, ncol=4, widths=PANEL_WIDTHS) +
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
