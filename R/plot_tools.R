# ==============================================================================
# plot_tools.R
# Display functions only — no computation, no physics.
#
# Supports two phase-space density types:
#   Wigner   — signed (negative possible), uses asymmetric diverging colormap
#   Semiclassical — non-negative, uses single-direction grayscale
#
# Layout architecture:
#   Three columns: heatmap, cross-section/caustic, kernel resolution.
#   Row labels are attached as patchwork plot tags.
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

COLUMN_TITLE_LEFT                 <- "Phase-Space Cells"
COLUMN_TITLE_CENTER_WIGNER        <- "Wigner Negativity"
COLUMN_TITLE_CENTER_SEMICLASSICAL <- "Semiclassical Caustics"
COLUMN_TITLE_RIGHT_HUSIMI         <- "Husimi Resolution"
COLUMN_TITLE_RIGHT_SYMPLECTIC     <- "Symplectic Resolution"

PANEL_WIDTHS    <- c(1, 1, 1)
FIGURE_WIDTH_IN <- 7.5
ROW_HEIGHT_IN   <- 1.8
FIGURE_PAD_IN   <- 0.5

HEATMAP_LEFT_MARGIN_PT <- 80

# Wigner colormap: asymmetric diverging (gray45 negative -> white zero -> gray10 positive)
HEATMAP_COLOR_NEG  <- "gray45"
HEATMAP_COLOR_ZERO <- "white"
HEATMAP_COLOR_POS  <- "gray10"

# Semiclassical colormap: non-negative, single direction
HEATMAP_COLOR_LOW  <- "white"
HEATMAP_COLOR_HIGH <- "gray10"

ROW_LABEL_SIZE_PT <- 11
TAG_X_NPC <- -0.18
TAG_Y_NPC <-  0.5

# ------------------------------------------------------------------------------
# WIGNER HEATMAP (signed density, diverging colormap)
# ------------------------------------------------------------------------------

plot_wigner_heatmap <- function(dt_w2d, overlay_layers, df_traj=NULL,
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

  for (layer in overlay_layers) p <- p + layer

  p +
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
# SEMICLASSICAL HEATMAP (non-negative density, single-direction colormap)
# Renders the regularized energy shell. No classical-orbit overlay since
# the heatmap *is* the orbit (regularized).
# ------------------------------------------------------------------------------

plot_semiclassical_heatmap <- function(dt_w2d, overlay_layers,
                                       q_lim, p_lim,
                                       custom_breaks_q, custom_breaks_p,
                                       label_format, base_font="") {
  ax_x <- expression(italic(q)/italic(q)[0])
  ax_y <- expression(italic(p)/italic(p)[0])

  p <- ggplot(dt_w2d, aes(x=q, y=p)) +
    geom_raster(aes(fill=w_plot), interpolate=TRUE) +
    scale_fill_gradient(low=HEATMAP_COLOR_LOW,
                        high=HEATMAP_COLOR_HIGH,
                        limits=c(0,1), guide="none")

  for (layer in overlay_layers) p <- p + layer

  p +
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
# WIGNER CROSS-SECTION
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

# ------------------------------------------------------------------------------
# WKB CAUSTIC CROSS-SECTION (with infinity-arrow rendering at turning points)
#
# Renders the analytical WKB density 1/(T * sqrt(2*(E - V(q)))) which
# diverges at classical turning points. The divergence is rendered as a
# full-height ribbon at the turning point with arrow markers and an
# infinity symbol, following the convention from Berry's textbooks.
# ------------------------------------------------------------------------------

plot_wkb_caustic_cross_section <- function(dt, q_lim, y_lim, custom_breaks,
                                           label_format, q_minus, q_plus,
                                           base_font="") {
  ax_x <- expression(italic(q)/italic(q)[0])
  ax_y <- expression(italic(P)[WKB](italic(q)))

  # Determine clip threshold for the line: anything above 0.85*y_lim hides
  # behind the turning-point ribbon
  clip_y <- y_lim * 0.85

  # Interior of allowed region with finite caustic values
  dt_interior <- dt[is.finite(wkb_density) & wkb_density <= clip_y &
                      q > q_minus & q < q_plus]
  # Capped at clip_y for path rendering near turning points
  dt_capped <- dt[is.finite(wkb_density) & q > q_minus & q < q_plus]
  dt_capped[, wkb_capped := pmin(wkb_density, clip_y)]

  ggplot() +
    geom_hline(yintercept=0, color="black", linewidth=0.3) +
    # Fill the bowl
    geom_ribbon(data=dt_capped, aes(x=q, ymin=0, ymax=wkb_capped),
                fill="gray85", color=NA) +
    # Vertical infinity bars at turning points (full height to top of plot)
    annotate("segment", x=q_minus, xend=q_minus, y=0, yend=y_lim,
             color="black", linewidth=0.4,
             arrow=arrow(length=unit(0.12,"cm"), ends="last", type="closed")) +
    annotate("segment", x=q_plus, xend=q_plus, y=0, yend=y_lim,
             color="black", linewidth=0.4,
             arrow=arrow(length=unit(0.12,"cm"), ends="last", type="closed")) +
    # Infinity annotations near top
    annotate("text", x=q_minus, y=y_lim*0.92, label="infinity", parse=TRUE,
             color="black", size=4, hjust=-0.3, family=base_font) +
    annotate("text", x=q_plus, y=y_lim*0.92, label="infinity", parse=TRUE,
             color="black", size=4, hjust=1.3, family=base_font) +
    # Caustic line in the interior
    geom_path(data=dt_interior, aes(x=q, y=wkb_density),
              color="black", linewidth=0.4) +
    coord_cartesian(xlim=q_lim, ylim=c(0,y_lim), expand=FALSE) +
    scale_x_continuous(breaks=custom_breaks, labels=label_format) +
    theme_bw(base_family=base_font) +
    theme(panel.grid.minor=element_blank(),
          axis.text=element_text(size=8),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          aspect.ratio=1, plot.margin=margin(2,4,2,4)) +
    labs(x=ax_x, y=ax_y)
}

# ------------------------------------------------------------------------------
# HUSIMI / SYMPLECTIC CROSS-SECTIONS (Wigner pipeline)
# ------------------------------------------------------------------------------

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

plot_symplectic_cross_section <- function(dt, q_lim, y_lim, custom_breaks, label_format, base_font="") {
  ax_x <- expression(italic(q)/italic(q)[0])
  ax_y <- expression(italic(P)[italic(delta*q)](italic(q)*","*0))
  ggplot(dt, aes(x=q, y=P_sympl)) +
    geom_hline(yintercept=0, color="black", linewidth=0.3) +
    geom_ribbon(aes(ymin=0, ymax=pmax(P_sympl,0)), fill="gray85", color=NA) +
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
# SEMICLASSICAL RIGHT-COLUMN: rho_{delta q}(q)
# 1D position density obtained by marginalizing the symplectic-convolved
# 2D shell over p. Always non-negative.
# ------------------------------------------------------------------------------

plot_semiclassical_resolution <- function(dt, q_lim, y_lim, custom_breaks, label_format, base_font="") {
  ax_x <- expression(italic(q)/italic(q)[0])
  ax_y <- expression(rho[italic(delta*q)](italic(q)))
  ggplot(dt, aes(x=q, y=rho_sympl)) +
    geom_hline(yintercept=0, color="black", linewidth=0.3) +
    geom_ribbon(aes(ymin=0, ymax=pmax(rho_sympl,0)), fill="gray85", color=NA) +
    geom_path(color="black", linewidth=0.4) +
    coord_cartesian(xlim=q_lim, ylim=c(0,y_lim), expand=FALSE) +
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
# GRID ASSEMBLY (now takes title_center as well as title_right)
# ------------------------------------------------------------------------------

assemble_grid <- function(rows, title_center, title_right, base_font="") {
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
                                  title_center,
                                  title_right)
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
