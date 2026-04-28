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
# Detects all classical turning points from the data (positions where the
# WKB density transitions between positive-finite and zero/Inf along q),
# then draws:
#   - Bowl ribbon fills only in classically-allowed regions
#   - Vertical infinity arrows at every turning point
#   - "infinity" text annotation outside each vertical line
#
# Handles both:
#   - 2 turning points (above-barrier / single-well states): one bowl
#   - 4 turning points (sub-barrier double-well states):    two bowls
# ------------------------------------------------------------------------------

#' Detect turning point positions from a wkb_density vector on a q grid.
#' Turning points are where the density transitions between positive-finite
#' (classically allowed) and zero (forbidden) cells along q.
#' Returns sorted q-values of all detected transitions.
detect_wkb_turning_points <- function(q, wkb_density) {
  # A cell is "allowed" if density > 0 and finite
  allowed <- is.finite(wkb_density) & wkb_density > 0
  # Detect transitions: where allowed[i] != allowed[i+1]
  transitions <- which(diff(allowed) != 0)
  if (length(transitions) == 0) return(numeric(0))
  # Place the turning point at the boundary between cells
  tp_q <- (q[transitions] + q[transitions + 1]) / 2
  sort(tp_q)
}

#' Identify contiguous classically-allowed regions for ribbon clipping.
#' Returns a list of (start_idx, end_idx) pairs, one per allowed segment.
identify_allowed_segments <- function(wkb_density) {
  allowed <- is.finite(wkb_density) & wkb_density > 0
  if (!any(allowed)) return(list())
  rle_a <- rle(allowed)
  ends   <- cumsum(rle_a$lengths)
  starts <- ends - rle_a$lengths + 1
  segments <- list()
  for (k in seq_along(rle_a$values)) {
    if (rle_a$values[k]) {
      segments[[length(segments)+1]] <- list(start=starts[k], end=ends[k])
    }
  }
  segments
}

plot_wkb_caustic_cross_section <- function(dt, q_lim, y_lim, custom_breaks,
                                           label_format, base_font="", ...) {
  ax_x <- expression(italic(q)/italic(q)[0])
  ax_y <- expression(italic(P)[WKB](italic(q)))

  # Auto-detect turning points from the data
  turning_points <- detect_wkb_turning_points(dt$q, dt$wkb_density)

  # Identify allowed segments for ribbon clipping
  segments <- identify_allowed_segments(dt$wkb_density)

  # The line and ribbon should run all the way to y_lim where they meet
  # the vertical infinity arrows seamlessly
  clip_y <- y_lim

  # Build the plot: start with the zero baseline
  p <- ggplot() +
    geom_hline(yintercept=0, color="black", linewidth=0.3)

  # Add ribbon and line for each allowed segment separately, so the
  # forbidden regions render as white space
  for (seg in segments) {
    dt_seg <- dt[seg$start:seg$end]
    dt_seg[, wkb_capped := pmin(wkb_density, clip_y)]
    p <- p +
      geom_ribbon(data=dt_seg[is.finite(wkb_density)],
                  aes(x=q, ymin=0, ymax=wkb_capped),
                  fill="gray85", color=NA) +
      geom_path(data=dt_seg[is.finite(wkb_density) & wkb_density <= clip_y],
                aes(x=q, y=wkb_density),
                color="black", linewidth=0.4)
  }

  # Vertical infinity arrows at every turning point
  for (tp_q in turning_points) {
    p <- p + annotate("segment", x=tp_q, xend=tp_q, y=0, yend=y_lim,
                      color="black", linewidth=0.4,
                      arrow=arrow(length=unit(0.12,"cm"),
                                  ends="last", type="closed"))
  }

  # Infinity text annotations: place each just outside its vertical line.
  # For the leftmost arrow, label goes left of it; rightmost, right of it.
  # For inner pairs (sub-barrier case), labels go on the outer side of
  # each so they appear in the forbidden barrier region.
  if (length(turning_points) >= 2) {
    n_tp <- length(turning_points)
    for (k in seq_along(turning_points)) {
      tp_q <- turning_points[k]
      # Determine "outside" direction:
      # leftmost (k==1): outside = left  -> hjust = 1.3
      # rightmost (k==n): outside = right -> hjust = -0.3
      # inner-left (k==2 of 4): outside = right (toward forbidden zone) -> hjust = -0.3
      # inner-right (k==3 of 4): outside = left  -> hjust = 1.3
      if (k == 1)            hj <- 1.3
      else if (k == n_tp)    hj <- -0.3
      else if (k == 2 && n_tp == 4) hj <- -0.3
      else if (k == 3 && n_tp == 4) hj <- 1.3
      else                   hj <- 1.3
      p <- p + annotate("text", x=tp_q, y=y_lim*0.92,
                        label="infinity", parse=TRUE,
                        color="black", size=4, hjust=hj, family=base_font)
    }
  }

  p +
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
# SEMICLASSICAL RIGHT-COLUMN (Husimi variant): rho_Q(q)
# 1D position density obtained by marginalizing the Husimi-convolved
# 2D shell over p. Always non-negative.
# ------------------------------------------------------------------------------

plot_semiclassical_husimi_resolution <- function(dt, q_lim, y_lim, custom_breaks, label_format, base_font="") {
  ax_x <- expression(italic(q)/italic(q)[0])
  ax_y <- expression(rho[italic(Q)](italic(q)))
  ggplot(dt, aes(x=q, y=rho_husimi)) +
    geom_hline(yintercept=0, color="black", linewidth=0.3) +
    geom_ribbon(aes(ymin=0, ymax=pmax(rho_husimi,0)), fill="gray85", color=NA) +
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
# GRID ASSEMBLY
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
