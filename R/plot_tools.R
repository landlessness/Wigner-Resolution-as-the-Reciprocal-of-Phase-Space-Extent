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

# Compass-figure-specific labels (used by assemble_compass.R).
# The compass figure has a different column structure than the original
# n-indexed Morse/harmonic figures, so it gets its own constants.
COMPASS_COLUMN_TITLE_LEFT   <- "Phase-Space Cells"
COMPASS_COLUMN_TITLE_CENTER <- "Cross Sections"
COMPASS_COLUMN_TITLE_RIGHT  <- "Resolved"
COMPASS_ROW_LABEL_TOP       <- "Wigner"
COMPASS_ROW_LABEL_MIDDLE    <- "Husimi"
COMPASS_ROW_LABEL_BOTTOM    <- "Symplectic"

PANEL_WIDTHS    <- c(1, 1, 1)
FIGURE_WIDTH_IN <- 7.5
ROW_HEIGHT_IN   <- 1.8
FIGURE_PAD_IN   <- 0.5

HEATMAP_LEFT_MARGIN_PT <- 80

# Wigner colormap: asymmetric diverging (gray20 negative -> white zero -> gray5 positive)
HEATMAP_COLOR_NEG  <- "gray20"
HEATMAP_COLOR_ZERO <- "white"
HEATMAP_COLOR_POS  <- "gray5"

# Cross-section ribbon fill — used to gently highlight Wigner negativity
# in the 1D cross-section panels. Decoupled from the heatmap saturation
# constants because the visual weight needs are different: the heatmap
# packs many pixels at saturation, while a cross-section ribbon is a
# single sliver and needs to be subtle to avoid overwhelming the curve.
CROSS_RIBBON_NEG_FILL <- "gray80"

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
                       color="gray30", linewidth=0.5, linetype="solid")

  for (layer in overlay_layers) p <- p + layer

  p +
    coord_fixed(xlim=q_lim, ylim=p_lim, expand=FALSE) +
    scale_x_continuous(breaks=custom_breaks_q, labels=label_format) +
    scale_y_continuous(breaks=custom_breaks_p, labels=label_format) +
    theme_bw(base_family=base_font) +
    theme(panel.grid.minor=element_blank(),
          panel.background=element_rect(fill="white"),
          axis.text=element_text(size=8),
          plot.margin=margin(2, 2, 2, 2)) +
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
          plot.margin=margin(2, 2, 2, 2)) +
    labs(x=ax_x, y=ax_y)
}

# ------------------------------------------------------------------------------
# CROSS-SECTION FILL CONVENTION
#
# The 1D cross-section panels visually act as a saturation legend for the
# 2D heatmap. Two solid fill colors:
#   positive ribbon (between zero and a positive curve value): HEATMAP_COLOR_POS
#   negative ribbon (between zero and a negative curve value): HEATMAP_COLOR_NEG
# These are exactly the saturation endpoints of the heatmap diverging
# colormap. A cross-section ribbon at full positive height appears the same
# color as a 2D pixel at full positive saturation.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# WIGNER CROSS-SECTION (signed, two-color fill matching heatmap endpoints)
# ------------------------------------------------------------------------------

plot_wigner_cross_section <- function(dt, q_lim, y_lim, custom_breaks, label_format, base_font="") {
  ax_x <- expression(italic(q)/italic(q)[0])
  ax_y <- expression(italic(W)[italic(n)](italic(q)*","*0))
  # Recover the data peak from y_lim (the convention is y_lim = peak*1.1).
  peak <- y_lim / 1.1
  y_breaks <- c(-peak, 0, peak)
  y_labels <- sprintf("%.2f", y_breaks)
  ggplot(dt, aes(x=q, y=W_raw)) +
    geom_ribbon(aes(ymin=pmin(W_raw,0), ymax=0),
                fill=CROSS_RIBBON_NEG_FILL, color=NA) +
    geom_path(color="black", linewidth=0.4) +
    coord_cartesian(xlim=q_lim, ylim=c(-y_lim,y_lim), expand=FALSE) +
    scale_x_continuous(breaks=custom_breaks, labels=label_format) +
    scale_y_continuous(breaks=y_breaks, labels=y_labels) +
    theme_bw(base_family=base_font) +
    theme(panel.grid.minor=element_blank(),
          axis.text=element_text(size=8),
          aspect.ratio=1, plot.margin=margin(2,2,2,2)) +
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
          aspect.ratio=1, plot.margin=margin(2,2,2,2)) +
    labs(x=ax_x, y=ax_y)
}

# ------------------------------------------------------------------------------
# HUSIMI / SYMPLECTIC CROSS-SECTIONS (Wigner pipeline)
# ------------------------------------------------------------------------------

plot_husimi_cross_section <- function(dt, q_lim, y_lim, custom_breaks, label_format, base_font="") {
  ax_x <- expression(italic(q)/italic(q)[0])
  ax_y <- expression(italic(Q)(italic(q)*","*0))
  # Recover the data peak from y_lim (the convention is y_lim = peak*1.1).
  peak <- y_lim / 1.1
  y_breaks <- c(0, peak)
  y_labels <- sprintf("%.2g", y_breaks)
  ggplot(dt, aes(x=q, y=Q_husimi)) +
    geom_path(color="black", linewidth=0.4) +
    coord_cartesian(xlim=q_lim, ylim=c(-y_lim,y_lim), expand=FALSE) +
    scale_x_continuous(breaks=custom_breaks, labels=label_format) +
    scale_y_continuous(breaks=y_breaks, labels=y_labels) +
    theme_bw(base_family=base_font) +
    theme(panel.grid.minor=element_blank(),
          axis.text=element_text(size=8),
          aspect.ratio=1, plot.margin=margin(2,2,2,2)) +
    labs(x=ax_x, y=ax_y)
}

plot_symplectic_cross_section <- function(dt, q_lim, y_lim, custom_breaks, label_format, base_font="") {
  ax_x <- expression(italic(q)/italic(q)[0])
  ax_y <- expression(italic(P)[italic(delta*q)](italic(q)*","*0))
  # Recover the data peak from y_lim (the convention is y_lim = peak*1.1).
  peak <- y_lim / 1.1
  y_breaks <- c(0, peak)
  y_labels <- sprintf("%.2g", y_breaks)
  ggplot(dt, aes(x=q, y=P_sympl)) +
    geom_path(color="black", linewidth=0.4) +
    coord_cartesian(xlim=q_lim, ylim=c(-y_lim,y_lim), expand=FALSE) +
    scale_x_continuous(breaks=custom_breaks, labels=label_format) +
    scale_y_continuous(breaks=y_breaks, labels=y_labels) +
    theme_bw(base_family=base_font) +
    theme(panel.grid.minor=element_blank(),
          axis.text=element_text(size=8),
          aspect.ratio=1, plot.margin=margin(2,2,2,2)) +
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
    geom_ribbon(aes(ymin=0, ymax=pmax(rho_sympl,0)),
                fill=HEATMAP_COLOR_POS, color=NA) +
    geom_path(color="black", linewidth=0.4) +
    coord_cartesian(xlim=q_lim, ylim=c(0,y_lim), expand=FALSE) +
    scale_x_continuous(breaks=custom_breaks, labels=label_format) +
    theme_bw(base_family=base_font) +
    theme(panel.grid.minor=element_blank(),
          axis.text=element_text(size=8),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          aspect.ratio=1, plot.margin=margin(2,2,2,2)) +
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
    geom_ribbon(aes(ymin=0, ymax=pmax(rho_husimi,0)),
                fill=HEATMAP_COLOR_POS, color=NA) +
    geom_path(color="black", linewidth=0.4) +
    coord_cartesian(xlim=q_lim, ylim=c(0,y_lim), expand=FALSE) +
    scale_x_continuous(breaks=custom_breaks, labels=label_format) +
    theme_bw(base_family=base_font) +
    theme(panel.grid.minor=element_blank(),
          axis.text=element_text(size=8),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          aspect.ratio=1, plot.margin=margin(2,2,2,2)) +
    labs(x=ax_x, y=ax_y)
}

# ------------------------------------------------------------------------------
# COLORMAP LEGEND PANEL
#
# Builds a small ggplot rendering the diverging Wigner colormap as a vertical
# colorbar. Used to fill an otherwise-empty panel slot in assembled figures
# (e.g., the upper-right of compass.pdf) and to tell the reader what the
# heatmap shading actually represents in physical units.
#
# Tick labels show the physical scale (max_abs at top, 0 in the middle,
# -max_abs at the bottom) so a reader looking at any heatmap can read off
# the value at any point.
# ------------------------------------------------------------------------------

#' Build a vertical colormap-legend ggplot panel for the Wigner pipeline.
#'
#' @param max_abs   Physical peak magnitude for tick labels at +/- max_abs.
#' @param y_label   Optional axis label (e.g., "W(q,p)"). Pass NULL for none.
#' @param base_font Font family for the panel.
#' @return A ggplot object suitable for placement in a patchwork layout.
plot_compass_legend <- function(max_abs, y_label=NULL, base_font="") {
  # 200 horizontal strips covering [-1, 1] in normalized coordinates.
  n_strips <- 200
  y_norm   <- seq(-1, 1, length.out=n_strips)
  # geom_rect xmin/xmax/ymin/ymax build a thin column.
  dy <- diff(y_norm)[1]
  rects <- data.table(
    xmin = 0,                 xmax = 1,
    ymin = y_norm - dy/2,     ymax = y_norm + dy/2,
    fill_value = y_norm
  )

  # Map normalized [-1, 1] back to physical labels.
  fmt <- function(z) sprintf("%.2f", z * max_abs)
  tick_norm   <- c(-1, 0, 1)
  tick_labels <- c(fmt(-1), "0", fmt(1))

  ax_y_expr <- if (is.null(y_label)) NULL else parse(text=y_label)

  ggplot(rects) +
    geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                  fill=fill_value), color=NA) +
    scale_fill_gradient2(low=HEATMAP_COLOR_NEG,
                         mid=HEATMAP_COLOR_ZERO,
                         high=HEATMAP_COLOR_POS,
                         midpoint=0, limits=c(-1, 1), guide="none") +
    scale_x_continuous(limits=c(0, 1), breaks=NULL, expand=c(0, 0)) +
    scale_y_continuous(limits=c(-1, 1), breaks=tick_norm, labels=tick_labels,
                       expand=c(0, 0), position="right") +
    coord_fixed(ratio=0.25) +    # tall narrow bar; ratio = width/height
    theme_bw(base_family=base_font) +
    theme(panel.grid=element_blank(),
          panel.background=element_rect(fill="white"),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.text.y=element_text(size=8),
          plot.margin=margin(2, 2, 2, 2)) +
    labs(x=NULL, y=ax_y_expr)
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

#' Attach a plain-text row label (rotated 90 deg) to the leftmost panel of
#' a compass figure row. Used by assemble_compass.R for "Wigner", "Husimi",
#' "Symplectic" labels that are too long to display horizontally.
#'
#' Differs from attach_row_tag() in three ways:
#'   - label is plain text, not parsed as a math expression
#'   - rotated 90 deg (reads bottom-to-top on the left side)
#'   - smaller left margin (the rotated text is narrow)
#'
#' @param p           Ggplot object (typically the leftmost panel of a row).
#' @param label_str   Plain text label string.
#' @param base_font   Font family.
#' @return Modified ggplot object.
attach_compass_row_tag <- function(p, label_str, base_font="") {
  p +
    labs(tag = label_str) +
    theme(
      plot.tag = element_text(size=ROW_LABEL_SIZE_PT, family=base_font,
                              hjust=0.5, vjust=0.5, angle=90),
      plot.tag.position = c(-0.22, 0.5),
      plot.margin = margin(2, 4, 2, 36)
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
