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
COLUMN_TITLE_LEFT_ORBIT           <- "Classical Orbit"
COLUMN_TITLE_CENTER_WIGNER        <- "Wigner Negativity"
COLUMN_TITLE_CENTER_SEMICLASSICAL <- "Semiclassical Caustics"
COLUMN_TITLE_RIGHT_HUSIMI         <- "Husimi Resolution"
COLUMN_TITLE_RIGHT_SYMPLECTIC     <- "Symplectic Resolution"
COLUMN_TITLE_RIGHT_AIRY           <- "Airy Approximation"

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

# Left margin reserved for the row-label tag. Sized to fit the longest
# expected label (e.g. "25.38 A_0") at ROW_LABEL_SIZE_PT in a serif body
# font, with a small gutter to the panel.
HEATMAP_LEFT_MARGIN_PT <- 65

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

# Semiclassical colormap: non-negative, single direction. The "high" stop
# is a mid-gray rather than near-black so that the QoA overlay lines (drawn
# in true black) read clearly on top of the energy shell. The shell is
# the *trajectory* — secondary to the QoA structure, which is the figure's
# physical message — so it visually recedes.
HEATMAP_COLOR_LOW  <- "white"
HEATMAP_COLOR_HIGH <- "gray55"

# Semiclassical 1D-density ribbon fill (middle WKB caustic and right
# symplectic-resolved density). Matches the conventional gray-fill /
# black-line look for non-negative probability densities; both panels
# share this single constant so they stay visually consistent.
SEMICLASSICAL_RIBBON_FILL <- "gray60"

ROW_LABEL_SIZE_PT <- 11
# Row-label tag is RIGHT-aligned (hjust=1 in attach_row_tag), with its
# right edge anchored at TAG_X_NPC (npc within the panel). With right
# alignment, the *trailing* element of the label (e.g. "A_0") column-
# aligns vertically across rows regardless of label length, while the
# leading numeric digit grows leftward into the margin gutter.
TAG_X_NPC <- -0.05
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
  ax_y <- expression(group("|", italic(psi)[WKB](italic(q)), "|")^2)

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
                  fill=SEMICLASSICAL_RIBBON_FILL, color=NA) +
      geom_path(data=dt_seg[is.finite(wkb_density) & wkb_density <= clip_y],
                aes(x=q, y=wkb_density),
                color="black", linewidth=0.4)
  }

  # Vertical infinity arrows at turning points. Special handling: in the
  # 4-turning-point sub-barrier case, if the inner pair is closer than a
  # fraction of the panel width, two ∞ text labels would clash visually.
  # In that case we still draw arrows at every turning point (the arrows
  # mark physical caustic locations and cannot be moved), but collapse
  # the inner pair's two labels into a single "∞" glyph centered between
  # them in the forbidden barrier region. This is a typography choice for
  # readability, not a physics one.
  inner_label_clash <- length(turning_points) == 4 &&
    (turning_points[3] - turning_points[2]) < 0.20 * diff(q_lim)

  # Arrows: always at every turning point (physical caustic locations).
  for (tp_q in turning_points) {
    p <- p + annotate("segment", x=tp_q, xend=tp_q, y=0, yend=y_lim,
                      color="black", linewidth=0.4,
                      arrow=arrow(length=unit(0.12,"cm"),
                                  ends="last", type="closed"))
  }

  # Labels: usually at every turning point with hjust pushing the label
  # outside (left or right) so it doesn't sit on the arrow itself.
  # Exception: when the inner pair is too close to fit two labels, we
  # emit a single centered label between them and skip the per-tp labels
  # for tp[2] and tp[3].
  if (inner_label_clash) {
    # Outer two labels as usual
    p <- p + annotate("text", x=turning_points[1], y=y_lim*0.92,
                      label="infinity", parse=TRUE,
                      color="black", size=4, hjust=1.3, family=base_font)
    p <- p + annotate("text", x=turning_points[4], y=y_lim*0.92,
                      label="infinity", parse=TRUE,
                      color="black", size=4, hjust=-0.3, family=base_font)
    # One centered ∞ label for the merged inner pair, sitting between
    # the two inner arrows in the forbidden barrier region.
    barrier_q <- 0.5 * (turning_points[2] + turning_points[3])
    p <- p + annotate("text", x=barrier_q, y=y_lim*0.92,
                      label="infinity", parse=TRUE,
                      color="black", size=4, hjust=0.5, family=base_font)
  } else if (length(turning_points) >= 2) {
    n_tp <- length(turning_points)
    for (k in seq_along(turning_points)) {
      tp_q <- turning_points[k]
      if (k == 1)                       hj <- 1.3
      else if (k == n_tp)               hj <- -0.3
      else if (k == 2 && n_tp == 4)     hj <- -0.3
      else if (k == 3 && n_tp == 4)     hj <- 1.3
      else                              hj <- 1.3
      p <- p + annotate("text", x=tp_q, y=y_lim*0.92,
                        label="infinity", parse=TRUE,
                        color="black", size=4, hjust=hj, family=base_font)
    }
  }

  # Y-tick labels for the (possibly oscillating) WKB density:
  # - 0 at bottom (the curve dips to 0 at every cos^2 zero for the
  #   oscillating density; for the smooth caustic, 0 is also the
  #   forbidden-region value).
  # - the peak amplitude of the curve below y_lim (the largest finite
  #   value the renderer will actually show — gives the reader a
  #   concrete amplitude scale before the curve runs into the
  #   infinity arrows at the turning points).
  finite_dens <- dt$wkb_density[is.finite(dt$wkb_density) & dt$wkb_density > 0
                                & dt$wkb_density < y_lim]
  peak_below  <- if (length(finite_dens) > 0) max(finite_dens, na.rm=TRUE) else NA
  if (is.finite(peak_below)) {
    y_breaks <- c(0, peak_below)
    y_labels <- c("0", sprintf("%.2g", peak_below))
  } else {
    y_breaks <- c(0)
    y_labels <- c("0")
  }

  p +
    coord_cartesian(xlim=q_lim, ylim=c(0,y_lim), expand=FALSE) +
    scale_x_continuous(breaks=custom_breaks, labels=label_format) +
    scale_y_continuous(breaks=y_breaks, labels=y_labels) +
    theme_bw(base_family=base_font) +
    theme(panel.grid.minor=element_blank(),
          axis.text=element_text(size=8),
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

#' Right-column symplectic-resolved density, with optional comparator
#' overlay curves.
#'
#' Each overlay is drawn on top of the ribbon and main symplectic curve.
#' Used to compare the symplectic resolution against ground-truth
#' (Schroedinger exact) and prior-art (Airy uniform) densities.
#'
#' @param dt           data.table with columns q, rho_sympl
#' @param q_lim,y_lim  axis limits
#' @param custom_breaks,label_format  x-axis ticks
#' @param base_font    font family
#' @param overlays     OPTIONAL list of named lists, each with fields:
#'                       data:      data.frame with columns q, rho
#'                       linetype:  ggplot linetype string ("dashed", ...)
#'                       color:     line color
#'                       linewidth: line weight
plot_semiclassical_resolution <- function(dt, q_lim, y_lim, custom_breaks,
                                          label_format, base_font="",
                                          overlays=NULL,
                                          y_label=NULL) {
  ax_x <- expression(italic(q)/italic(q)[0])
  # Default y-label is the marginal rho_{delta q}(q). Callers plotting a
  # cross-section P_{delta q}(q, 0) — e.g. the Wigner-symplectic right
  # column — should pass an explicit cross-section y_label.
  ax_y <- if (is.null(y_label)) expression(rho[italic(delta*q)](italic(q))) else y_label
  # Y-tick label uses the maximum of all rendered curves (symplectic +
  # any overlays), so the displayed peak number actually reflects what
  # the reader sees on the plot.
  rho_peak <- max(dt$rho_sympl, na.rm=TRUE)
  if (!is.null(overlays)) {
    for (ov in overlays) rho_peak <- max(rho_peak, max(ov$data$rho, na.rm=TRUE))
  }
  y_breaks <- c(0, rho_peak)
  y_labels <- c("0", sprintf("%.2g", rho_peak))
  p <- ggplot(dt, aes(x=q, y=rho_sympl)) +
    geom_ribbon(aes(ymin=0, ymax=pmax(rho_sympl,0)),
                fill=SEMICLASSICAL_RIBBON_FILL, color=NA)

  # Overlay fills are drawn after the symplectic ribbon (so they sit on
  # top of it via alpha-blending) but before the symplectic line and the
  # overlay lines (so those line elements remain on top of all fills).
  # Each overlay may optionally specify `fill` (color) and `fill_alpha`
  # (transparency 0-1); when omitted, no fill is drawn for that overlay.
  if (!is.null(overlays)) {
    for (ov in overlays) {
      if (!is.null(ov$fill)) {
        fa <- if (is.null(ov$fill_alpha)) 0.3 else ov$fill_alpha
        p <- p + geom_ribbon(data=ov$data,
                             aes(x=q, ymin=0, ymax=pmax(rho,0)),
                             inherit.aes=FALSE,
                             fill=ov$fill, color=NA, alpha=fa)
      }
    }
  }

  p <- p + geom_path(color="black", linewidth=0.4)

  if (!is.null(overlays)) {
    for (ov in overlays) {
      lt <- if (is.null(ov$linetype)) "solid" else ov$linetype
      p <- p + geom_path(data=ov$data, aes(x=q, y=rho),
                         inherit.aes=FALSE,
                         color=ov$color,
                         linewidth=ov$linewidth,
                         linetype=lt)
    }
  }

  p +
    coord_cartesian(xlim=q_lim, ylim=c(0,y_lim), expand=FALSE) +
    scale_x_continuous(breaks=custom_breaks, labels=label_format) +
    scale_y_continuous(breaks=y_breaks, labels=y_labels) +
    theme_bw(base_family=base_font) +
    theme(panel.grid.minor=element_blank(),
          axis.text=element_text(size=8),
          aspect.ratio=1, plot.margin=margin(2,2,2,2)) +
    labs(x=ax_x, y=ax_y)
}

# ------------------------------------------------------------------------------
# SPLIT RIGHT-COLUMN: SYMPLECTIC (top, 80%) + COMPARATOR (bottom, 20%)
#
# Vertically stacks two ribbon plots sharing the q-axis:
#   - top sub-panel: symplectic curve (the result), filled with the
#     standard SEMICLASSICAL_RIBBON_FILL
#   - bottom sub-panel: comparator curve (Husimi or Airy), filled with a
#     lighter gray to signal subordinate visual role
#
# Each sub-panel has its own y-axis with independent scaling — the
# comparison story is about shape, not amplitude (the prose handles
# amplitude commentary). The bottom sub-panel renders the x-axis text
# and title; the top sub-panel hides them so the two share a single
# x-axis legend at the bottom.
#
# Returns a patchwork composition. Behaves like a ggplot when passed
# through wrap_plots() in assemble_grid(); add_column_titles() and the
# x/y title suppression helpers operate on it correctly.
# ------------------------------------------------------------------------------

#' Split right-column panel: symplectic on top, comparator on bottom.
#'
#' Implemented as a single ggplot with `facet_grid` (rather than patchwork)
#' so all the existing layout machinery — column titles, axis-title
#' suppression, aspect.ratio — operates on it without special handling.
#' Each facet has its own y-axis range (via scales="free_y"); space="free"
#' allocates row heights proportional to `panel_heights`.
#'
#' @param dt_sympl    data.table with columns q, rho_sympl (the main result).
#' @param dt_comp     data.table with columns q, rho_comp (the comparator).
#' @param q_lim       Shared x-axis limits.
#' @param custom_breaks,label_format  X-axis tick spec.
#' @param y_label_top    Expression or string for the top sub-panel y-label.
#' @param y_label_bottom Expression or string for the bottom sub-panel y-label.
#' @param base_font   Font family.
#' @param comp_fill   Fill color for the comparator ribbon (default
#'                    "gray85" — lighter than SEMICLASSICAL_RIBBON_FILL).
#' @return A ggplot.
plot_semiclassical_resolution_split <- function(dt_sympl, dt_comp,
                                                q_lim,
                                                custom_breaks, label_format,
                                                y_label_top, y_label_bottom,
                                                base_font="",
                                                comp_fill="gray85") {
  ax_x <- expression(italic(q)/italic(q)[0])

  # Combine into a single data frame with a facet-key column. Each row's
  # rho holds the appropriate density (symplectic or comparator). The
  # facet column is an ordered factor with the symplectic facet first
  # (top); facet_grid will respect this ordering.
  df <- rbind(
    data.frame(q=dt_sympl$q, rho=dt_sympl$rho_sympl,
               facet="sympl", stringsAsFactors=FALSE),
    data.frame(q=dt_comp$q,  rho=dt_comp$rho_comp,
               facet="comp",  stringsAsFactors=FALSE)
  )
  df$facet <- factor(df$facet, levels=c("sympl", "comp"))

  # Per-facet fill mapping. We use a manual scale_fill_manual so each
  # facet gets its own ribbon fill color while sharing the same plot.
  df$fill_key <- df$facet

  # Per-facet y-axis labels via labeller. facet_grid uses strip text for
  # facet labels; we hide the strip and use a y-axis title strategy
  # instead — but facet_grid doesn't support per-facet y-axis titles
  # natively. Workaround: place the labels via strip.text (rotated) on
  # the left, with the strip styled to look like a y-axis label.
  facet_labels <- c(sympl=as.character(deparse(y_label_top)),
                    comp =as.character(deparse(y_label_bottom)))

  ggplot(df, aes(x=q, y=rho, fill=fill_key)) +
    geom_ribbon(aes(ymin=0, ymax=pmax(rho,0)), color=NA) +
    geom_path(color="black", linewidth=0.4) +
    facet_grid(facet ~ .,
               scales="free_y", space="fixed",
               switch="y",
               labeller=as_labeller(c(sympl=" ", comp=" "))) +
    scale_fill_manual(values=c(sympl=SEMICLASSICAL_RIBBON_FILL,
                                comp =comp_fill),
                      guide="none") +
    coord_cartesian(xlim=q_lim, expand=FALSE) +
    scale_x_continuous(breaks=custom_breaks, labels=label_format) +
    theme_bw(base_family=base_font) +
    theme(panel.grid.minor=element_blank(),
          axis.text=element_text(size=8),
          strip.background=element_blank(),
          strip.placement="outside",
          strip.text.y.left=element_blank(),
          panel.spacing.y=unit(2, "pt"),
          aspect.ratio=NULL,
          plot.margin=margin(2,2,2,2)) +
    labs(x=ax_x, y=NULL)
}

# ------------------------------------------------------------------------------
# SEMICLASSICAL RIGHT-COLUMN (Husimi variant): rho_Q(q)
# 1D position density obtained by marginalizing the Husimi-convolved
# 2D shell over p. Always non-negative.
# ------------------------------------------------------------------------------

plot_semiclassical_husimi_resolution <- function(dt, q_lim, y_lim, custom_breaks, label_format, base_font="") {
  ax_x <- expression(italic(q)/italic(q)[0])
  ax_y <- expression(rho[italic(Q)](italic(q)))
  rho_peak <- max(dt$rho_husimi, na.rm=TRUE)
  y_breaks <- c(0, rho_peak)
  y_labels <- c("0", sprintf("%.2g", rho_peak))
  ggplot(dt, aes(x=q, y=rho_husimi)) +
    geom_ribbon(aes(ymin=0, ymax=pmax(rho_husimi,0)),
                fill=SEMICLASSICAL_RIBBON_FILL, color=NA) +
    geom_path(color="black", linewidth=0.4) +
    coord_cartesian(xlim=q_lim, ylim=c(0,y_lim), expand=FALSE) +
    scale_x_continuous(breaks=custom_breaks, labels=label_format) +
    scale_y_continuous(breaks=y_breaks, labels=y_labels) +
    theme_bw(base_family=base_font) +
    theme(panel.grid.minor=element_blank(),
          axis.text=element_text(size=8),
          aspect.ratio=1, plot.margin=margin(2,2,2,2)) +
    labs(x=ax_x, y=ax_y)
}

# ------------------------------------------------------------------------------
# SEMICLASSICAL RIGHT-COLUMN (Airy variant): rho_Airy(q)
# 1D position density from the Berry-Mount/Langer uniform Airy
# approximation. Always non-negative; replaces the WKB caustic divergence
# at turning points with a finite Airy-function lobe. Same overlay
# protocol as plot_semiclassical_resolution: pass `overlays` to draw a
# Schroedinger ground-truth or other comparator on top.
# ------------------------------------------------------------------------------

plot_airy_resolution <- function(dt, q_lim, y_lim, custom_breaks,
                                 label_format, base_font="",
                                 overlays=NULL) {
  ax_x <- expression(italic(q)/italic(q)[0])
  ax_y <- expression(rho[Airy](italic(q)))
  rho_peak <- max(dt$rho_airy, na.rm=TRUE)
  if (!is.null(overlays)) {
    for (ov in overlays) rho_peak <- max(rho_peak, max(ov$data$rho, na.rm=TRUE))
  }
  y_breaks <- c(0, rho_peak)
  y_labels <- c("0", sprintf("%.2g", rho_peak))
  p <- ggplot(dt, aes(x=q, y=rho_airy)) +
    geom_ribbon(aes(ymin=0, ymax=pmax(rho_airy,0)),
                fill=SEMICLASSICAL_RIBBON_FILL, color=NA) +
    geom_path(color="black", linewidth=0.4)

  if (!is.null(overlays)) {
    for (ov in overlays) {
      lt <- if (is.null(ov$linetype)) "solid" else ov$linetype
      p <- p + geom_path(data=ov$data, aes(x=q, y=rho),
                         inherit.aes=FALSE,
                         color=ov$color,
                         linewidth=ov$linewidth,
                         linetype=lt)
    }
  }

  p +
    coord_cartesian(xlim=q_lim, ylim=c(0,y_lim), expand=FALSE) +
    scale_x_continuous(breaks=custom_breaks, labels=label_format) +
    scale_y_continuous(breaks=y_breaks, labels=y_labels) +
    theme_bw(base_family=base_font) +
    theme(panel.grid.minor=element_blank(),
          axis.text=element_text(size=8),
          aspect.ratio=1, plot.margin=margin(2,2,2,2)) +
    labs(x=ax_x, y=ax_y)
}

# ------------------------------------------------------------------------------
# CLASSICAL ORBIT PHASE-SPACE PANEL
# Black orbit on white. No heatmap. Same coord/axis convention as
# plot_semiclassical_heatmap so it slots into the same left-column
# position. Used by figures whose left column should show the
# deterministic classical orbit.
#
# Optional overlay_layers (e.g. the QoA cells from
# symplectic_overlay_layers()) draw on top of the orbit. The symplectic
# figure uses overlay_layers to show the action-capacity geometry; the
# Airy figure passes overlay_layers=NULL since the Airy method has no
# such geometry to show.
# ------------------------------------------------------------------------------

plot_classical_orbit_phase_space <- function(df_traj,
                                             q_lim, p_lim,
                                             custom_breaks_q, custom_breaks_p,
                                             label_format, base_font="",
                                             overlay_layers=NULL,
                                             orbit_color="black",
                                             orbit_linewidth=0.4) {
  ax_x <- expression(italic(q)/italic(q)[0])
  ax_y <- expression(italic(p)/italic(p)[0])

  p <- ggplot(df_traj, aes(x=q, y=p, group=group)) +
    geom_path(color=orbit_color, linewidth=orbit_linewidth, linetype="solid")

  if (!is.null(overlay_layers)) {
    for (layer in overlay_layers) p <- p + layer
  }

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

assemble_grid <- function(rows, title_center, title_right, base_font="",
                          title_left=COLUMN_TITLE_LEFT) {
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
                                  title_left,
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
