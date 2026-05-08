# ==============================================================================
# plot_heisenberg_cells.R
# Bare-bones SVG export for OmniGraffle: four rectangles representing the
# extended cell A and three squeezed cells (a_p at theta=0, a at theta=pi/4,
# a_q at theta=pi/2) in the action-capacity (QoA) ellipse family.
#
# This file produces an SVG containing ONLY four rectangle outlines on a
# transparent background. No axes, no gridlines, no labels, no fills, no
# panel border, no margins. The output is intended to be imported into
# OmniGraffle, where everything other than the four shapes is discarded
# anyway.
#
# Each rectangle's width and height are set from the ELLIPSE semi-axes of
# the corresponding cell:
#   width  = 2 * (long  semi-axis along the orientation)
#   height = 2 * (short semi-axis perpendicular to the orientation)
# i.e., the bounding box of the corresponding ellipse, drawn rotated to
# match the cell's orientation theta. (This is the "Heisenberg-cell as
# rectangle" convention common in textbook drawings.)
#
# Areas scale as 4 * L * S (rectangle area), not pi * L * S (ellipse area),
# so the rectangle areas come out 4/pi ~ 1.273 times the ellipse areas.
# The four shapes remain mutually self-consistent: every squeezed-cell
# rectangle has the SAME area (= 4 * hbar) by the polar-dual identity, and
# the extended-cell rectangle's area is 4 * Delta_q * Delta_p, with the
# same A/A_0 ratio as in the ellipse picture.
#
# Geometry (with hbar = 1):
#   Delta_q = 2.0,  Delta_p = 1.5  ->  A/A_0 = Delta_q*Delta_p / hbar = 3
#   so the outer rectangle is 3x the area of each inner rectangle. Tuned
#   smaller than the ellipse-version's 6x ratio so the squeezed cells fill
#   more of the extended cell visually.
#
# Long semi-axis at orientation theta (inscription / tangency):
#   L(theta) = ( cos^2(theta)/Delta_q^2 + sin^2(theta)/Delta_p^2 )^{-1/2}
# Short semi-axis from polar-dual identity:
#   S(theta) = hbar / L(theta)
#
# Author: Brian S. Mulloy
# ==============================================================================

library(here)
library(ggplot2)

dir_figures <- here("figures")
if (!dir.exists(dir_figures)) dir.create(dir_figures, recursive=TRUE)
file_output_svg <- file.path(dir_figures, "heisenberg_cells.svg")

# ------------------------------------------------------------------------------
# CELL DIMENSIONS
# ------------------------------------------------------------------------------

hbar    <- 1.0
Delta_q <- 2.0
Delta_p <- 1.5

#' Long semi-axis of the inscribed squeezed cell at orientation theta.
sq_L <- function(theta) {
  1 / sqrt( (cos(theta)/Delta_q)^2 + (sin(theta)/Delta_p)^2 )
}

#' Short semi-axis from polar-dual identity: S = hbar / L.
sq_S <- function(theta) {
  hbar / sq_L(theta)
}

# ------------------------------------------------------------------------------
# RECTANGLE PATH BUILDER
#
# Each cell is rendered as a rectangle whose width is twice the long
# semi-axis and whose height is twice the short semi-axis, rotated so
# that the long side lies along (cos theta, sin theta).
#
# The four corners (in the cell's own frame, before rotation) are
#   (+L, +S), (+L, -S), (-L, -S), (-L, +S)
# rotated by theta into (q, p).
# ------------------------------------------------------------------------------

#' Closed path (5 points: 4 corners + return-to-start) for a rectangle
#' with long semi-axis L along (cos theta, sin theta), short semi-axis S
#' perpendicular.
rect_path <- function(L, S, theta_orient, group_id) {
  ct <- cos(theta_orient)
  st <- sin(theta_orient)
  # Corners in body frame: long axis = x_body, short axis = y_body
  body_corners <- matrix(c( L,  S,
                            L, -S,
                            -L, -S,
                            -L,  S,
                            L,  S),    # close the path
                         ncol=2, byrow=TRUE)
  # Rotate body frame into (q, p): q = L_x cos - S_y sin, p = L_x sin + S_y cos
  q_vals <- body_corners[,1]*ct - body_corners[,2]*st
  p_vals <- body_corners[,1]*st + body_corners[,2]*ct
  data.frame(q=q_vals, p=p_vals, group=group_id)
}

# Extended cell: axis-aligned, "long axis" is along q, "short axis" is along p.
# Width 2*Delta_q, height 2*Delta_p.
df_A <- rect_path(L=Delta_q, S=Delta_p, theta_orient=0, group_id="A")

# Squeezed cells: each at its own orientation, with L and S from the cell
# functions above so the polar-dual identity L*S = hbar holds exactly.
df_a_p <- rect_path(L=sq_L(0),    S=sq_S(0),    theta_orient=0,    group_id="a_p")
df_a   <- rect_path(L=sq_L(pi/4), S=sq_S(pi/4), theta_orient=pi/4, group_id="a")
df_a_q <- rect_path(L=sq_L(pi/2), S=sq_S(pi/2), theta_orient=pi/2, group_id="a_q")

# ------------------------------------------------------------------------------
# DISPLAY WINDOW
#
# Just enough to contain the outer rectangle with a hairline of breathing
# room. OmniGraffle re-frames as needed; this is purely so the SVG
# viewBox isn't degenerate.
# ------------------------------------------------------------------------------

q_lim <- c(-1.05 * Delta_q,  1.05 * Delta_q)
p_lim <- c(-1.05 * Delta_p,  1.05 * Delta_p)

# ------------------------------------------------------------------------------
# PLOT (bare-bones)
#
# Four geom_path calls. theme_void() strips all decoration: no axes, no
# gridlines, no panel background, no panel border, no plot background, no
# margins. The result is an SVG containing only the four shape outlines.
# ------------------------------------------------------------------------------

p <- ggplot() +
  geom_path(data=df_A,   aes(x=q, y=p, group=group),
            color="black", linewidth=0.4) +
  geom_path(data=df_a_p, aes(x=q, y=p, group=group),
            color="black", linewidth=0.4) +
  geom_path(data=df_a_q, aes(x=q, y=p, group=group),
            color="black", linewidth=0.4) +
  geom_path(data=df_a,   aes(x=q, y=p, group=group),
            color="black", linewidth=0.4) +
  coord_fixed(xlim=q_lim, ylim=p_lim, expand=FALSE) +
  theme_void() +
  theme(plot.background  = element_rect(fill="transparent", color=NA),
        panel.background = element_rect(fill="transparent", color=NA),
        plot.margin      = margin(0, 0, 0, 0))

# ------------------------------------------------------------------------------
# SAVE (SVG, transparent background)
# ------------------------------------------------------------------------------

panel_in <- 3.5
ggsave(filename=file_output_svg, plot=p, device="svg",
       width=panel_in, height=panel_in,
       bg="transparent",
       limitsize=FALSE)

cat("Done.", file_output_svg, "\n")
