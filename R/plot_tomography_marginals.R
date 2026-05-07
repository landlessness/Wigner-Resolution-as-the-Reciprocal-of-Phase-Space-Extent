# ==============================================================================
# plot_tomography_marginals.R
# 1D marginal comparison from the new tilde_W reconstruction:
#   Position: rho_q(q) from tilde_W vs |psi(q)|^2
#   Momentum: rho_p(p) from tilde_W vs |psi-hat(p)|^2
#
# Each panel is an 80/20 stack: marginal comparison on top, residual
# bar chart below. The marginal panel is rendered via
# plot_semiclassical_resolution() from plot_tools.R so the model ribbon
# (gray60 fill + black line) is identical to the right column of
# plot_wigner.R / plot_cats.R. The exact density is added as an overlay
# styled like the prior Husimi overlay (gray85 fill alpha 0.5, gray30
# line). The residual bar chart is bespoke.
# ==============================================================================

library(here)
library(data.table)
library(ggplot2)
library(patchwork)

source(here("R", "plot_tools.R"))
source(here("R", "state_builder.R"))

DATA_RDS    <- here("data", "tomography_data.rds")
COMPUTE_R   <- here("R", "tomography_compute.R")
OUTPUT_PDF  <- here("figures", "tomography_marginals.pdf")

results <- load_tomography_cache(DATA_RDS, COMPUTE_R)

STATE_LABELS <- c(
  squeezed_vacuum = "Squeezed vacuum",
  harmonic_n1     = "Harmonic n=1",
  morse_n8        = "Morse n=8",
  double_well_n5  = "Double-well n=5",
  cat_2           = "2-cat",
  cat_3           = "3-cat",
  cat_4_square    = "4-cat square",
  cat_compass     = "Zurek compass"
)
STATE_ORDER <- names(STATE_LABELS)

latex_font <- "CMU Serif"

# Match the exact-density overlay style to the prior Husimi convention.
EXACT_OVERLAY_COLOR     <- "gray30"
EXACT_OVERLAY_LINEWIDTH <- 0.35
EXACT_OVERLAY_FILL      <- "gray85"
EXACT_OVERLAY_ALPHA     <- 0.5

# Match the residual bar style to the model ribbon color.
RESIDUAL_FILL  <- SEMICLASSICAL_RIBBON_FILL   # "gray60", same as model
ZERO_LINE_COL  <- "gray30"

# 80/20 height split.
PRIMARY_HEIGHT  <- 4
RESIDUAL_HEIGHT <- 1

# ------------------------------------------------------------------------------
# UTILITIES
# ------------------------------------------------------------------------------

l1_pct <- function(exact, approx, d_grid) {
  total <- sum(exact) * d_grid
  if (total <= 0) total <- 1
  100 * sum(abs(exact - approx)) * d_grid / total
}

# ------------------------------------------------------------------------------
# BUILD ONE STACKED PANEL (marginal + residual)
# ------------------------------------------------------------------------------

build_panel <- function(x_grid, exact, approx, custom_breaks,
                         x_label, y_label_main, title) {

  d_grid   <- diff(x_grid)[1]
  residual <- exact - approx
  l1       <- l1_pct(exact, approx, d_grid)

  # Y limits.
  y_max_main  <- max(c(exact, approx), na.rm=TRUE) * 1.1
  if (!is.finite(y_max_main) || y_max_main <= 0) y_max_main <- 1
  y_max_resid <- max(abs(residual), na.rm=TRUE) * 1.2
  if (!is.finite(y_max_resid) || y_max_resid <= 0) y_max_resid <- 1

  # Main panel via plot_semiclassical_resolution + exact overlay.
  dt_model <- data.table(q=x_grid, rho_sympl=approx)
  dt_exact <- data.frame(q=x_grid, rho=exact)
  exact_overlay <- list(
    list(
      data       = dt_exact,
      color      = EXACT_OVERLAY_COLOR,
      linewidth  = EXACT_OVERLAY_LINEWIDTH,
      fill       = EXACT_OVERLAY_FILL,
      fill_alpha = EXACT_OVERLAY_ALPHA
    )
  )

  p_main <- plot_semiclassical_resolution(
      dt_model, q_lim=range(x_grid), y_lim=y_max_main,
      custom_breaks=custom_breaks,
      label_format=function(v) sprintf("%.1f", v),
      base_font=latex_font,
      overlays=exact_overlay,
      y_label=y_label_main)

  # Add a small L1 annotation in the upper-left, plus the title.
  p_main <- p_main +
    annotate("text",
             x=min(x_grid) + 0.05*(max(x_grid)-min(x_grid)),
             y=0.95*y_max_main,
             label=sprintf("L1 = %.1f%%", l1),
             hjust=0, vjust=1, size=2.5, color="gray30",
             family=latex_font) +
    ggtitle(title) +
    theme(plot.title=element_text(size=10, hjust=0.5,
                                   family=latex_font,
                                   margin=margin(b=2)),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          plot.margin=margin(2, 2, 0, 2))

  # Residual panel: bar chart, exact - approx.
  dt_resid <- data.table(x=x_grid, y=residual)
  p_resid <- ggplot(dt_resid, aes(x=x, y=y)) +
    geom_hline(yintercept=0, color=ZERO_LINE_COL, linewidth=0.3) +
    geom_col(fill=RESIDUAL_FILL, width=d_grid) +
    coord_cartesian(xlim=range(x_grid),
                    ylim=c(-y_max_resid, y_max_resid), expand=FALSE) +
    scale_x_continuous(breaks=custom_breaks,
                       labels=function(v) sprintf("%.1f", v)) +
    scale_y_continuous(breaks=c(-y_max_resid, 0, y_max_resid),
                       labels=c(sprintf("%.2g", -y_max_resid),
                                "0",
                                sprintf("%.2g", y_max_resid))) +
    labs(x=x_label, y="resid") +
    theme_bw(base_family=latex_font) +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          axis.text=element_text(size=7),
          axis.title.x=element_text(size=8),
          axis.title.y=element_text(size=7),
          plot.margin=margin(0, 2, 2, 2))

  p_main / p_resid + plot_layout(heights=c(PRIMARY_HEIGHT, RESIDUAL_HEIGHT))
}

# ------------------------------------------------------------------------------
# BUILD ALL 16 PANELS
# ------------------------------------------------------------------------------

panels <- list()
for (s in STATE_ORDER) {
  r  <- results[[s]]
  qg <- r$q_grid_tomo
  pg <- r$p_grid_tomo

  panels[[paste0(s, "_q")]] <- build_panel(
    x_grid       = qg,
    exact        = r$exact_q,
    approx       = r$rho_q,
    custom_breaks= r$custom_breaks_q,
    x_label      = expression(italic(q)/italic(q)[0]),
    y_label_main = expression(italic(P)[italic(delta*q)](italic(q))),
    title        = paste0(STATE_LABELS[[s]], " (q)")
  )
  panels[[paste0(s, "_p")]] <- build_panel(
    x_grid       = pg,
    exact        = r$exact_p,
    approx       = r$rho_p,
    custom_breaks= r$custom_breaks_p,
    x_label      = expression(italic(p)/italic(p)[0]),
    y_label_main = expression(italic(P)[italic(delta*p)](italic(p))),
    title        = paste0(STATE_LABELS[[s]], " (p)")
  )
}

# Layout: 4 rows x 4 columns. Each row holds 2 states' (q, p) pair-pairs.
build_row <- function(s1, s2) {
  panels[[paste0(s1, "_q")]] |
    panels[[paste0(s1, "_p")]] |
    panels[[paste0(s2, "_q")]] |
    panels[[paste0(s2, "_p")]]
}

fig <- build_row("squeezed_vacuum", "harmonic_n1") /
       build_row("morse_n8",        "double_well_n5") /
       build_row("cat_2",           "cat_3") /
       build_row("cat_4_square",    "cat_compass")

dir.create(dirname(OUTPUT_PDF), showWarnings=FALSE, recursive=TRUE)
ggsave(OUTPUT_PDF, fig, width=14, height=8, units="in", device=cairo_pdf)
cat(sprintf("Wrote %s\n", OUTPUT_PDF))
