# ==============================================================================
# plot_marginals.R
# Figure: position and momentum marginals of P_{delta q} and P_{delta p}
# vs |psi(q)|^2 and |psi-hat(p)|^2 across eight non-Gaussian states.
#
# Layout: 4 rows x 4 columns of state-pairs. Each "state-pair" is two
# adjacent panels: q-marginal on the left, p-marginal on the right.
# Eight states fill 8 state-pairs, arranged as 4 across x 2 down, but
# rendered as a 16-panel grid.
#
# Within each panel: ~80% height for the marginal comparison (exact line,
# convolved fill), ~20% height below for the residual (approx - exact)
# on its own small axis.
#
# Reads:  data/marginals_data.rds  (produced by table_marginals.R)
# Writes: figures/marginals.pdf
#
# Author: Brian S. Mulloy
# ==============================================================================

library(here)
library(ggplot2)
library(patchwork)
library(data.table)

# ------------------------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------------------------

OUTPUT_PDF       <- here("figures", "marginals.pdf")
DATA_RDS         <- here("data", "marginals_data.rds")

# ------------------------------------------------------------------------------
# ENSURE DATA EXISTS
#
# If marginals_data.rds is missing, run the driver first. Delete the file
# to force regeneration on the next run.
# ------------------------------------------------------------------------------

if (!file.exists(DATA_RDS)) {
  cat(sprintf("%s not found; running table_marginals.R to generate it...\n",
              DATA_RDS))
  source(here("R", "table_marginals.R"))
}

PLOT_WIDTH_IN    <- 14         # full page width in inches (PRL figure*)
PLOT_HEIGHT_IN   <- 8          # top half of page

PRIMARY_HEIGHT   <- 4          # relative height of marginal panel
RESIDUAL_HEIGHT  <- 1          # relative height of residual panel  (4:1 = 80/20)

# Style constants matching plot_tools.R (existing pipeline).
# Model rho_{delta q}: gray60 ribbon + black line, just like the right-column
# panels in plot_wigner.R / plot_cats.R.
MODEL_FILL       <- "gray60"
MODEL_LINE       <- "black"
MODEL_LINEWIDTH  <- 0.4

# Exact |psi|^2 / |psi-hat|^2: rendered as the Husimi overlay was rendered
# (gray85 fill alpha 0.5, gray30 line, linewidth 0.35).
EXACT_FILL       <- "gray85"
EXACT_FILL_ALPHA <- 0.5
EXACT_LINE       <- "gray30"
EXACT_LINEWIDTH  <- 0.35

# Residual bars: same gray as the model ribbon for visual coherence.
RESIDUAL_FILL    <- "gray60"
ZERO_LINE_COLOR  <- "gray30"

BASE_FONT_SIZE   <- 8

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

# Order matches table_marginals.R driver loop.
STATE_ORDER <- c("squeezed_vacuum", "harmonic_n1", "morse_n8", "double_well_n5",
                 "cat_2", "cat_3", "cat_4_square", "cat_compass")

# ------------------------------------------------------------------------------
# LOAD
# ------------------------------------------------------------------------------

cat(sprintf("Reading %s\n", DATA_RDS))
results <- readRDS(DATA_RDS)

# ------------------------------------------------------------------------------
# BUILD ONE PANEL (marginal + residual stack)
#
# Color scheme matches the existing pipeline:
#   - Model rho_{delta q}: gray60 ribbon + black line (as in the right
#     column of plot_wigner.R / plot_cats.R).
#   - Exact |psi|^2: rendered as the prior Husimi overlay was rendered
#     (gray85 fill alpha 0.5, gray30 line). The exact density now plays
#     the role the prior-art overlay used to play.
#   - Residual: bar chart of (exact - model). Positive bars mark places
#     where smoothing reduced the value (peaks); negative bars mark
#     where smoothing redistributed mass into (valleys, tails).
#
# Drawing order on the marginal panel:
#   1. model ribbon (gray60)
#   2. exact ribbon (gray85 alpha 0.5) — alpha-blends on top of model
#   3. model line (black) — sits above all fills
#   4. exact line (gray30) — sits above all fills
# ------------------------------------------------------------------------------

build_panel <- function(x_grid, exact, approx, l1_pct,
                         x_label, y_label, title) {
  d_grid <- diff(x_grid)[1]
  # Sign convention: exact - model. Positive = model under-estimates
  # (typically at peaks, since the model is a smoothed version).
  residual <- exact - approx
  dt_model    <- data.table(x=x_grid, y=approx)
  dt_exact    <- data.table(x=x_grid, y=exact)
  dt_residual <- data.table(x=x_grid, y=residual)

  y_max_main  <- max(c(exact, approx), na.rm=TRUE) * 1.1
  if (!is.finite(y_max_main) || y_max_main <= 0) y_max_main <- 1
  y_max_resid <- max(abs(residual), na.rm=TRUE) * 1.2
  if (!is.finite(y_max_resid) || y_max_resid <= 0) y_max_resid <- 1

  # y breaks: 0 and the peak, just like plot_semiclassical_resolution.
  rho_peak  <- max(approx, na.rm=TRUE)
  y_breaks  <- c(0, rho_peak)
  y_labels  <- c("0", sprintf("%.2g", rho_peak))

  # Primary panel: model ribbon + exact ribbon overlay + lines on top.
  p_main <- ggplot() +
    # 1. Model ribbon (gray60, the existing P_{delta q} style).
    geom_ribbon(data=dt_model,
                aes(x=x, ymin=0, ymax=pmax(y, 0)),
                fill=MODEL_FILL, color=NA) +
    # 2. Exact ribbon overlay (Husimi-style: gray85, alpha 0.5).
    geom_ribbon(data=dt_exact,
                aes(x=x, ymin=0, ymax=pmax(y, 0)),
                fill=EXACT_FILL, color=NA, alpha=EXACT_FILL_ALPHA) +
    # 3. Model line (black).
    geom_path(data=dt_model, aes(x=x, y=y),
              color=MODEL_LINE, linewidth=MODEL_LINEWIDTH) +
    # 4. Exact line (gray30, Husimi-overlay style).
    geom_path(data=dt_exact, aes(x=x, y=y),
              color=EXACT_LINE, linewidth=EXACT_LINEWIDTH) +
    annotate("text",
             x=min(x_grid) + 0.05*(max(x_grid)-min(x_grid)),
             y=0.95*y_max_main,
             label=sprintf("L1 = %.2f%%", l1_pct),
             hjust=0, vjust=1,
             size=BASE_FONT_SIZE*0.32, color="gray30") +
    coord_cartesian(xlim=range(x_grid), ylim=c(0, y_max_main), expand=FALSE) +
    scale_y_continuous(breaks=y_breaks, labels=y_labels) +
    labs(title=title, y=y_label) +
    theme_bw(base_size=BASE_FONT_SIZE) +
    theme(
      plot.title       = element_text(size=BASE_FONT_SIZE,
                                       face="plain", hjust=0.5,
                                       margin=margin(b=2)),
      axis.title.x     = element_blank(),
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.title.y     = element_text(size=BASE_FONT_SIZE-1),
      axis.text.y      = element_text(size=BASE_FONT_SIZE-2),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.margin      = margin(2, 2, 0, 2)
    )

  # Residual panel: bar chart, exact - model.
  # Bars width set to grid spacing so adjacent bars touch.
  p_resid <- ggplot(dt_residual, aes(x=x, y=y)) +
    geom_hline(yintercept=0, color=ZERO_LINE_COLOR, linewidth=0.3) +
    geom_col(fill=RESIDUAL_FILL, width=d_grid) +
    coord_cartesian(xlim=range(x_grid),
                    ylim=c(-y_max_resid, y_max_resid), expand=FALSE) +
    scale_y_continuous(breaks=c(-y_max_resid, 0, y_max_resid),
                       labels=c(sprintf("%.2g", -y_max_resid),
                                "0",
                                sprintf("%.2g", y_max_resid))) +
    labs(x=x_label, y="resid") +
    theme_bw(base_size=BASE_FONT_SIZE) +
    theme(
      axis.title.x     = element_text(size=BASE_FONT_SIZE-1),
      axis.text.x      = element_text(size=BASE_FONT_SIZE-2),
      axis.title.y     = element_text(size=BASE_FONT_SIZE-2),
      axis.text.y      = element_text(size=BASE_FONT_SIZE-3),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.margin      = margin(0, 2, 2, 2)
    )

  # Stack with 80/20 height split.
  p_main / p_resid + plot_layout(heights=c(PRIMARY_HEIGHT, RESIDUAL_HEIGHT))
}

# ------------------------------------------------------------------------------
# BUILD ALL 16 PANELS (2 per state, q and p)
# ------------------------------------------------------------------------------

cat("Building panels...\n")

panels <- list()
for (state_name in STATE_ORDER) {
  r <- results[[state_name]]
  label <- STATE_LABELS[[state_name]]

  panel_q <- build_panel(
    x_grid   = r$q_display,
    exact    = r$exact_q,
    approx   = r$pdq_q_marg,
    l1_pct   = 100 * r$res_q$l1,
    x_label  = expression(italic(q)/italic(q)[0]),
    y_label  = expression(group("|", italic(psi)(italic(q)), "|")^2),
    title    = paste0(label, " (q)")
  )
  panel_p <- build_panel(
    x_grid   = r$p_display,
    exact    = r$exact_p,
    approx   = r$pdp_p_marg,
    l1_pct   = 100 * r$res_p$l1,
    x_label  = expression(italic(p)/italic(p)[0]),
    y_label  = expression(group("|", italic(hat(psi))(italic(p)), "|")^2),
    title    = paste0(label, " (p)")
  )
  panels[[paste0(state_name, "_q")]] <- panel_q
  panels[[paste0(state_name, "_p")]] <- panel_p
}

# ------------------------------------------------------------------------------
# ASSEMBLE INTO 4 ROWS x 4 COLUMNS
#
# Row 1: states 1-2  (each occupies (q,p) = 2 columns)
# Row 2: states 3-4
# Row 3: states 5-6
# Row 4: states 7-8
#
# That gives a 4x4 grid of (q,p) pairs. Eight states fill it exactly.
# Wait: 8 states / 2 states per row = 4 rows, 4 columns. Correct.
# ------------------------------------------------------------------------------

cat("Assembling figure...\n")

build_row <- function(s1, s2) {
  panels[[paste0(s1, "_q")]] |
    panels[[paste0(s1, "_p")]] |
    panels[[paste0(s2, "_q")]] |
    panels[[paste0(s2, "_p")]]
}

row1 <- build_row("squeezed_vacuum", "harmonic_n1")
row2 <- build_row("morse_n8",        "double_well_n5")
row3 <- build_row("cat_2",           "cat_3")
row4 <- build_row("cat_4_square",    "cat_compass")

figure <- row1 / row2 / row3 / row4

# ------------------------------------------------------------------------------
# SAVE
# ------------------------------------------------------------------------------

dir.create(dirname(OUTPUT_PDF), showWarnings=FALSE, recursive=TRUE)
ggsave(OUTPUT_PDF, figure,
       width=PLOT_WIDTH_IN, height=PLOT_HEIGHT_IN, units="in",
       device=cairo_pdf)

cat(sprintf("Wrote %s\n", OUTPUT_PDF))
