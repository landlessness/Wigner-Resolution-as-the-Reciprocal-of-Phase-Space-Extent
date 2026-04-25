# ==============================================================================
# plot_qho_grid.R  —  QHO symplectic density grid, A_0 action units
# ==============================================================================

library(data.table)
library(ggplot2)
library(ggforce)
library(patchwork)
library(gsl)

source(here::here("R", "symplectic_tools.R"))

# --- 1. Configuration ---
latex_font  <- "CMU Serif"
dir_figures <- "./figures"
if (!dir.exists(dir_figures)) dir.create(dir_figures, recursive = TRUE)
file_output_pdf <- file.path(dir_figures, "qho_grid.pdf")

# --- 2. Data Selection ---
# Action leads. Quantum number is parenthetical.
target_n_levels <- c(0, 1, 2, 3, 24)
dt_selected <- data.table(
  quantum_n = target_n_levels,
  action_A  = 2 * target_n_levels + 1
)

# --- 3. Plot Function ---
plot_qho_grid <- function(dt_meta, base_font = "") {
  plot_list <- list()
  num_rows  <- nrow(dt_meta)

  ax_x     <- expression(italic(q)/italic(q)[0])
  ax_y     <- expression(italic(p)/italic(p)[0])
  ax_y_cdf <- expression("Pr(" * italic(Q) <= italic(q) * ")")
  ax_y_den <- expression(italic(P)[delta*italic(q)](italic(q)))

  for (i in seq_len(num_rows)) {
    n_val <- dt_meta$quantum_n[i]
    alpha <- dt_meta$action_A[i]  # A/A_0 = 2n+1 for QHO

    # Robertson-Schrödinger covariance for QHO eigenstate n
    # sigma_qq = sigma_pp = alpha/2 in natural units (hbar=1)
    rs <- qho_covariance(n_val)

    cat(sprintf("\nn=%d | A/A_0=%.0f | Delta_q=%.4f | delta_q=%.4f | RS satisfied: %s\n",
                n_val, rs$A_over_A0, rs$Delta_q, rs$delta_q,
                ifelse(rs$rs_satisfied, "YES", "NO")))

    Delta_q <- rs$Delta_q
    Delta_p <- rs$Delta_p
    delta_q <- rs$delta_q
    delta_p <- rs$delta_p

    ell_lim       <- Delta_q * 1.15
    plot_lim      <- max(Delta_q * 1.15, Delta_q + 2.5)
    custom_breaks <- c(-round(Delta_q, 1), 0, round(Delta_q, 1))
    label_format  <- function(x) sprintf("%.1f", x)

    # Row label
    p_label <- ggplot() + theme_void() +
      coord_cartesian(xlim=c(-1,1), ylim=c(0,1), clip="off") +
      annotate("text", x=0, y=0.5,
               label  = sprintf("%.0f*A[0]", alpha),
               parse  = TRUE,
               family = base_font,
               size   = 4.5,
               hjust  = 0.5)

    # Display grid (ell_lim) and computation grid (plot_lim)
    q_display <- seq(-plot_lim, plot_lim, length.out = max(400, 15 * alpha))
    dq_disp   <- diff(q_display)[1]

    # 2D Wigner on display grid for phase-space panel
    q_ell     <- seq(-ell_lim, ell_lim, length.out = 400)
    p_ell_seq <- seq(-ell_lim, ell_lim, length.out = 400)
    dt_w2d    <- as.data.table(expand.grid(q = q_ell, p = p_ell_seq))
    dt_w2d[, w := qho_wigner(n_val, q, p)]
    dt_w2d[, w_plot := sign(w) * abs(w)^0.4]
    max_w <- max(abs(dt_w2d$w_plot), na.rm = TRUE)
    if (max_w > 0) dt_w2d[, w_plot := w_plot / max_w]

    # Symplectic density via shared pipeline
    result <- compute_symplectic_density(
      n          = n_val,
      wigner_fn  = qho_wigner,
      kernel_fn  = squeezed_kernel_q,
      rs         = rs,
      q_display  = q_display
    )

    cat(sprintf("  Wigner norm: %.6f | max_negative: %.2e | tolerance: %.2e\n",
                result$w_norm, result$max_negative, result$tolerance))

    P_q <- result$P_q

    # Verify display-grid integral (diagnostic only, no renormalization)
    disp_integral <- sum(P_q) * dq_disp
    cat(sprintf("  Display grid integral (informational): %.6f\n", disp_integral))

    dt_density <- data.table(
      q       = q_display,
      density = P_q,
      cdf     = cumsum(P_q) * dq_disp
    )

    # Ellipse overlays — derived from RS, not hardcoded
    df_circles <- data.frame(x0=0, y0=0, r_A=Delta_q)
    df_cigars  <- data.frame(
      x0=0, y0=0,
      aq_a=delta_q, aq_b=Delta_p,
      ap_a=Delta_q, ap_b=delta_p
    )

    # Column 1: Phase-space Wigner distribution with symplectic blobs
    p_ell <- ggplot(dt_w2d, aes(x=q, y=p)) +
      geom_circle(data=df_circles, aes(x0=x0, y0=y0, r=r_A),
                  inherit.aes=FALSE, fill="white", color=NA) +
      geom_raster(aes(fill=w_plot), interpolate=TRUE) +
      scale_fill_gradient2(low="gray10", mid="white", high="gray40",
                           midpoint=0, limits=c(-1,1), guide="none") +
      geom_circle(data=df_circles, aes(x0=x0, y0=y0, r=r_A),
                  inherit.aes=FALSE, color="black", linewidth=0.3) +
      geom_ellipse(data=df_cigars, aes(x0=x0, y0=y0, a=aq_a, b=aq_b, angle=0),
                   inherit.aes=FALSE, color="black", linewidth=0.5) +
      geom_ellipse(data=df_cigars, aes(x0=x0, y0=y0, a=ap_a, b=ap_b, angle=0),
                   inherit.aes=FALSE, color="black", linewidth=0.5) +
      coord_fixed(xlim=c(-ell_lim, ell_lim),
                  ylim=c(-ell_lim, ell_lim), expand=FALSE) +
      scale_x_continuous(breaks=custom_breaks, labels=label_format) +
      scale_y_continuous(breaks=custom_breaks, labels=label_format) +
      theme_bw(base_family=base_font) +
      theme(panel.grid.minor=element_blank(),
            panel.background=element_rect(fill="white"),
            axis.text=element_text(size=8),
            plot.margin=margin(2,4,2,4)) +
      labs(x=ax_x, y=ax_y)

    # Column 2: Quantization map (CDF of symplectic density)
    p_stair <- ggplot(dt_density, aes(x=q, y=cdf)) +
      geom_line(color="black", linewidth=0.8) +
      coord_cartesian(xlim=c(-ell_lim, ell_lim), ylim=c(0,1), expand=FALSE) +
      scale_x_continuous(breaks=custom_breaks, labels=label_format) +
      theme_bw(base_family=base_font) +
      theme(panel.grid.minor=element_blank(),
            axis.text=element_text(size=8),
            aspect.ratio=1,
            plot.margin=margin(2,4,2,4)) +
      labs(x=ax_x, y=ax_y_cdf)

    # Column 3: Symplectic density P_delta_q(q)
    density_peak <- max(dt_density$density, na.rm=TRUE)
    y_lim_den    <- density_peak * 1.15

    p_den <- ggplot() +
      geom_ribbon(data=dt_density, aes(x=q, ymin=0, ymax=density),
                  fill="gray85", color=NA) +
      geom_path(data=dt_density, aes(x=q, y=density),
                color="black", linewidth=0.4) +
      coord_cartesian(xlim=c(-ell_lim, ell_lim),
                      ylim=c(0, y_lim_den), expand=FALSE) +
      scale_x_continuous(breaks=custom_breaks, labels=label_format) +
      theme_bw(base_family=base_font) +
      theme(panel.grid.minor=element_blank(),
            axis.text=element_text(size=8),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            aspect.ratio=1,
            plot.margin=margin(2,4,2,4)) +
      labs(x=ax_x, y=ax_y_den)

    # Headers on first row only
    if (i == 1) {
      p_label <- p_label + labs(title=" ") +
        theme(plot.title=element_text(size=11, hjust=0.5))
      p_ell   <- p_ell   + labs(title="Quantum of Action") +
        theme(plot.title=element_text(size=11, hjust=0.5))
      p_stair <- p_stair + labs(title="Quantization Map") +
        theme(plot.title=element_text(size=11, hjust=0.5))
      p_den   <- p_den   + labs(title="Symplectic Density") +
        theme(plot.title=element_text(size=11, hjust=0.5))
    }

    # X-axis titles on bottom row only
    if (i != num_rows) {
      p_ell   <- p_ell   + theme(axis.title.x=element_blank())
      p_stair <- p_stair + theme(axis.title.x=element_blank())
      p_den   <- p_den   + theme(axis.title.x=element_blank())
    }

    # Y-axis titles on middle row only
    if (i != 3) {
      p_ell   <- p_ell   + theme(axis.title.y=element_blank())
      p_stair <- p_stair + theme(axis.title.y=element_blank())
      p_den   <- p_den   + theme(axis.title.y=element_blank())
    }

    plot_list <- c(plot_list, list(p_label, p_ell, p_stair, p_den))
  }

  wrap_plots(plot_list, ncol=4, widths=c(0.25, 1, 1, 1))
}

# --- 4. Execution ---
cat("Computing QHO symplectic densities via Robertson-Schrödinger pipeline...\n")

p_final <- plot_qho_grid(dt_meta=dt_selected, base_font=latex_font)
p_final <- p_final + theme(plot.margin=margin(10,10,10,10))

fig_width  <- 7.0
fig_height <- nrow(dt_selected) * 1.8 + 0.5

ggsave(
  filename  = file_output_pdf,
  plot      = p_final,
  device    = cairo_pdf,
  width     = fig_width,
  height    = fig_height,
  limitsize = FALSE
)

cat("Done.", file_output_pdf, "\n")
