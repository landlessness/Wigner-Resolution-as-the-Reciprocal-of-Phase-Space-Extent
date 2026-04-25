# ==============================================================================
# plot_qho_berry_grid.R  —  Berry semiclassical grid, A_0 action units
# Three columns:
#   Left:   Berry 2D heatmap with quantum of action
#   Center: Raw Berry 1D density as ribbon — diverging at turning points
#   Right:  Berry convolved with squeezed kernel — finite everywhere
# ==============================================================================

library(here)
library(data.table)
library(ggplot2)
library(ggforce)
library(patchwork)

source(here::here("R", "symplectic_tools.R"))

# --- 1. Configuration ---
latex_font  <- "CMU Serif"
dir_figures <- here::here("figures")
if (!dir.exists(dir_figures)) dir.create(dir_figures, recursive = TRUE)
file_output_pdf <- file.path(dir_figures, "qho_berry_grid.pdf")

# --- 2. Data Selection ---
target_n_levels <- c(0, 1, 2, 3, 24)
dt_selected <- data.table(
  quantum_n = target_n_levels,
  action_A  = 2 * target_n_levels + 1
)

# --- 3. Berry Physics Functions ---

chord_area <- function(r, R) {
  R^2 * acos(pmin(r/R, 1)) - r * sqrt(pmax(R^2 - r^2, 0))
}

berry_wigner <- function(n, q, p) {
  R   <- sqrt(2*n + 1)
  r   <- sqrt(q^2 + p^2)
  val <- numeric(length(r))
  idx_in <- r < R
  if (any(idx_in)) {
    r_in  <- r[idx_in]
    Area  <- chord_area(r_in, R)
    denom <- pmax((R^2 - r_in^2)^0.25, 1e-6)
    val[idx_in] <- (1 / (pi * denom)) * cos(Area - pi/4)
  }
  idx_out <- r >= R
  if (any(idx_out)) {
    r_out    <- r[idx_out]
    Area_out <- r_out * sqrt(r_out^2 - R^2) - R^2 * acosh(r_out/R)
    denom    <- pmax((r_out^2 - R^2)^0.25, 1e-6)
    val[idx_out] <- (1 / (2*pi * denom)) * exp(-Area_out)
  }
  return(val)
}

berry_density_1d <- function(n, q_seq) {
  R   <- sqrt(2*n + 1)
  val <- rep(0, length(q_seq))
  idx_in <- abs(q_seq) < R
  if (any(idx_in)) {
    q_in  <- q_seq[idx_in]
    Area  <- chord_area(abs(q_in), R)
    denom <- pmax(sqrt(R^2 - q_in^2), 1e-6)
    val[idx_in] <- (2 / (pi * denom)) * cos(Area/2 - pi/4)^2
  }
  idx_out <- abs(q_seq) >= R
  if (any(idx_out)) {
    q_out    <- q_seq[idx_out]
    Area_out <- abs(q_out) * sqrt(q_out^2 - R^2) - R^2 * acosh(abs(q_out)/R)
    denom    <- pmax(sqrt(q_out^2 - R^2), 1e-6)
    val[idx_out] <- (1 / (2*pi * denom)) * exp(-2*Area_out)
  }
  return(val)
}

convolve_berry_1d <- function(berry_raw, q_seq, delta_q) {
  dq    <- diff(q_seq)[1]
  n_pts <- length(q_seq)
  G_1d  <- exp(-q_seq^2 / delta_q^2)
  G_1d  <- G_1d / (sum(G_1d) * dq)
  P_conv <- convolve(berry_raw, rev(G_1d), type="open")
  half   <- floor(n_pts / 2)
  start  <- half + 1
  P_out  <- P_conv[start:(start + n_pts - 1)]
  P_out  <- pmax(P_out, 0)
  P_out  / (sum(P_out) * dq)
}

# --- 4. Plot Function ---
plot_berry_grid <- function(dt_meta, base_font = "") {
  plot_list <- list()
  num_rows  <- nrow(dt_meta)

  ax_x     <- expression(italic(q)/italic(q)[0])
  ax_y     <- expression(italic(p)/italic(p)[0])
  ax_y_raw <- expression(italic(P)[Berry](italic(q)))
  ax_y_den <- expression(italic(P)[delta*italic(q)](italic(q)))

  for (i in seq_len(num_rows)) {
    n_val <- dt_meta$quantum_n[i]
    alpha <- dt_meta$action_A[i]

    rs <- qho_covariance(n_val)

    cat(sprintf(
      "\nn=%d | A/A_0=%.0f | Delta_q=%.4f | delta_q=%.4f | RS: %s | SP: %s\n",
      n_val, rs$A_over_A0, rs$Delta_q, rs$delta_q,
      ifelse(rs$rs_satisfied, "OK", "FAIL"),
      ifelse(rs$sp_satisfied, "OK", "FAIL")
    ))

    Delta_q <- rs$Delta_q
    Delta_p <- rs$Delta_p
    delta_q <- rs$delta_q
    delta_p <- rs$delta_p

    ell_lim       <- Delta_q * 1.25
    plot_lim      <- max(Delta_q * 1.25, Delta_q + 2.5)
    break_val     <- round(Delta_q, 1)
    break_val     <- min(break_val, floor(ell_lim * 10) / 10)
    custom_breaks <- c(-break_val, 0, break_val)
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

    # Grids
    grid_res  <- max(400, 15 * alpha)
    q_display <- seq(-plot_lim, plot_lim, length.out = grid_res)
    dq_disp   <- diff(q_display)[1]

    # Berry 2D heatmap
    q_ell     <- seq(-ell_lim, ell_lim, length.out = 400)
    p_ell_seq <- seq(-ell_lim, ell_lim, length.out = 400)
    dt_w2d    <- as.data.table(expand.grid(q = q_ell, p = p_ell_seq))
    dt_w2d[, w := berry_wigner(n_val, q, p)]
    dt_w2d[, w_plot := sign(w) * abs(w)^0.4]
    max_w <- max(abs(dt_w2d$w_plot), na.rm = TRUE)
    if (max_w > 0) dt_w2d[, w_plot := w_plot / max_w]

    # Berry 1D raw and convolved
    berry_raw    <- berry_density_1d(n_val, q_display)
    P_berry_conv <- convolve_berry_1d(berry_raw, q_display, delta_q)

    cat(sprintf("  Berry conv peak: %.4f at q=%.4f\n",
                max(P_berry_conv), q_display[which.max(P_berry_conv)]))

    # Scale berry_raw for display to match convolved density interior amplitude
    interior_mask <- abs(q_display) < Delta_q * 0.6 &
      berry_raw < quantile(berry_raw[abs(q_display) < Delta_q], 0.7, na.rm=TRUE)
    if (any(interior_mask) && any(P_berry_conv[interior_mask] > 0)) {
      scale_factor <- mean(P_berry_conv[interior_mask], na.rm=TRUE) /
        max(mean(berry_raw[interior_mask], na.rm=TRUE), 1e-10)
    } else {
      scale_factor <- 1 / max(berry_raw[abs(q_display) < Delta_q * 0.9], na.rm=TRUE)
    }
    berry_raw_display <- berry_raw * scale_factor

    density_peak <- max(P_berry_conv, na.rm=TRUE)
    y_lim_den    <- density_peak * 1.4
    clip_y       <- density_peak * 1.2

    dt_density <- data.table(
      q                 = q_display,
      density           = P_berry_conv,
      berry_raw_display = berry_raw_display
    )

    dt_center <- dt_density[abs(q) < Delta_q & berry_raw_display <= clip_y]

    # Ellipse overlays
    df_circles <- data.frame(x0=0, y0=0, r_A=Delta_q)
    df_cigars  <- data.frame(
      x0=0, y0=0,
      aq_a=delta_q, aq_b=Delta_p,
      ap_a=Delta_q, ap_b=delta_p
    )

    # Column 1: Berry 2D heatmap with symplectic blobs
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

    # Column 2: Raw Berry 1D as ribbon — shows divergence at turning points
    # Line spans exactly from -Delta_q to +Delta_q with compact support
    dt_center   <- dt_density[q >= -Delta_q & q <= Delta_q &
                                berry_raw_display <= clip_y]

    dt_interior <- copy(dt_density[q >= -Delta_q & q <= Delta_q])
    dt_interior[, fill_height := pmin(berry_raw_display, y_lim_den * 0.88)]

    p_raw <- ggplot() +
      # Ribbon fills from zero to the Berry density — uncapped, reaches infinity symbols
      geom_area(data=dt_interior, aes(x=q, y=fill_height),
                fill="gray85", color=NA) +
      # Line runs full interior with arrows at both turning points
      geom_path(data=dt_center[order(q)], aes(x=q, y=berry_raw_display),
                color="black", linewidth=0.4,
                arrow=arrow(length=unit(0.12,"cm"), ends="both", type="closed")) +
      annotate("text", x=-Delta_q, y=y_lim_den*0.92,
               label="infinity", parse=TRUE, color="gray30", size=4.5) +
      annotate("text", x= Delta_q, y=y_lim_den*0.92,
               label="infinity", parse=TRUE, color="gray30", size=4.5) +
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
      labs(x=ax_x, y=ax_y_raw)

    # Column 3: Convolved Berry — finite everywhere
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
      p_ell   <- p_ell   + labs(title="Phase Space") +
        theme(plot.title=element_text(size=11, hjust=0.5))
      p_raw   <- p_raw   + labs(title="Berry Density") +
        theme(plot.title=element_text(size=11, hjust=0.5))
      p_den   <- p_den   + labs(title="Symplectic Density") +
        theme(plot.title=element_text(size=11, hjust=0.5))
    }

    # X-axis titles on bottom row only
    if (i != num_rows) {
      p_ell <- p_ell + theme(axis.title.x=element_blank())
      p_raw <- p_raw + theme(axis.title.x=element_blank())
      p_den <- p_den + theme(axis.title.x=element_blank())
    }

    # Y-axis titles on middle row only
    if (i != 3) {
      p_ell <- p_ell + theme(axis.title.y=element_blank())
      p_raw <- p_raw + theme(axis.title.y=element_blank())
      p_den <- p_den + theme(axis.title.y=element_blank())
    }

    plot_list <- c(plot_list, list(p_label, p_ell, p_raw, p_den))
  }

  wrap_plots(plot_list, ncol=4, widths=c(0.25, 1, 1, 1))
}

# --- 5. Execution ---
cat("Computing Berry semiclassical grid with squeezed convolution...\n")

p_final <- plot_berry_grid(dt_meta=dt_selected, base_font=latex_font)
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
