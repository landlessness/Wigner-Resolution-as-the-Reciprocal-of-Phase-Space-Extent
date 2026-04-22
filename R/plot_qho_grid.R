# ==============================================================================
# plot_qho_grid.R
# ==============================================================================

library(data.table)
library(ggplot2)
library(ggforce)
library(patchwork)
library(stats)

# --- 1. Configuration ---
latex_font <- "CMU Serif"
dir_figures <- "./figures"
if (!dir.exists(dir_figures)) dir.create(dir_figures, recursive = TRUE)

# Output filename
file_output_pdf  <- file.path(dir_figures, "qho_grid.pdf")

# --- 2. Data Selection ---
# Matches textbook progression (0,1,2,3) and jumps to a clear semiclassical state (24)
target_n_levels <- c(0, 1, 2, 3, 24)

dt_selected <- data.table(
  quantum_n = target_n_levels,
  action_A = 2 * target_n_levels + 1
)

#' Render the Wigner-Husimi Kinematic Grid (Pure QHO Chord Geometry)
#'
#' @param dt_meta data.table containing quantum_n and action_A.
#' @param style String, either "manuscript" or "readme".
#' @param base_font String for the font family (e.g., "CMU Serif").
#' @param overarching_title Optional title for the readme style.
plot_qho_grid <- function(dt_meta, style = c("manuscript", "readme"),
                          base_font = "", overarching_title = NULL) {
  style <- match.arg(style)
  plot_list <- list()
  num_rows <- nrow(dt_meta)
  is_ms <- style == "manuscript"

  # Axis Labels
  ax_x_ell <- if(is_ms) expression("Position, " * italic(q)/italic(q)[0]) else "Position, q/q_0"
  ax_y_ell <- if(is_ms) expression("Momentum, " * italic(p)/italic(p)[0]) else "Momentum, p/p_0"
  ax_x_stair <- if(is_ms) expression("Position, " * italic(q)/italic(q)[0]) else "Position, q/q_0"
  ax_y_stair <- "Cumulative Probability"
  ax_x_den <- if(is_ms) expression("Position, " * italic(q)/italic(q)[0]) else "Position, q/q_0"
  ax_y_den <- "Probability Density P(q)"

  # --- 1. PURE GEOMETRIC FUNCTIONS (No Schrödinger Wave Equations) ---

  # Calculate the geometric area of the circular segment cut by the chord
  chord_area <- function(r, R) {
    return(R^2 * acos(r/R) - r * sqrt(R^2 - r^2))
  }

  # QHO Geometric 2D Wigner Landscape
  qho_wigner <- function(n, q, p) {
    R <- sqrt(2*n + 1)
    r <- sqrt(q^2 + p^2)
    val <- numeric(length(r))

    idx_in <- r < R
    if (any(idx_in)) {
      r_in <- r[idx_in]
      Area <- chord_area(r_in, R)
      # 1e-6 cap to let the singularity soar safely
      denom <- pmax((R^2 - r_in^2)^0.25, 1e-6)
      val[idx_in] <- (1 / (pi * denom)) * cos(Area - pi/4)
    }

    # Evanescent geometric tail
    idx_out <- r >= R
    if (any(idx_out)) {
      r_out <- r[idx_out]
      Area_out <- r_out * sqrt(r_out^2 - R^2) - R^2 * acosh(r_out/R)
      denom <- pmax((r_out^2 - R^2)^0.25, 1e-6)
      val[idx_out] <- (1 / (2 * pi * denom)) * exp(-Area_out)
    }
    return(val)
  }

  # Semiclassical WKB Spatial Density (1D projection of the chord geometry)
  qho_wkb_density <- function(n, q) {
    R <- sqrt(2*n + 1)
    val <- numeric(length(q))

    idx_in <- abs(q) < R
    if (any(idx_in)) {
      q_in <- q[idx_in]
      Area <- chord_area(abs(q_in), R)
      # 1e-6 cap to preserve the true mathematical singularity before integration
      class_prob <- 1 / (pi * pmax(sqrt(R^2 - q_in^2), 1e-6))
      val[idx_in] <- 2 * class_prob * cos((0.5 * Area) - pi/4)^2
    }

    idx_out <- abs(q) >= R
    if (any(idx_out)) {
      q_out <- q[idx_out]
      Area_out <- abs(q_out) * sqrt(q_out^2 - R^2) - R^2 * acosh(abs(q_out)/R)
      class_prob <- 1 / (2 * pi * pmax(sqrt(q_out^2 - R^2), 1e-6))
      val[idx_out] <- class_prob * exp(-2 * Area_out)
    }
    return(val)
  }

  # --- 2. MAIN RENDERING LOOP ---
  for (i in seq_len(num_rows)) {
    current_row <- dt_meta[i]
    n_val   <- current_row$quantum_n
    a_ratio <- current_row$action_A

    Delta_q <- sqrt(a_ratio)
    Delta_p <- sqrt(a_ratio)
    delta_q <- 1 / Delta_p
    delta_p <- 1 / Delta_q

    plot_lim_x <- Delta_q * 1.15
    custom_breaks <- c(-round(Delta_q, 1), 0, round(Delta_q, 1))
    label_format <- function(x) sprintf("%.1f", x)

    # Row Label: Centered and properly parsed to prevent clipping
    p_label <- ggplot() + theme_void() + coord_cartesian(xlim = c(-1, 1), ylim = c(0, 1), clip = "off")
    row_label_str <- sprintf("%.1f*A[0]", a_ratio)
    p_label <- p_label + annotate("text", x=0, y=0.5, label=row_label_str, parse=TRUE, family=base_font, size=4.5, hjust=0.5)

    # Generate Grid Data
    grid_res <- 1200
    q_seq <- seq(-plot_lim_x, plot_lim_x, length.out = grid_res)
    p_seq <- seq(-plot_lim_x, plot_lim_x, length.out = grid_res)

    dt_wigner2d <- as.data.table(expand.grid(q = q_seq, p = p_seq))
    dt_wigner2d[, w := qho_wigner(n_val, q, p)]

    # --- CONTRAST BOOST ---
    dt_wigner2d[, w_plot := sign(w) * (abs(w)^0.4)]
    max_w <- max(abs(dt_wigner2d$w_plot), na.rm = TRUE)
    if(max_w > 0) dt_wigner2d[, w_plot := w_plot / max_w]

    # --- 3. THE KINEMATIC SCAN (Direct Spatial Integration) ---
    raw_geometric_density <- qho_wkb_density(n_val, q_seq)
    step_size <- diff(q_seq)[1]
    sigma_val <- delta_q / sqrt(2)

    # Direct Riemann sum avoids FFT ringing on the sharp caustics
    aligned_density <- sapply(q_seq, function(q_target) {
      kernel <- dnorm(q_target - q_seq, mean = 0, sd = sigma_val)
      sum(raw_geometric_density * kernel) * step_size
    })

    dt_density <- data.table(q = q_seq)
    dt_density[, density := aligned_density]
    dt_density[, density := density / (sum(density) * step_size)]
    dt_density[, cdf := cumsum(density) * step_size]
    dt_density[, raw_density := raw_geometric_density]

    # Geometric Outlines
    df_circles <- data.frame(x0 = 0, y0 = 0, r_A = Delta_q, r_a = delta_q)
    df_cigars <- data.frame(x0 = 0, y0 = 0, aq_a = delta_q, aq_b = Delta_p, ap_a = Delta_q, ap_b = delta_p)

    # Column 1: Geometric Action Landscape & Apparatus
    p_ell <- ggplot(dt_wigner2d, aes(x=q, y=p)) +
      geom_circle(data = df_circles, aes(x0=x0, y0=y0, r=r_A), inherit.aes=FALSE, fill="white", color=NA) +
      geom_raster(aes(fill = w_plot), interpolate = TRUE) +
      scale_fill_gradient2(low = "gray40", mid = "white", high = "black", midpoint = 0, limits=c(-1, 1), guide="none") +
      geom_circle(data = df_circles, aes(x0=x0, y0=y0, r=r_A), inherit.aes=FALSE, color="black", linetype="solid", linewidth=0.3) +
      geom_ellipse(data = df_cigars, aes(x0=x0, y0=y0, a=aq_a, b=aq_b, angle=0), inherit.aes=FALSE, color="black", linetype="solid", linewidth=0.5) +
      geom_ellipse(data = df_cigars, aes(x0=x0, y0=y0, a=ap_a, b=ap_b, angle=0), inherit.aes=FALSE, color="black", linetype="solid", linewidth=0.5) +
      geom_circle(data = df_circles, aes(x0=x0, y0=y0, r=r_a), inherit.aes=FALSE, color="black", linetype="solid", linewidth=0.3) +
      coord_fixed(xlim=c(-plot_lim_x, plot_lim_x), ylim=c(-plot_lim_x, plot_lim_x), expand=FALSE) +
      scale_x_continuous(breaks = custom_breaks, labels = label_format) +
      scale_y_continuous(breaks = custom_breaks, labels = label_format) +
      theme_bw(base_family = base_font) +
      theme(panel.grid.minor=element_blank(), panel.background = element_rect(fill = "white"), axis.text=element_text(size=8), plot.margin = margin(2, 4, 2, 4)) +
      labs(x = ax_x_ell, y = ax_y_ell)

    # Column 2: Quantization Map
    p_stair <- ggplot(dt_density, aes(x=q, y=cdf)) +
      geom_line(color="black", linewidth=0.8) +
      coord_cartesian(xlim=c(-plot_lim_x, plot_lim_x), ylim=c(0, 1), expand=FALSE) +
      scale_x_continuous(breaks = custom_breaks, labels = label_format) +
      theme_bw(base_family=base_font) +
      theme(panel.grid.minor=element_blank(), axis.text=element_text(size=8), aspect.ratio=1, plot.margin = margin(2, 4, 2, 4)) +
      labs(x = ax_x_stair, y = ax_y_stair)

    # Column 3: Symplectic Distribution
    y_lim_den <- max(dt_density$density) * 1.5
    clip_y <- y_lim_den * 0.88

    # Isolate the three geometric paths heading toward the turning point asymptotes
    dt_left   <- dt_density[q < -Delta_q & raw_density <= clip_y]
    dt_center <- dt_density[q > -Delta_q & q < Delta_q & raw_density <= clip_y]
    dt_right  <- dt_density[q > Delta_q & raw_density <= clip_y]

    p_den <- ggplot() +
      # Base fill for the smoothed density
      geom_ribbon(data=dt_density, aes(x=q, ymin=0, ymax=density), fill="gray85", color=NA) +

      # Raw WKB Density with organic tangential arrows and TIGHT dash pattern ("22")
      geom_path(data=dt_left, aes(x=q, y=raw_density), linetype="22", color="gray40", linewidth=0.5,
                arrow = arrow(length = unit(0.12, "cm"), ends = "last", type="closed")) +
      geom_path(data=dt_center, aes(x=q, y=raw_density), linetype="22", color="gray40", linewidth=0.5,
                arrow = arrow(length = unit(0.12, "cm"), ends = "both", type="closed")) +
      geom_path(data=dt_right, aes(x=q, y=raw_density), linetype="22", color="gray40", linewidth=0.5,
                arrow = arrow(length = unit(0.12, "cm"), ends = "first", type="closed")) +

      # Infinity labels centered above the classical turning points
      annotate("text", x = -Delta_q, y = clip_y * 1.08, label = "infinity", parse = TRUE, color = "gray30", size = 4.5) +
      annotate("text", x = Delta_q, y = clip_y * 1.08, label = "infinity", parse = TRUE, color = "gray30", size = 4.5) +

      # Final convolved Symplectic Density line (Solid)
      geom_path(data=dt_density, aes(x=q, y=density), color="black", linewidth=0.4) +

      coord_cartesian(xlim=c(-plot_lim_x, plot_lim_x), ylim=c(0, y_lim_den), expand=FALSE) +
      scale_x_continuous(breaks = custom_breaks, labels = label_format) +
      theme_bw(base_family=base_font) +
      theme(panel.grid.minor=element_blank(), axis.text=element_text(size=8), axis.text.y=element_blank(), axis.ticks.y=element_blank(), aspect.ratio=1, plot.margin = margin(2, 4, 2, 4)) +
      labs(x = ax_x_den, y = ax_y_den)

    # Headers and Axis Label Management
    if (i == 1) {
      p_label <- p_label + labs(title = " ") + theme(plot.title=element_text(size=11, hjust=0.5, face="plain"))
      p_ell   <- p_ell   + labs(title = "Quantum of Action") + theme(plot.title=element_text(size=11, hjust=0.5, face=ifelse(is_ms, "plain", "bold")))
      p_stair <- p_stair + labs(title = "Quantization Map") + theme(plot.title=element_text(size=11, hjust=0.5, face=ifelse(is_ms, "plain", "bold")))
      p_den   <- p_den   + labs(title = "Symplectic Distribution") + theme(plot.title=element_text(size=11, hjust=0.5, face=ifelse(is_ms, "plain", "bold")))
    }

    # Hide X-axis titles on all but the bottom row
    if (i != num_rows) {
      p_ell   <- p_ell   + theme(axis.title.x = element_blank())
      p_stair <- p_stair + theme(axis.title.x = element_blank())
      p_den   <- p_den   + theme(axis.title.x = element_blank())
    }

    # Hide Y-axis titles on all but the middle row (Row 3) to prevent bleeding
    if (i != 3) {
      p_ell   <- p_ell   + theme(axis.title.y = element_blank())
      p_stair <- p_stair + theme(axis.title.y = element_blank())
      p_den   <- p_den   + theme(axis.title.y = element_blank())
    }

    plot_list <- c(plot_list, list(p_label, p_ell, p_stair, p_den))
  }

  final_plot <- wrap_plots(plot_list, ncol = 4, widths = c(0.25, 1, 1, 1))
  return(final_plot)
}

# --- 4. Execution ---
cat("Generating QHO geometric landscape and executing kinematic scans...\n")

p_final <- plot_qho_grid(
  dt_meta = dt_selected,
  style = "manuscript",
  base_font = latex_font
)

p_final <- p_final + theme(plot.margin = margin(10, 10, 10, 10))

# Adjusted width for 4-column layout mapping to PRL 2-column spread (~7.0 inches)
fig_width  <- 7.0
fig_height <- nrow(dt_selected) * 1.8 + 0.5

ggsave(
  filename = file_output_pdf,
  plot = p_final,
  device = cairo_pdf,
  width = fig_width,
  height = fig_height,
  limitsize = FALSE
)

cat("Success! Kinematic QHO geometric plots generated in", dir_figures, ":\n")
cat("PDF (Vector):", file_output_pdf, "\n")
