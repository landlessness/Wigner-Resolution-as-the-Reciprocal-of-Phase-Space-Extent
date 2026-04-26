# ==============================================================================
# plot_qho_berry_grid.R  —  Berry semiclassical cross-section grid
# Three columns:
#   Left:   W_Berry(q,p) heatmap with quantum of action
#   Center: W_Berry(q,0) cross-section — diverges at turning points
#   Right:  P_delta_q(q,0) cross-section of convolved 2D distribution
# ==============================================================================

library(here)
library(data.table)
library(ggplot2)
library(ggforce)
library(patchwork)

source(here::here("R", "symplectic_tools.R"))

latex_font      <- "CMU Serif"
dir_figures     <- here::here("figures")
if (!dir.exists(dir_figures)) dir.create(dir_figures, recursive=TRUE)
file_output_pdf <- file.path(dir_figures, "qho_berry_grid.pdf")

target_n_levels <- c(0, 1, 2, 3, 24)
dt_selected     <- data.table(quantum_n=target_n_levels, action_A=2*target_n_levels+1)

# Berry physics functions

chord_area <- function(r, R) R^2*acos(pmin(r/R,1)) - r*sqrt(pmax(R^2-r^2,0))

berry_wigner <- function(n, q, p) {
  R <- sqrt(2*n+1); r <- sqrt(q^2+p^2); val <- numeric(length(r))
  idx_in <- r < R
  if (any(idx_in)) { r_in <- r[idx_in]; Area <- chord_area(r_in,R); denom <- pmax((R^2-r_in^2)^0.25,1e-6); val[idx_in] <- (1/(pi*denom))*cos(Area-pi/4) }
  idx_out <- r >= R
  if (any(idx_out)) { r_out <- r[idx_out]; Area_out <- r_out*sqrt(r_out^2-R^2)-R^2*acosh(r_out/R); denom <- pmax((r_out^2-R^2)^0.25,1e-6); val[idx_out] <- (1/(2*pi*denom))*exp(-Area_out) }
  return(val)
}

# Berry cross-section W_Berry(q,0) at p=0 — diverges at turning points
berry_cross_section <- function(n, q_seq) {
  R <- sqrt(2*n+1); val <- rep(0, length(q_seq))
  idx_in <- abs(q_seq) < R
  if (any(idx_in)) { q_in <- q_seq[idx_in]; Area <- chord_area(abs(q_in),R); denom <- pmax(sqrt(R^2-q_in^2),1e-6); val[idx_in] <- (2/(pi*denom))*cos(Area/2-pi/4)^2 }
  idx_out <- abs(q_seq) >= R
  if (any(idx_out)) { q_out <- q_seq[idx_out]; Area_out <- abs(q_out)*sqrt(q_out^2-R^2)-R^2*acosh(abs(q_out)/R); denom <- pmax(sqrt(q_out^2-R^2),1e-6); val[idx_out] <- (1/(2*pi*denom))*exp(-2*Area_out) }
  return(val)
}

# Convolve Berry cross-section with 1D marginal of squeezed kernel
convolve_berry_cross <- function(berry_raw, q_seq, delta_q) {
  dq <- diff(q_seq)[1]; n_pts <- length(q_seq)
  G_1d <- exp(-q_seq^2/delta_q^2); G_1d <- G_1d/(sum(G_1d)*dq)
  P_conv <- convolve(berry_raw, rev(G_1d), type="open")
  half <- floor(n_pts/2); start <- half+1
  P_out <- P_conv[start:(start+n_pts-1)]
  P_out <- pmax(P_out,0); P_out/(sum(P_out)*dq)
}

plot_berry_grid <- function(dt_meta, base_font="") {
  plot_list <- list()
  num_rows  <- nrow(dt_meta)

  ax_x     <- expression(italic(q)/italic(q)[0])
  ax_y     <- expression(italic(p)/italic(p)[0])
  ax_y_raw <- expression(italic(W)[Berry](italic(q)*","*0))
  ax_y_den <- expression(italic(P)[delta*italic(q)](italic(q)*","*0))

  for (i in seq_len(num_rows)) {
    n_val <- dt_meta$quantum_n[i]
    alpha <- dt_meta$action_A[i]
    rs    <- qho_covariance(n_val)

    cat(sprintf("\nn=%d | A/A_0=%.0f | Delta_q=%.4f | delta_q=%.4f | RS: %s | SP: %s\n", n_val, rs$A_over_A0, rs$Delta_q, rs$delta_q, ifelse(rs$rs_satisfied,"OK","FAIL"), ifelse(rs$sp_satisfied,"OK","FAIL")))

    Delta_q <- rs$Delta_q; Delta_p <- rs$Delta_p
    delta_q <- rs$delta_q; delta_p <- rs$delta_p

    ell_lim       <- Delta_q * 1.25
    plot_lim      <- max(Delta_q*1.25, Delta_q+2.5)
    break_val     <- min(round(Delta_q,1), floor(ell_lim*10)/10)
    custom_breaks <- c(-break_val, 0, break_val)
    label_format  <- function(x) sprintf("%.1f", x)

    p_label <- ggplot() + theme_void() + coord_cartesian(xlim=c(-1,1), ylim=c(0,1), clip="off") + annotate("text", x=0, y=0.5, label=sprintf("%.0f*A[0]", alpha), parse=TRUE, family=base_font, size=4.5, hjust=0.5)

    grid_res  <- max(400, 15*alpha)
    q_display <- seq(-plot_lim, plot_lim, length.out=grid_res)
    dq_disp   <- diff(q_display)[1]

    # Berry 2D heatmap
    q_ell     <- seq(-ell_lim, ell_lim, length.out=400)
    p_ell_seq <- seq(-ell_lim, ell_lim, length.out=400)
    dt_w2d    <- as.data.table(expand.grid(q=q_ell, p=p_ell_seq))
    dt_w2d[, w := berry_wigner(n_val, q, p)]
    dt_w2d[, w_plot := sign(w)*abs(w)^0.4]
    max_w <- max(abs(dt_w2d$w_plot), na.rm=TRUE)
    if (max_w > 0) dt_w2d[, w_plot := w_plot/max_w]

    # Center: W_Berry(q,0) cross-section — diverges at turning points
    berry_raw <- berry_cross_section(n_val, q_display)

    # Right: P_delta_q(q,0) — convolved, finite everywhere
    P_berry_conv <- convolve_berry_cross(berry_raw, q_display, delta_q)

    cat(sprintf("  Berry conv peak: %.4f at q=%.4f\n", max(P_berry_conv), q_display[which.max(P_berry_conv)]))

    # Scale berry_raw for display to match convolved interior amplitude
    interior_mask <- abs(q_display) < Delta_q*0.6 & berry_raw < quantile(berry_raw[abs(q_display) < Delta_q], 0.7, na.rm=TRUE)
    if (any(interior_mask) && any(P_berry_conv[interior_mask] > 0)) {
      scale_factor <- mean(P_berry_conv[interior_mask], na.rm=TRUE) / max(mean(berry_raw[interior_mask], na.rm=TRUE), 1e-10)
    } else {
      scale_factor <- 1 / max(berry_raw[abs(q_display) < Delta_q*0.9], na.rm=TRUE)
    }
    berry_raw_display <- berry_raw * scale_factor

    density_peak <- max(P_berry_conv, na.rm=TRUE)
    y_lim_den    <- density_peak * 1.4
    clip_y       <- density_peak * 1.2

    dt_density <- data.table(q=q_display, density=P_berry_conv, berry_raw_display=berry_raw_display)

    dt_forced <- rbind(dt_density[abs(q) < Delta_q], data.table(q=-Delta_q, density=0, berry_raw_display=y_lim_den*0.88), data.table(q=Delta_q, density=0, berry_raw_display=y_lim_den*0.88))[order(q)]

    dt_center <- dt_density[abs(q) < Delta_q & berry_raw_display <= clip_y]

    # Symplectic ellipse overlays — for QHO r_system = Delta_q so boundaries coincide
    ell_data <- symplectic_ellipse_data(rs, r_system=Delta_q)

    # Column 1: Berry 2D heatmap
    # Bottom layers (white fill) drawn before raster, top layers (outlines) after
    p_ell <- ggplot(dt_w2d, aes(x=q, y=p)) +
      symplectic_ellipse_layers_bottom(ell_data) +
      geom_raster(aes(fill=w_plot), interpolate=TRUE) +
      scale_fill_gradient2(low="gray10", mid="white", high="gray40", midpoint=0, limits=c(-1,1), guide="none") +
      symplectic_ellipse_layers_top(ell_data) +
      coord_fixed(xlim=c(-ell_lim,ell_lim), ylim=c(-ell_lim,ell_lim), expand=FALSE) +
      scale_x_continuous(breaks=custom_breaks, labels=label_format) +
      scale_y_continuous(breaks=custom_breaks, labels=label_format) +
      theme_bw(base_family=base_font) +
      theme(panel.grid.minor=element_blank(), panel.background=element_rect(fill="white"), axis.text=element_text(size=8), plot.margin=margin(2,4,2,4)) +
      labs(x=ax_x, y=ax_y)

    # Column 2: W_Berry(q,0) cross-section — diverges at turning points
    p_raw <- ggplot() +
      geom_ribbon(data=dt_forced, aes(x=q, ymin=0, ymax=pmin(berry_raw_display, y_lim_den*0.88)), fill="gray85", color=NA) +
      geom_path(data=dt_center[order(q)], aes(x=q, y=berry_raw_display), color="black", linewidth=0.4, arrow=arrow(length=unit(0.12,"cm"), ends="both", type="closed")) +
      annotate("text", x=-Delta_q, y=y_lim_den*0.92, label="infinity", parse=TRUE, color="gray30", size=4.5) +
      annotate("text", x= Delta_q, y=y_lim_den*0.92, label="infinity", parse=TRUE, color="gray30", size=4.5) +
      coord_cartesian(xlim=c(-ell_lim,ell_lim), ylim=c(0,y_lim_den), expand=FALSE) +
      scale_x_continuous(breaks=custom_breaks, labels=label_format) +
      theme_bw(base_family=base_font) +
      theme(panel.grid.minor=element_blank(), axis.text=element_text(size=8), axis.text.y=element_blank(), axis.ticks.y=element_blank(), aspect.ratio=1, plot.margin=margin(2,4,2,4)) +
      labs(x=ax_x, y=ax_y_raw)

    # Column 3: P_delta_q(q,0) — convolved 2D at p=0, finite everywhere
    p_den <- ggplot() +
      geom_ribbon(data=dt_density, aes(x=q, ymin=0, ymax=density), fill="gray85", color=NA) +
      geom_path(data=dt_density, aes(x=q, y=density), color="black", linewidth=0.4) +
      coord_cartesian(xlim=c(-ell_lim,ell_lim), ylim=c(0,y_lim_den), expand=FALSE) +
      scale_x_continuous(breaks=custom_breaks, labels=label_format) +
      theme_bw(base_family=base_font) +
      theme(panel.grid.minor=element_blank(), axis.text=element_text(size=8), axis.text.y=element_blank(), axis.ticks.y=element_blank(), aspect.ratio=1, plot.margin=margin(2,4,2,4)) +
      labs(x=ax_x, y=ax_y_den)

    if (i==1) {
      p_label <- p_label + labs(title=" ") + theme(plot.title=element_text(size=11, hjust=0.5))
      p_ell   <- p_ell   + labs(title="Quantum of Action") + theme(plot.title=element_text(size=11, hjust=0.5))
      p_raw   <- p_raw   + labs(title="Semiclassical Caustics") + theme(plot.title=element_text(size=11, hjust=0.5))
      p_den   <- p_den   + labs(title="Finite-Action Resolution") + theme(plot.title=element_text(size=11, hjust=0.5))
    }

    if (i != num_rows) {
      p_ell <- p_ell + theme(axis.title.x=element_blank())
      p_raw <- p_raw + theme(axis.title.x=element_blank())
      p_den <- p_den + theme(axis.title.x=element_blank())
    }

    if (i != 3) {
      p_ell <- p_ell + theme(axis.title.y=element_blank())
      p_raw <- p_raw + theme(axis.title.y=element_blank())
      p_den <- p_den + theme(axis.title.y=element_blank())
    }

    plot_list <- c(plot_list, list(p_label, p_ell, p_raw, p_den))
  }

  wrap_plots(plot_list, ncol=4, widths=c(0.25,1,1,1))
}

cat("Computing Berry cross-section grid...\n")
p_final <- plot_berry_grid(dt_meta=dt_selected, base_font=latex_font)
p_final <- p_final + theme(plot.margin=margin(10,10,10,10))
ggsave(filename=file_output_pdf, plot=p_final, device=cairo_pdf, width=7.0, height=nrow(dt_selected)*1.8+0.5, limitsize=FALSE)
cat("Done.", file_output_pdf, "\n")
