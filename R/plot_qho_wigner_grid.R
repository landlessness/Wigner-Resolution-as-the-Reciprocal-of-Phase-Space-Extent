# ==============================================================================
# plot_qho_wigner_grid.R  —  QHO Wigner cross-section grid
# Three columns:
#   Left:   W_n(q,p) heatmap with quantum of action
#   Center: W_n(q,0) cross-section — oscillates, goes negative
#   Right:  P_delta_q(q,0) cross-section — non-negative everywhere
# ==============================================================================

library(here)
library(data.table)
library(ggplot2)
library(ggforce)
library(patchwork)
library(gsl)

source(here::here("R", "symplectic_tools.R"))

latex_font      <- "CMU Serif"
dir_figures     <- here::here("figures")
if (!dir.exists(dir_figures)) dir.create(dir_figures, recursive=TRUE)
file_output_pdf <- file.path(dir_figures, "qho_wigner_grid.pdf")

target_n_levels <- c(0, 1, 2, 3, 24)
dt_selected     <- data.table(quantum_n=target_n_levels, action_A=2*target_n_levels+1)

plot_qho_wigner_grid <- function(dt_meta, base_font="") {
  plot_list <- list()
  num_rows  <- nrow(dt_meta)

  ax_y_cross <- expression(italic(W)[italic(n)](italic(q)*","*0))
  ax_y_conv  <- expression(italic(P)[delta*italic(q)](italic(q)*","*0))
  ax_x       <- expression(italic(q)/italic(q)[0])

  for (i in seq_len(num_rows)) {
    n_val <- dt_meta$quantum_n[i]
    alpha <- dt_meta$action_A[i]
    rs    <- qho_covariance(n_val)

    cat(sprintf("\nn=%d | A/A_0=%.0f | Delta_q=%.4f | delta_q=%.4f | RS: %s | SP: %s\n", n_val, rs$A_over_A0, rs$Delta_q, rs$delta_q, ifelse(rs$rs_satisfied,"OK","FAIL"), ifelse(rs$sp_satisfied,"OK","FAIL")))

    Delta_q <- rs$Delta_q

    # Display geometry from shared library
    geom          <- display_geometry(Delta_q)
    ell_lim       <- geom$ell_lim
    plot_lim      <- geom$plot_lim
    custom_breaks <- geom$custom_breaks
    label_format  <- geom$label_format

    # Row label from shared library
    p_label <- plot_row_label(sprintf("%.0f*A[0]", alpha), base_font=base_font)

    # Heatmap data from shared library
    dt_w2d <- wigner_heatmap_data(n_val, qho_wigner, ell_lim)

    # Center column: W_n(q,0) exact cross-section at p=0
    grid_res  <- max(400, 15*alpha)
    q_display <- seq(-plot_lim, plot_lim, length.out=grid_res)
    W_cross   <- qho_wigner(n_val, q_display, rep(0, length(q_display)))
    w_max     <- max(abs(W_cross), na.rm=TRUE)
    y_lim_w   <- w_max * 1.3

    # Right column: P_delta_q(q,0) — p=0 slice of convolved 2D distribution
    integ_lim <- max(Delta_q*2.0, 4.0)
    integ_res <- max(601, 20*ceiling(rs$A_over_A0))
    if (integ_res%%2==0) integ_res <- integ_res+1
    q_int <- seq(-integ_lim, integ_lim, length.out=integ_res)
    p_int <- seq(-integ_lim, integ_lim, length.out=integ_res)
    dq    <- diff(q_int)[1]; dp <- diff(p_int)[1]
    nq    <- length(q_int);  np_i <- length(p_int)

    W_mat <- outer(q_int, p_int, FUN=function(q,p) qho_wigner(n_val,q,p))
    K_mat <- outer(q_int, p_int, FUN=function(q,p) squeezed_kernel_q(q,p,rs))
    K_mat <- K_mat / (sum(K_mat)*dq*dp)

    ifftshift2d <- function(m) {
      nr <- nrow(m); nc <- ncol(m); sr <- floor(nr/2); sc <- floor(nc/2)
      rbind(cbind(m[(sr+1):nr,(sc+1):nc], m[(sr+1):nr,1:sc]), cbind(m[1:sr,(sc+1):nc], m[1:sr,1:sc]))
    }
    K_shift <- ifftshift2d(K_mat)
    P_mat   <- Re(fft(fft(W_mat)*fft(K_shift), inverse=TRUE)) / (nq*np_i) * dq*dp

    p0_idx      <- which.min(abs(p_int))
    P_cross_int <- P_mat[, p0_idx]
    P_cross     <- approx(q_int, P_cross_int, xout=q_display, rule=1)$y
    P_cross[is.na(P_cross)] <- 0

    # Scale convolved cross-section to raw Wigner amplitude for shape comparison
    p_max_display <- max(abs(P_cross), na.rm=TRUE)
    if (p_max_display > 0) P_cross_display <- P_cross / p_max_display * w_max else P_cross_display <- P_cross

    cat(sprintf("  W_n(0,0)=%.4f | scaled peak=%.4f\n", qho_wigner(n_val,0,0), max(abs(P_cross_display),na.rm=TRUE)))

    dt_cross <- data.table(q=q_display, W_raw=W_cross, P_conv=P_cross_display)

    # Symplectic ellipse overlays — for QHO system and Fermi boundaries coincide
    ell_data <- symplectic_ellipse_data(rs, r_system=Delta_q)

    # Column 1: heatmap from shared library
    p_ell <- plot_phase_space_heatmap(dt_w2d, ell_data, ell_lim, custom_breaks, label_format, base_font)

    # Column 2: W_n(q,0) — exact Wigner cross-section, goes negative
    p_cross_raw <- ggplot(dt_cross, aes(x=q, y=W_raw)) +
      geom_hline(yintercept=0, color="black", linewidth=0.3) +
      geom_ribbon(aes(ymin=pmin(W_raw,0), ymax=0), fill="gray60", alpha=0.6, color=NA) +
      geom_ribbon(aes(ymin=0, ymax=pmax(W_raw,0)), fill="gray85", color=NA) +
      geom_path(color="black", linewidth=0.4) +
      coord_cartesian(xlim=c(-ell_lim,ell_lim), ylim=c(-y_lim_w,y_lim_w), expand=FALSE) +
      scale_x_continuous(breaks=custom_breaks, labels=label_format) +
      theme_bw(base_family=base_font) +
      theme(panel.grid.minor=element_blank(), axis.text=element_text(size=8), axis.text.y=element_blank(), axis.ticks.y=element_blank(), aspect.ratio=1, plot.margin=margin(2,4,2,4)) +
      labs(x=ax_x, y=ax_y_cross)

    # Column 3: P_delta_q(q,0) — convolved 2D at p=0, non-negative
    p_cross_conv <- ggplot(dt_cross, aes(x=q, y=P_conv)) +
      geom_hline(yintercept=0, color="black", linewidth=0.3) +
      geom_ribbon(aes(ymin=0, ymax=pmax(P_conv,0)), fill="gray85", color=NA) +
      geom_path(color="black", linewidth=0.4) +
      coord_cartesian(xlim=c(-ell_lim,ell_lim), ylim=c(-y_lim_w,y_lim_w), expand=FALSE) +
      scale_x_continuous(breaks=custom_breaks, labels=label_format) +
      theme_bw(base_family=base_font) +
      theme(panel.grid.minor=element_blank(), axis.text=element_text(size=8), axis.text.y=element_blank(), axis.ticks.y=element_blank(), aspect.ratio=1, plot.margin=margin(2,4,2,4)) +
      labs(x=ax_x, y=ax_y_conv)

    if (i==1) {
      p_label      <- p_label      + labs(title=" ") + theme(plot.title=element_text(size=11, hjust=0.5))
      p_ell        <- p_ell        + labs(title="Phase-Space Cells") + theme(plot.title=element_text(size=11, hjust=0.5))
      p_cross_raw  <- p_cross_raw  + labs(title="Wigner Negativity") + theme(plot.title=element_text(size=11, hjust=0.5))
      p_cross_conv <- p_cross_conv + labs(title="Finite-Action Resolution") + theme(plot.title=element_text(size=11, hjust=0.5))
    }

    if (i != num_rows) {
      p_ell        <- p_ell        + theme(axis.title.x=element_blank())
      p_cross_raw  <- p_cross_raw  + theme(axis.title.x=element_blank())
      p_cross_conv <- p_cross_conv + theme(axis.title.x=element_blank())
    }

    if (i != 3) {
      p_ell        <- p_ell        + theme(axis.title.y=element_blank())
      p_cross_raw  <- p_cross_raw  + theme(axis.title.y=element_blank())
      p_cross_conv <- p_cross_conv + theme(axis.title.y=element_blank())
    }

    plot_list <- c(plot_list, list(p_label, p_ell, p_cross_raw, p_cross_conv))
  }

  wrap_plots(plot_list, ncol=4, widths=c(0.25,1,1,1))
}

cat("Computing Wigner cross-section grid...\n")
p_final <- plot_qho_wigner_grid(dt_meta=dt_selected, base_font=latex_font)
p_final <- p_final + theme(plot.margin=margin(10,10,10,10))
save_figure(p_final, file_output_pdf, nrow(dt_selected))
