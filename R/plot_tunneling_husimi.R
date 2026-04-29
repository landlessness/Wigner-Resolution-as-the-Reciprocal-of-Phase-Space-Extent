# ==============================================================================
# plot_tunneling_husimi.R
# Appendix: Husimi tunneling figure
#
# Single panel. q on x. log density on y.
#
# Black curves (the Husimi resolution):
#   - rho_Q^Wigner(q)  Wigner input through Husimi kernel    (solid)
#   - rho_Q^cl(q)      Classical input through Husimi kernel (dashed)
#
# Gray reference curves (other formalisms):
#   - |psi_0(q)|^2  exact Schrodinger             (solid)
#   - P_WKB(q)      classical caustic, zero in barrier (dotdash)
#
# Same lowest-eigenstate target as the symplectic figure for direct comparison.
# ==============================================================================

library(here)
library(data.table)
library(ggplot2)

source(here("R", "plot_tools.R"))
source(here("R", "double_well_potential.R"))
source(here("R", "wigner_tools.R"))
source(here("R", "wigner_state.R"))
source(here("R", "semiclassical_state.R"))
source(here("R", "math_tools.R"))
source(here("R", "husimi_kernel.R"))
source(here("R", "classical_action_tools.R"))

latex_font      <- "CMU Serif"
dir_figures     <- here("figures")
if (!dir.exists(dir_figures)) dir.create(dir_figures, recursive=TRUE)
file_output_pdf <- file.path(dir_figures, "tunneling_husimi.pdf")

n_target  <- 0  # lowest tunneling doublet partner

# ------------------------------------------------------------------------------
# Solve for psi_0 and key quantities
# ------------------------------------------------------------------------------

cat("Solving double-well Schrodinger equation...\n")
dw_soln <- solve_schrodinger(double_well_V,
                             q_min=DOUBLE_WELL_Q_MIN,
                             q_max=DOUBLE_WELL_Q_MAX,
                             dq=DOUBLE_WELL_DQ,
                             n_states=DOUBLE_WELL_N_STATES)

E_n     <- dw_soln$energies[n_target+1]
q_grid  <- dw_soln$q_grid
psi_vec <- dw_soln$psi_matrix[, n_target+1]
tp      <- double_well_turning_points(E_n)
roots   <- tp$roots
q_minus <- roots[1]
q_plus  <- roots[4]

rs <- numerical_covariance(psi_vec, q_grid)
cat(sprintf("\nn=%d | E=%.4f | turning points: %.3f, %.3f, %.3f, %.3f\n",
            n_target, E_n, roots[1], roots[2], roots[3], roots[4]))

# ------------------------------------------------------------------------------
# Display window (matches symplectic version for direct comparison)
# ------------------------------------------------------------------------------

q_extent <- max(abs(q_minus), abs(q_plus))
q_pad    <- q_extent * 0.15
q_lo     <- -(q_extent + q_pad)
q_hi     <-   q_extent + q_pad

V_min <- min(double_well_V(q_grid))
p_max <- sqrt(2 * (E_n - V_min))
p_lo  <- -(p_max * 1.3)
p_hi  <-   p_max * 1.3

q_display <- seq(q_lo, q_hi, length.out=800)

# ------------------------------------------------------------------------------
# Curve 1: |psi_0(q)|^2 (exact Schrodinger)
# ------------------------------------------------------------------------------

psi_density_full <- abs(psi_vec)^2
psi_density_full <- psi_density_full / sum(psi_density_full * (q_grid[2]-q_grid[1]))
psi_density <- approx(q_grid, psi_density_full, xout=q_display, rule=1)$y
psi_density[is.na(psi_density)] <- 0

# ------------------------------------------------------------------------------
# Curve 2: Husimi resolution of Wigner input (full marginal over p)
# ------------------------------------------------------------------------------

cat("Building Wigner state and Husimi-Wigner marginal...\n")
w_state <- build_wigner_state(psi_vec, q_grid,
                              q_lo, q_hi, p_lo, p_hi, q_display)
K_husimi_w <- husimi_kernel_matrix(w_state$q_int, w_state$p_int)
conv_w     <- fft_convolve_2d(w_state$W_matrix, K_husimi_w,
                              w_state$dq_int, w_state$dp_int)
rho_husimi_int <- rowSums(conv_w$P_mat) * w_state$dp_int
rho_husimi <- approx(w_state$q_int, rho_husimi_int, xout=q_display, rule=1)$y
rho_husimi[is.na(rho_husimi)] <- 0
rho_husimi <- rho_husimi / sum(rho_husimi * (q_display[2]-q_display[1]))

# ------------------------------------------------------------------------------
# Curve 3: Husimi resolution of classical input (full marginal over p)
# ------------------------------------------------------------------------------

cat("Building semiclassical state and Husimi-classical marginal...\n")
sc_state <- build_semiclassical_state(E_n, double_well_V,
                                      q_lo, q_hi, p_lo, p_hi, q_display)
K_husimi_sc <- husimi_kernel_matrix(sc_state$q_int, sc_state$p_int)
conv_cl     <- fft_convolve_2d(sc_state$W_matrix, K_husimi_sc,
                               sc_state$dq_int, sc_state$dp_int)
rho_cl_husimi_int <- rowSums(conv_cl$P_mat) * sc_state$dp_int
rho_cl_husimi <- approx(sc_state$q_int, rho_cl_husimi_int, xout=q_display, rule=1)$y
rho_cl_husimi[is.na(rho_cl_husimi)] <- 0
rho_cl_husimi <- rho_cl_husimi / sum(rho_cl_husimi * (q_display[2]-q_display[1]))

# ------------------------------------------------------------------------------
# Curve 4: WKB caustic density (analytical)
# ------------------------------------------------------------------------------

wkb_density <- sc_state$wkb_density

# ------------------------------------------------------------------------------
# Background: V(q) scaled into the plot region
# ------------------------------------------------------------------------------

V_q <- double_well_V(q_display)

# ------------------------------------------------------------------------------
# Build the data table and plot
# ------------------------------------------------------------------------------

y_peak  <- max(psi_density, rho_husimi, rho_cl_husimi, na.rm=TRUE)
y_floor <- y_peak * 1e-6
y_top   <- y_peak * 2.0

clip_low <- function(x, floor) ifelse(x < floor | !is.finite(x), NA_real_, x)

dt_curves <- rbind(
  data.table(q=q_display, density=clip_low(psi_density, y_floor),      curve="psi"),
  data.table(q=q_display, density=clip_low(rho_husimi, y_floor),       curve="rho_husimi"),
  data.table(q=q_display, density=clip_low(rho_cl_husimi, y_floor),    curve="rho_cl_husimi"),
  data.table(q=q_display, density=clip_low(wkb_density, y_floor),      curve="wkb")
)
dt_curves[, curve := factor(curve,
                            levels=c("psi","rho_husimi","rho_cl_husimi","wkb"))]

V_min_disp <- min(V_q)
V_max_disp <- max(V_q)
V_scaled   <- y_floor * (y_peak / y_floor)^((V_q - V_min_disp) /
                                              (V_max_disp - V_min_disp))
dt_potential <- data.table(q=q_display, V_scaled=V_scaled)

tp_lines <- data.table(q=roots)

custom_breaks_q <- round(c(roots[1], roots[2], 0, roots[3], roots[4]), 1)
label_format    <- function(x) sprintf("%.1f", x)

# ------------------------------------------------------------------------------
# Compose plot
#
# Layering, back to front:
#   1. faint V(q) potential
#   2. vertical turning-point guides
#   3. gray reference curves: psi (solid), WKB (dotdash)
#   4. black Husimi curves: rho_cl_husimi (dashed), rho_husimi (solid)
# Black curves drawn last so they sit on top.
# ------------------------------------------------------------------------------

p <- ggplot() +
  geom_path(data=dt_potential, aes(x=q, y=V_scaled),
            color="gray85", linewidth=0.4) +
  geom_vline(data=tp_lines, aes(xintercept=q),
             color="gray70", linewidth=0.3, linetype="dotted") +
  # Gray reference curves
  geom_path(data=dt_curves[curve=="psi"],
            aes(x=q, y=density), color="gray50",
            linewidth=0.5, linetype="solid", na.rm=TRUE) +
  geom_path(data=dt_curves[curve=="wkb"],
            aes(x=q, y=density), color="gray50",
            linewidth=0.5, linetype="dotdash", na.rm=TRUE) +
  # Black Husimi curves (drawn last for visual prominence)
  geom_path(data=dt_curves[curve=="rho_cl_husimi"],
            aes(x=q, y=density), color="black",
            linewidth=0.6, linetype="dashed", na.rm=TRUE) +
  geom_path(data=dt_curves[curve=="rho_husimi"],
            aes(x=q, y=density), color="black",
            linewidth=0.7, linetype="solid", na.rm=TRUE) +
  scale_x_continuous(breaks=custom_breaks_q, labels=label_format) +
  scale_y_log10(limits=c(y_floor, y_top), expand=c(0,0)) +
  coord_cartesian(xlim=c(q_lo, q_hi), expand=FALSE) +
  theme_bw(base_family=latex_font) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.y=element_line(color="gray95", linewidth=0.3),
        panel.grid.major.x=element_blank(),
        axis.text=element_text(size=9),
        axis.title=element_text(size=10),
        plot.margin=margin(8, 12, 8, 12)) +
  labs(x=expression(italic(q)/italic(q)[0]),
       y=expression(rho(italic(q))))

ggsave(filename=file_output_pdf, plot=p, device=cairo_pdf,
       width=7.0, height=4.0, limitsize=FALSE)
cat("Done. ", file_output_pdf, "\n")
