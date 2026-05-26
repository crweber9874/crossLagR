# ============================================================================
# Bias in the CLPM
# Static figures: AR bias as a function of trait covariance and trait variance
# ============================================================================

library(MASS)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

# --- DGP: First-principles RI-CLPM (standardized) --------------------------

sim_riclpm <- function(waves, n, beta_x, beta_y, omega_xy, omega_yx,
                        var_p, var_q, cov_pq, var_BX, var_BY, cov_BXBY) {

  # Traits
  traits <- mvrnorm(n, mu = c(0, 0),
                    Sigma = matrix(c(var_BX, cov_BXBY, cov_BXBY, var_BY), 2))
  eta_x <- traits[, 1]; eta_y <- traits[, 2]

  # Within-person dynamics
  S <- matrix(c(var_p, cov_pq, cov_pq, var_q), 2)
  wp_x <- wp_y <- matrix(NA, n, waves)
  w1 <- mvrnorm(n, c(0, 0), S)
  wp_x[, 1] <- w1[, 1]; wp_y[, 1] <- w1[, 2]

  for (t in 2:waves) {
    e <- mvrnorm(n, c(0, 0), S)
    wp_x[, t] <- beta_x * wp_x[, t-1] + omega_yx * wp_y[, t-1] + e[, 1]
    wp_y[, t] <- beta_y * wp_y[, t-1] + omega_xy * wp_x[, t-1] + e[, 2]
  }

  # Observed = trait + state, then standardize each wave
  ox <- scale(wp_x + eta_x)
  oy <- scale(wp_y + eta_y)

  dat <- data.frame(ox, oy)
  names(dat) <- c(paste0("x", 1:waves), paste0("y", 1:waves))
  dat
}

# --- Estimator: demeaned FE (OLS on within-person deviations) ---------------

estimate_fe <- function(dat, waves) {
  xm <- as.matrix(dat[, paste0("x", 1:waves)])
  ym <- as.matrix(dat[, paste0("y", 1:waves)])
  xw <- xm - rowMeans(xm)
  yw <- ym - rowMeans(ym)

  rows <- lapply(2:waves, function(t)
    data.frame(wx = xw[, t], wy = yw[, t],
               wx_l = xw[, t-1], wy_l = yw[, t-1]))
  d <- bind_rows(rows)

  cx <- coef(lm(wx ~ wx_l + wy_l, d))
  cy <- coef(lm(wy ~ wx_l + wy_l, d))
  c(fe_ar_x = cx[["wx_l"]], fe_ar_y = cy[["wy_l"]],
    fe_cl_xy = cy[["wx_l"]], fe_cl_yx = cx[["wy_l"]])
}

# --- Estimator: pooled OLS (raw lagged regression) -------------------------

estimate_ols <- function(dat, waves) {
  xm <- as.matrix(dat[, paste0("x", 1:waves)])
  ym <- as.matrix(dat[, paste0("y", 1:waves)])

  rows <- lapply(2:waves, function(t)
    data.frame(x = xm[, t], y = ym[, t],
               xl = xm[, t-1], yl = ym[, t-1]))
  d <- bind_rows(rows)

  cx <- coef(lm(x ~ xl + yl, d))
  cy <- coef(lm(y ~ xl + yl, d))
  c(ols_ar_x = cx[["xl"]], ols_ar_y = cy[["yl"]],
    ols_cl_xy = cy[["xl"]], ols_cl_yx = cx[["yl"]])
}

# --- Monte Carlo wrapper ---------------------------------------------------

run_mc <- function(trials, waves, n, beta_x, beta_y, omega_xy, omega_yx,
                    var_p, var_q, cov_pq, var_BX, var_BY, cov_BXBY) {
  map_dfr(1:trials, function(i) {
    d <- sim_riclpm(waves, n, beta_x, beta_y, omega_xy, omega_yx,
                     var_p, var_q, cov_pq, var_BX, var_BY, cov_BXBY)
    fe  <- estimate_fe(d, waves)
    ols <- estimate_ols(d, waves)
    as.data.frame(t(c(fe, ols, trial = i)))
  })
}

# ============================================================================
# Figure 1: Vary trait COVARIANCE (cov_BXBY from -0.5 to 0.5)
# ============================================================================

cat("Figure 1: Varying trait covariance...\n")

grid_cov <- expand.grid(cov_BXBY = seq(-0.5, 0.5, by = 0.05))

results_cov <- map_dfr(1:nrow(grid_cov), function(i) {
  cat("\r  cov =", grid_cov$cov_BXBY[i])
  run_mc(trials = 1000, waves = 5, n = 500,
         beta_x = 0.3, beta_y = 0.3, omega_xy = 0, omega_yx = 0,
         var_p = 1, var_q = 1, cov_pq = 0.1,
         var_BX = 1, var_BY = 1, cov_BXBY = grid_cov$cov_BXBY[i]) %>%
    mutate(cov_BXBY = grid_cov$cov_BXBY[i])
})
cat("\n")

save(results_cov, file = "../data/bias_by_trait_covariance.RData")

# ============================================================================
# Figure 2: Vary trait VARIANCE (var_BX = var_BY from 0 to 2)
# ============================================================================

cat("Figure 2: Varying trait variance...\n")

grid_var <- expand.grid(var_B = seq(0, 2, by = 0.1))

results_var <- map_dfr(1:nrow(grid_var), function(i) {
  cat("\r  var =", grid_var$var_B[i])
  run_mc(trials = 1000, waves = 5, n = 500,
         beta_x = 0.3, beta_y = 0.3, omega_xy = 0, omega_yx = 0,
         var_p = 1, var_q = 1, cov_pq = 0.1,
         var_BX = grid_var$var_B[i], var_BY = grid_var$var_B[i],
         cov_BXBY = 0.3) %>%
    mutate(var_B = grid_var$var_B[i])
})
cat("\n")

save(results_var, file = "../data/bias_by_trait_variance.RData")

# ============================================================================
# Plot theme
# ============================================================================

theme_bias <- function() {
  theme_minimal(base_size = 13, base_family = "Helvetica") +
    theme(
      plot.title       = element_text(face = "bold", size = 15, hjust = 0),
      plot.subtitle    = element_text(color = "grey40", size = 10, hjust = 0,
                                       margin = margin(b = 12)),
      plot.caption     = element_text(color = "grey50", size = 8, hjust = 1),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
      axis.title       = element_text(face = "bold", size = 11),
      axis.text        = element_text(size = 10),
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold", size = 10),
      legend.text      = element_text(size = 9),
      strip.text       = element_text(face = "bold", size = 11),
      plot.margin      = margin(15, 15, 10, 15)
    )
}

# ============================================================================
# Figure 1 plot
# ============================================================================

summ_cov <- results_cov %>%
  pivot_longer(cols = c(ols_ar_x, ols_ar_y, fe_ar_x, fe_ar_y),
               names_to = "param", values_to = "est") %>%
  mutate(
    Estimator = ifelse(grepl("^ols", param), "Pooled OLS", "Fixed Effects"),
    Variable  = ifelse(grepl("_x$", param), "AR(X)", "AR(Y)")
  ) %>%
  group_by(cov_BXBY, Estimator, Variable) %>%
  summarise(mean = mean(est), lo = quantile(est, 0.025),
            hi = quantile(est, 0.975), .groups = "drop")

p1 <- ggplot(summ_cov, aes(x = cov_BXBY, y = mean,
                             color = Estimator, fill = Estimator)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0.3, linetype = "dashed", color = "grey30",
             linewidth = 0.6) +
  annotate("text", x = -0.48, y = 0.31, label = "True AR = 0.30",
           hjust = 0, size = 3.2, color = "grey30", fontface = "italic") +
  facet_wrap(~Variable) +
  scale_color_manual(values = c("Pooled OLS" = "#E63946", "Fixed Effects" = "#457B9D")) +
  scale_fill_manual(values = c("Pooled OLS" = "#E63946", "Fixed Effects" = "#457B9D")) +
  labs(
    title    = "Autoregressive Bias by Trait Covariance",
    subtitle = expression(paste("DGP: RI-CLPM | ", beta[xx] == beta[yy] == 0.30,
                                 " | Var(trait) = 1 | 1,000 MC draws | N = 500 | T = 5")),
    x = expression(paste("Trait Covariance  ", Cov(eta[x], eta[y]))),
    y = "Mean AR Estimate",
    caption  = "Ribbon = 95% Monte Carlo interval"
  ) +
  theme_bias()

ggsave("../data/fig1_bias_by_trait_covariance.pdf", p1, width = 9, height = 5)
ggsave("../data/fig1_bias_by_trait_covariance.png", p1, width = 9, height = 5, dpi = 300)

cat("Saved Figure 1\n")

# ============================================================================
# Figure 2 plot
# ============================================================================

summ_var <- results_var %>%
  pivot_longer(cols = c(ols_ar_x, ols_ar_y, fe_ar_x, fe_ar_y),
               names_to = "param", values_to = "est") %>%
  mutate(
    Estimator = ifelse(grepl("^ols", param), "Pooled OLS", "Fixed Effects"),
    Variable  = ifelse(grepl("_x$", param), "AR(X)", "AR(Y)")
  ) %>%
  group_by(var_B, Estimator, Variable) %>%
  summarise(mean = mean(est), lo = quantile(est, 0.025),
            hi = quantile(est, 0.975), .groups = "drop")

p2 <- ggplot(summ_var, aes(x = var_B, y = mean,
                             color = Estimator, fill = Estimator)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0.3, linetype = "dashed", color = "grey30",
             linewidth = 0.6) +
  annotate("text", x = 0.05, y = 0.31, label = "True AR = 0.30",
           hjust = 0, size = 3.2, color = "grey30", fontface = "italic") +
  facet_wrap(~Variable) +
  scale_color_manual(values = c("Pooled OLS" = "#E63946", "Fixed Effects" = "#457B9D")) +
  scale_fill_manual(values = c("Pooled OLS" = "#E63946", "Fixed Effects" = "#457B9D")) +
  scale_x_continuous(breaks = seq(0, 2, 0.5)) +
  labs(
    title    = "Autoregressive Bias by Trait Variance",
    subtitle = expression(paste("DGP: RI-CLPM | ", beta[xx] == beta[yy] == 0.30,
                                 " | Cov(trait) = 0.30 | 1,000 MC draws | N = 500 | T = 5")),
    x = expression(paste("Trait Variance  ", Var(eta[x]) == Var(eta[y]))),
    y = "Mean AR Estimate",
    caption  = "Ribbon = 95% Monte Carlo interval"
  ) +
  theme_bias()

ggsave("../data/fig2_bias_by_trait_variance.pdf", p2, width = 9, height = 5)
ggsave("../data/fig2_bias_by_trait_variance.png", p2, width = 9, height = 5, dpi = 300)

cat("Saved Figure 2\nDone.\n")
