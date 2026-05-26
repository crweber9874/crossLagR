#' @title sim_trait_state
#' @description Simulate panel data from a trait-plus-state data-generating process.
#'
#' @param waves Integer. Number of waves (time points). Must be at least 2.
#' @param n Integer. Sample size (number of individuals).
#' @param beta_x Numeric. Autoregressive effect for X.
#' @param beta_y Numeric. Autoregressive effect for Y.
#' @param omega_xy Numeric. Cross-lagged effect from X(t-1) to Y(t).
#' @param omega_yx Numeric. Cross-lagged effect from Y(t-1) to X(t).
#' @param var_p Numeric. Innovation variance for within-person X process.
#' @param var_q Numeric. Innovation variance for within-person Y process.
#' @param cov_pq Numeric. Innovation covariance between within-person X and Y processes.
#' @param var_BX Numeric. Between-person trait variance for X.
#' @param var_BY Numeric. Between-person trait variance for Y.
#' @param cov_BXBY Numeric. Between-person trait covariance between X and Y.
#' @param cor_trait_initial Numeric. Correlation between traits and wave-1
#'   within-person states. Default is 0 (independent, which satisfies the RI-CLPM
#'   identification assumption). Non-zero values violate this assumption and can
#'   reveal bias in models that assume independence. Applied symmetrically:
#'   cor(eta_x, wp_x1) = cor(eta_y, wp_y1) = cor_trait_initial.
#'
#' @return A wide-format data frame with columns \code{x1, ..., x[waves]} and
#'   \code{y1, ..., y[waves]}.
#'
#' @details
#' The function simulates two stable traits (between-person random intercept components)
#' and two within-person AR(1) processes with optional cross-lagged effects. Observed scores
#' are the sum of trait and state components.
#'
#' When \code{cor_trait_initial = 0} (default), traits and wave-1 states are drawn
#' independently. When non-zero, they are drawn jointly from a multivariate normal,
#' so individuals with higher trait values tend to also have higher (or lower) initial
#' state values. This violates the RI-CLPM identification constraint that random
#' intercepts are uncorrelated with wave-1 within-person factors.
#'
#' @examples
#' dat <- sim_trait_state(
#'   waves = 5,
#'   n = 500,
#'   beta_x = 0.4,
#'   beta_y = 0.3,
#'   omega_xy = 0.1,
#'   omega_yx = 0.1,
#'   var_p = 1,
#'   var_q = 1,
#'   cov_pq = 0.2,
#'   var_BX = 0.6,
#'   var_BY = 0.6,
#'   cov_BXBY = 0.3
#' )
#' head(dat)
#'
#' @export
sim_trait_state <- function(waves = 5,
                            n = 1000,
                            beta_x = 0,
                            beta_y = 0,
                            omega_xy = 0,
                            omega_yx = 0,
                            var_p = 1,
                            var_q = 1,
                            cov_pq = 0,
                            var_BX = 0,
                            var_BY = 0,
                            cov_BXBY = 0,
                            cor_trait_initial = 0) {

  if (!is.numeric(waves) || length(waves) != 1 || waves != as.integer(waves) || waves < 2) {
    stop("Error: 'waves' must be a single integer >= 2.")
  }

  if (!is.numeric(n) || length(n) != 1 || n != as.integer(n) || n <= 0) {
    stop("Error: 'n' must be a single positive integer.")
  }

  numeric_params <- list(
    beta_x = beta_x,
    beta_y = beta_y,
    omega_xy = omega_xy,
    omega_yx = omega_yx,
    var_p = var_p,
    var_q = var_q,
    cov_pq = cov_pq,
    var_BX = var_BX,
    var_BY = var_BY,
    cov_BXBY = cov_BXBY,
    cor_trait_initial = cor_trait_initial
  )

  for (param_name in names(numeric_params)) {
    value <- numeric_params[[param_name]]
    if (!is.numeric(value) || length(value) != 1 || is.na(value)) {
      stop(paste0("Error: '", param_name, "' must be a single non-missing numeric value."))
    }
  }

  if (var_p <= 0 || var_q <= 0 || var_BX < 0 || var_BY < 0) {
    stop("Error: 'var_p' and 'var_q' must be > 0, and 'var_BX' and 'var_BY' must be >= 0.")
  }

  if (abs(cor_trait_initial) >= 1) {
    stop("Error: 'cor_trait_initial' must be between -1 and 1 (exclusive).")
  }

  # Within-person innovation covariance matrix (used for waves 2+ and wave 1 if independent)
  innovation_sigma <- matrix(c(var_p, cov_pq, cov_pq, var_q), 2)
  wp_x <- wp_y <- matrix(NA_real_, nrow = n, ncol = waves)

  if (cor_trait_initial == 0 || (var_BX == 0 && var_BY == 0)) {
    # Independent: draw traits and wave-1 states separately
    traits <- MASS::mvrnorm(
      n,
      mu = c(0, 0),
      Sigma = matrix(c(var_BX, cov_BXBY, cov_BXBY, var_BY), 2)
    )
    eta_x <- traits[, 1]
    eta_y <- traits[, 2]

    w1 <- MASS::mvrnorm(n, c(0, 0), innovation_sigma)
    wp_x[, 1] <- w1[, 1]
    wp_y[, 1] <- w1[, 2]
  } else {
    # Joint draw: traits correlated with wave-1 within-person states
    # Order: eta_x, eta_y, wp_x1, wp_y1
    # Same-variable trait-state covariances
    cov_etax_wpx1 <- cor_trait_initial * sqrt(var_BX * var_p)
    cov_etay_wpy1 <- cor_trait_initial * sqrt(var_BY * var_q)
    # Cross-variable trait-state covariances (scaled by both trait and state correlations)
    cov_etax_wpy1 <- cor_trait_initial * sqrt(var_BX * var_q) *
      (cov_BXBY / sqrt(var_BX * var_BY)) * (cov_pq / sqrt(var_p * var_q))
    cov_etay_wpx1 <- cor_trait_initial * sqrt(var_BY * var_p) *
      (cov_BXBY / sqrt(var_BX * var_BY)) * (cov_pq / sqrt(var_p * var_q))

    joint_sigma <- matrix(c(
      var_BX,          cov_BXBY,        cov_etax_wpx1,   cov_etax_wpy1,
      cov_BXBY,        var_BY,          cov_etay_wpx1,   cov_etay_wpy1,
      cov_etax_wpx1,   cov_etay_wpx1,   var_p,           cov_pq,
      cov_etax_wpy1,   cov_etay_wpy1,   cov_pq,          var_q
    ), nrow = 4, byrow = TRUE)

    # Check positive definiteness
    eig <- eigen(joint_sigma, symmetric = TRUE, only.values = TRUE)$values
    if (any(eig <= 0)) {
      stop("Error: joint covariance matrix of traits and wave-1 states is not positive ",
           "definite. Try a smaller 'cor_trait_initial' value (current: ",
           cor_trait_initial, ").")
    }

    joint <- MASS::mvrnorm(n, mu = rep(0, 4), Sigma = joint_sigma)
    eta_x <- joint[, 1]
    eta_y <- joint[, 2]
    wp_x[, 1] <- joint[, 3]
    wp_y[, 1] <- joint[, 4]
  }

  for (t in 2:waves) {
    e <- MASS::mvrnorm(n, c(0, 0), innovation_sigma)
    wp_x[, t] <- beta_x * wp_x[, t - 1] + omega_yx * wp_y[, t - 1] + e[, 1]
    wp_y[, t] <- beta_y * wp_y[, t - 1] + omega_xy * wp_x[, t - 1] + e[, 2]
  }

  # Observed = trait + state
  ox <- wp_x + eta_x
  oy <- wp_y + eta_y

  dat <- data.frame(ox, oy)
  names(dat) <- c(paste0("x", 1:waves), paste0("y", 1:waves))
  dat
}
