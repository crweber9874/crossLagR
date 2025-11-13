#' Simulate Data from a Random Intercept Cross-Lagged Panel Model (RI-CLPM)
#'
#' Generates simulated longitudinal data following a Random Intercept Cross-Lagged
#' Panel Model structure. The model includes stable between-person traits (random
#' intercepts) and within-person dynamics (autoregressive and cross-lagged effects).
#'
#' @param n Integer. Number of individuals/units to simulate. Default is 500.
#' @param waves Integer. Number of time points/waves. Default is 4.
#' @param ar_x Numeric. Autoregressive effect for X (stability of X over time).
#'   Default is 0.5.
#' @param ar_y Numeric. Autoregressive effect for Y (stability of Y over time).
#'   Default is 0.5.
#' @param cl_xy Numeric. Cross-lagged effect from X to Y (X_{t-1} -> Y_t).
#'   Default is 0.3.
#' @param cl_yx Numeric. Cross-lagged effect from Y to X (Y_{t-1} -> X_t).
#'   Default is 0.2.
#' @param cor_traits Numeric. Correlation between stable traits (random intercepts).
#'   Range: [-1, 1]. Default is 0.5.
#' @param innov_var Numeric. Variance of wave-specific innovations (within-person
#'   deviations). Default is 0.5.
#' @param cor_xy Numeric. Within-time correlation between X and Y innovations.
#'   Range: [-1, 1]. Default is 0.3.
#' @param seed Integer. Random seed for reproducibility. If NULL, no seed is set.
#'   Default is NULL.
#'
#' @return A tibble in wide format with columns:
#'   \item{id}{Individual identifier (1 to n)}
#'   \item{traitx}{Stable between-person trait for X}
#'   \item{traity}{Stable between-person trait for Y}
#'   \item{x1, x2, ..., xW}{Observed values of X at each wave}
#'   \item{y1, y2, ..., yW}{Observed values of Y at each wave}
#'
#' @details
#' The RI-CLPM separates between-person and within-person variance:
#'
#' **Between-person level:**
#' - Stable traits are drawn from a bivariate normal distribution with correlation
#'   \code{cor_traits}
#'
#' **Within-person level (for t >= 2):**
#' \deqn{X_t = \mu_x + ar_x(X_{t-1} - \mu_x) + cl_yx(Y_{t-1} - \mu_y) + \epsilon_{x,t}}
#' \deqn{Y_t = \mu_y + ar_y(Y_{t-1} - \mu_y) + cl_xy(X_{t-1} - \mu_x) + \epsilon_{y,t}}
#'
#' where \eqn{\mu_x} and \eqn{\mu_y} are individual-specific traits (random intercepts),
#' and innovations have \code{Cor(\epsilon_{x,t}, \epsilon_{y,t}) = cor_xy}.
#'
#' @examples
#' # Basic simulation with defaults
#' data <- baseRICLPM_sim()
#'
#' # Custom parameters
#' data <- baseRICLPM_sim(
#'   n = 1000,
#'   waves = 5,
#'   ar_x = 0.6,
#'   ar_y = 0.4,
#'   cl_xy = 0.2,
#'   cl_yx = 0.15,
#'   cor_traits = 0.4,
#'   seed = 42
#' )
#'
#' # Strong cross-lagged effects, weak stability
#' data <- baseRICLPM_sim(
#'   ar_x = 0.3,
#'   ar_y = 0.3,
#'   cl_xy = 0.5,
#'   cl_yx = 0.5,
#'   seed = 123
#' )
#'
#' @export
baseRICLPM_sim <- function(n = 500,
                           waves = 4,
                           ar_x = 0.5,
                           ar_y = 0.5,
                           cl_xy = 0.3,
                           cl_yx = 0.2,
                           cor_traits = 0.5,
                           innov_var = 0.5,
                           cor_xy = 0.3,
                           seed = NULL) {

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Generate correlated traits (random intercepts)
  trait_mean <- c(0, 0)
  trait_cov <- matrix(c(1, cor_traits,
                        cor_traits, 1),
                      nrow = 2, byrow = TRUE)

  traits <- MASS::mvrnorm(n = n,
                          mu = trait_mean,
                          Sigma = trait_cov)
  trait_x <- traits[, 1]
  trait_y <- traits[, 2]

  # Initialize data matrices
  X <- matrix(NA, nrow = n, ncol = waves)
  Y <- matrix(NA, nrow = n, ncol = waves)

  # Wave 1: Traits + correlated innovations
  innov_cov <- matrix(c(innov_var, cor_xy * innov_var,
                        cor_xy * innov_var, innov_var), 2, 2)
  wave1_innov <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = innov_cov)

  X[, 1] <- trait_x + wave1_innov[, 1]
  Y[, 1] <- trait_y + wave1_innov[, 2]

  # Waves 2+: RI-CLPM dynamics with correlated innovations
  for (t in 2:waves) {
    innovations <- MASS::mvrnorm(n, c(0, 0),
                                 matrix(c(innov_var, cor_xy * innov_var,
                                          cor_xy * innov_var, innov_var), 2, 2))

    # Within-person deviations from trait
    X[, t] <- trait_x + ar_x * (X[, t-1] - trait_x) +
      cl_yx * (Y[, t-1] - trait_y) + innovations[, 1]
    Y[, t] <- trait_y + ar_y * (Y[, t-1] - trait_y) +
      cl_xy * (X[, t-1] - trait_x) + innovations[, 2]
  }

  # Convert to wide-format tibble
  data <- tibble::tibble(
    id = rep(1:n, each = waves),
    wave = rep(1:waves, times = n)
  ) |>
    dplyr::mutate(
      x = c(t(X)),
      y = c(t(Y)),
      trait_x = rep(trait_x, each = waves),
      trait_y = rep(trait_y, each = waves)
    ) |>
    tidyr::pivot_wider(
      names_from = wave,
      values_from = c(x, y),
      names_prefix = ""
    ) |>
    dplyr::rename_with(~ stringr::str_replace_all(.x, "_", ""), everything())

  return(data)
}
