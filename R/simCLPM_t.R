#' @title simCLPM_t
#' @description Simulate data from a Cross-Lagged Panel Model that includes
#'   contemporaneous causal effects between x and y within each wave (in
#'   addition to the usual lagged effects).
#'
#' @details
#' The structural model at wave \eqn{t > 1} is
#' \deqn{x_t = \beta_x x_{t-1} + \omega_{yx} y_{t-1} + \gamma_{yx} y_t + \epsilon_{x,t}}
#' \deqn{y_t = \beta_y y_{t-1} + \omega_{xy} x_{t-1} + \gamma_{xy} x_t + \epsilon_{y,t}}
#' where \eqn{\gamma_{xy}} (\code{contemp_xy}) is the simultaneous effect of
#' \eqn{x_t} on \eqn{y_t} and \eqn{\gamma_{yx}} (\code{contemp_yx}) is the
#' simultaneous effect of \eqn{y_t} on \eqn{x_t}. Stacking
#' \eqn{z_t = (x_t, y_t)'}, this is
#' \deqn{z_t = \Gamma z_t + B z_{t-1} + \epsilon_t,}
#' which has reduced form
#' \deqn{z_t = (I - \Gamma)^{-1} B z_{t-1} + (I - \Gamma)^{-1} \epsilon_t.}
#' The function simulates from this reduced form directly. Identification of
#' the simultaneous system requires \eqn{|I - \Gamma| \neq 0}, equivalently
#' \eqn{\gamma_{xy} \gamma_{yx} \neq 1}, and stationarity of the reduced-form
#' AR matrix \eqn{M = (I - \Gamma)^{-1} B}.
#'
#' @param waves Number of waves (\eqn{\geq 2}).
#' @param beta_x Autoregressive effect for x.
#' @param beta_y Autoregressive effect for y.
#' @param omega_xy Lagged effect \eqn{x_{t-1} \to y_t}.
#' @param omega_yx Lagged effect \eqn{y_{t-1} \to x_t}.
#' @param contemp_xy Contemporaneous effect \eqn{x_t \to y_t}. Default 0.
#' @param contemp_yx Contemporaneous effect \eqn{y_t \to x_t}. Default 0.
#' @param var_x,var_y Innovation variances of the structural shocks
#'   \eqn{\epsilon_{x,t}}, \eqn{\epsilon_{y,t}} (\eqn{t > 1}).
#' @param cov_xy Innovation covariance between \eqn{\epsilon_{x,t}} and
#'   \eqn{\epsilon_{y,t}}. Default 0 -- under non-zero contemporaneous
#'   effects the within-wave observed covariance is mostly carried by
#'   \eqn{\Gamma}, so leaving residual covariance at 0 is the usual choice.
#' @param var_x1,var_y1,cov_xy1 First-wave covariance structure. Default
#'   to \code{var_x}, \code{var_y}, \code{cov_xy}.
#' @param mean_x1,mean_y1 First-wave means. Default 0.
#' @param sample_size Number of persons. Default 1000.
#' @param seed Optional integer seed.
#'
#' @return List with elements:
#'   \itemize{
#'     \item \code{model}: human-readable description of the structural model.
#'     \item \code{data}: simulated data, wide format with columns
#'       \code{x1, y1, x2, y2, ...}.
#'     \item \code{parameters}: the parameter values used.
#'     \item \code{reduced_form}: the reduced-form AR matrix \code{M} and
#'       reduced-form innovation covariance \code{Sigma_reduced}.
#'   }
#' @examples
#' sim <- simCLPM_t(waves = 5, beta_x = 0.3, beta_y = 0.3,
#'                  omega_xy = 0.05, omega_yx = 0.05,
#'                  contemp_xy = 0.10, contemp_yx = 0.15,
#'                  sample_size = 500, seed = 1)
#' head(sim$data)
#' sim$reduced_form$M
#' @export
simCLPM_t <- function(waves = 5,
                      beta_x = 0.30,
                      beta_y = 0.30,
                      omega_xy = 0.05,
                      omega_yx = 0.05,
                      contemp_xy = 0,
                      contemp_yx = 0,
                      var_x = 1,
                      var_y = 1,
                      cov_xy = 0,
                      var_x1 = NULL,
                      var_y1 = NULL,
                      cov_xy1 = NULL,
                      mean_x1 = 0,
                      mean_y1 = 0,
                      sample_size = 1000,
                      seed = NULL) {

  if (!is.numeric(waves) || waves < 2 || waves != as.integer(waves)) {
    stop("`waves` must be an integer >= 2.", call. = FALSE)
  }
  scalars <- list(beta_x = beta_x, beta_y = beta_y,
                  omega_xy = omega_xy, omega_yx = omega_yx,
                  contemp_xy = contemp_xy, contemp_yx = contemp_yx,
                  var_x = var_x, var_y = var_y, cov_xy = cov_xy,
                  mean_x1 = mean_x1, mean_y1 = mean_y1,
                  sample_size = sample_size)
  for (nm in names(scalars)) {
    if (!is.numeric(scalars[[nm]]) || length(scalars[[nm]]) != 1L) {
      stop(sprintf("`%s` must be a single numeric value.", nm), call. = FALSE)
    }
  }
  if (var_x <= 0 || var_y <= 0) {
    stop("`var_x` and `var_y` must be positive.", call. = FALSE)
  }

  if (is.null(var_x1))  var_x1  <- var_x
  if (is.null(var_y1))  var_y1  <- var_y
  if (is.null(cov_xy1)) cov_xy1 <- cov_xy

  # Structural matrices: z = (x, y)'
  Gamma <- matrix(c(0,          contemp_yx,
                    contemp_xy, 0),         nrow = 2, byrow = TRUE)
  Bmat  <- matrix(c(beta_x,   omega_yx,
                    omega_xy, beta_y),      nrow = 2, byrow = TRUE)
  Sigma_eps <- matrix(c(var_x, cov_xy,
                        cov_xy, var_y),     nrow = 2)

  # Identification: (I - Gamma) must be invertible
  det_check <- 1 - contemp_xy * contemp_yx
  if (abs(det_check) < .Machine$double.eps^0.5) {
    stop("(I - Gamma) is singular: contemp_xy * contemp_yx ~= 1.",
         call. = FALSE)
  }
  I2 <- diag(2)
  IGinv <- solve(I2 - Gamma)
  M     <- IGinv %*% Bmat
  Sigma_reduced <- IGinv %*% Sigma_eps %*% t(IGinv)

  # Stationarity check on reduced-form M
  eig <- abs(eigen(M, only.values = TRUE)$values)
  if (any(eig >= 1)) {
    warning(sprintf(
      "Reduced-form AR matrix is non-stationary (max |eigenvalue| = %.3f).",
      max(eig)), call. = FALSE)
  }

  # Validate first-wave covariance is PD
  Sigma_1 <- matrix(c(var_x1, cov_xy1, cov_xy1, var_y1), nrow = 2)
  if (any(eigen(Sigma_1, only.values = TRUE)$values <= 0)) {
    stop("First-wave covariance matrix is not positive definite.",
         call. = FALSE)
  }
  if (any(eigen(Sigma_reduced, only.values = TRUE)$values <= 0)) {
    stop("Reduced-form innovation covariance is not positive definite.",
         call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)

  L1 <- chol(Sigma_1)
  Lr <- chol(Sigma_reduced)
  mu1 <- c(mean_x1, mean_y1)

  N <- sample_size
  z <- array(NA_real_, dim = c(N, waves, 2))
  z[, 1, ] <- matrix(stats::rnorm(2 * N), N, 2) %*% L1 +
    matrix(mu1, N, 2, byrow = TRUE)
  for (t in 2:waves) {
    eps_t <- matrix(stats::rnorm(2 * N), N, 2) %*% Lr
    z[, t, ] <- z[, t - 1, ] %*% t(M) + eps_t
  }

  # Wide format: x1, y1, x2, y2, ...
  out_data <- data.frame(matrix(NA_real_, nrow = N, ncol = 2 * waves))
  col_names <- character(2 * waves)
  for (t in seq_len(waves)) {
    out_data[, 2 * t - 1] <- z[, t, 1]
    out_data[, 2 * t]     <- z[, t, 2]
    col_names[2 * t - 1]  <- paste0("x", t)
    col_names[2 * t]      <- paste0("y", t)
  }
  colnames(out_data) <- col_names

  model_string <- paste0(
    "# simCLPM_t structural model (non-recursive within-wave)\n",
    "# z_t = Gamma * z_t + B * z_{t-1} + eps_t,  z = (x, y)'\n",
    sprintf("# Gamma = [[0, %.4g], [%.4g, 0]]\n", contemp_yx, contemp_xy),
    sprintf("# B     = [[%.4g, %.4g], [%.4g, %.4g]]\n",
            beta_x, omega_yx, omega_xy, beta_y),
    sprintf("# Sigma_eps = [[%.4g, %.4g], [%.4g, %.4g]]\n",
            var_x, cov_xy, cov_xy, var_y),
    sprintf("# Wave-1 cov = [[%.4g, %.4g], [%.4g, %.4g]],  means = (%.4g, %.4g)\n",
            var_x1, cov_xy1, cov_xy1, var_y1, mean_x1, mean_y1)
  )

  list(
    model = model_string,
    data  = out_data,
    parameters = list(
      waves = waves, beta_x = beta_x, beta_y = beta_y,
      omega_xy = omega_xy, omega_yx = omega_yx,
      contemp_xy = contemp_xy, contemp_yx = contemp_yx,
      var_x = var_x, var_y = var_y, cov_xy = cov_xy,
      var_x1 = var_x1, var_y1 = var_y1, cov_xy1 = cov_xy1,
      mean_x1 = mean_x1, mean_y1 = mean_y1,
      sample_size = sample_size
    ),
    reduced_form = list(M = M, Sigma_reduced = Sigma_reduced)
  )
}
