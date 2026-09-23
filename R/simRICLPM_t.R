#' @title simRICLPM_t
#' @description Simulate data from a Random Intercept CLPM that includes
#'   contemporaneous causal effects between x and y within each wave (in
#'   addition to the usual lagged effects and stable random intercepts).
#'
#' @details
#' Observed variables decompose as
#' \deqn{x_{it} = I_{x,i} + \xi_{x,it}, \quad y_{it} = I_{y,i} + \xi_{y,it}}
#' with stable random intercepts \eqn{(I_{x,i}, I_{y,i}) \sim N(\mu_B, \Sigma_B)}.
#' The within-person process \eqn{\xi_t} follows the same non-recursive
#' simultaneous-equation structure as \code{\link{simCLPM_t}}:
#' \deqn{\xi_t = \Gamma \xi_t + B \xi_{t-1} + \epsilon_t,}
#' so the reduced form is \eqn{\xi_t = (I-\Gamma)^{-1} B \xi_{t-1} + (I-\Gamma)^{-1} \epsilon_t}.
#' Identification requires \eqn{|I - \Gamma| \neq 0}, and stationarity of
#' the reduced-form within AR matrix.
#'
#' @param waves Number of waves (\eqn{\geq 2}).
#' @param beta_x,beta_y Within-person autoregressions.
#' @param omega_xy,omega_yx Within-person lagged cross-lags.
#' @param contemp_xy,contemp_yx Contemporaneous effects on the
#'   within-person process.
#' @param var_x,var_y,cov_xy Within-person innovation (co)variances.
#' @param var_BX,var_BY,cov_BXBY Random-intercept (co)variances.
#' @param mean_BX,mean_BY Random-intercept means.
#' @param sample_size Number of persons. Default 1000.
#' @param seed Optional integer seed.
#'
#' @return List with elements:
#'   \itemize{
#'     \item \code{model}: human-readable description of the structural model.
#'     \item \code{data}: simulated data, wide format with columns
#'       \code{x1, y1, x2, y2, ...}.
#'     \item \code{parameters}: the parameter values used.
#'     \item \code{reduced_form}: the reduced-form within AR matrix \code{M},
#'       reduced-form within innovation covariance \code{Sigma_reduced}, and
#'       random-intercept covariance \code{Sigma_B}.
#'   }
#' @examples
#' sim <- simRICLPM_t(waves = 5, beta_x = 0.3, beta_y = 0.3,
#'                    omega_xy = 0.05, omega_yx = 0.05,
#'                    contemp_xy = 0.10, contemp_yx = 0.15,
#'                    var_BX = 1, var_BY = 1, cov_BXBY = 0.5,
#'                    sample_size = 500, seed = 1)
#' head(sim$data)
#' @export
simRICLPM_t <- function(waves = 5,
                        beta_x = 0.30,
                        beta_y = 0.30,
                        omega_xy = 0.05,
                        omega_yx = 0.05,
                        contemp_xy = 0,
                        contemp_yx = 0,
                        var_x = 1,
                        var_y = 1,
                        cov_xy = 0,
                        var_BX = 1,
                        var_BY = 1,
                        cov_BXBY = 0.5,
                        mean_BX = 0,
                        mean_BY = 0,
                        sample_size = 1000,
                        seed = NULL) {

  if (!is.numeric(waves) || waves < 2 || waves != as.integer(waves)) {
    stop("`waves` must be an integer >= 2.", call. = FALSE)
  }
  scalars <- list(beta_x = beta_x, beta_y = beta_y,
                  omega_xy = omega_xy, omega_yx = omega_yx,
                  contemp_xy = contemp_xy, contemp_yx = contemp_yx,
                  var_x = var_x, var_y = var_y, cov_xy = cov_xy,
                  var_BX = var_BX, var_BY = var_BY, cov_BXBY = cov_BXBY,
                  mean_BX = mean_BX, mean_BY = mean_BY,
                  sample_size = sample_size)
  for (nm in names(scalars)) {
    if (!is.numeric(scalars[[nm]]) || length(scalars[[nm]]) != 1L) {
      stop(sprintf("`%s` must be a single numeric value.", nm), call. = FALSE)
    }
  }
  if (var_x <= 0 || var_y <= 0 || var_BX <= 0 || var_BY <= 0) {
    stop("Variance parameters must be positive.", call. = FALSE)
  }

  # Structural matrices
  Gamma <- matrix(c(0,          contemp_yx,
                    contemp_xy, 0),         nrow = 2, byrow = TRUE)
  Bmat  <- matrix(c(beta_x,   omega_yx,
                    omega_xy, beta_y),      nrow = 2, byrow = TRUE)
  Sigma_eps <- matrix(c(var_x, cov_xy,
                        cov_xy, var_y),     nrow = 2)
  Sigma_B   <- matrix(c(var_BX, cov_BXBY,
                        cov_BXBY, var_BY),  nrow = 2)

  if (abs(1 - contemp_xy * contemp_yx) < .Machine$double.eps^0.5) {
    stop("(I - Gamma) is singular: contemp_xy * contemp_yx ~= 1.",
         call. = FALSE)
  }
  I2 <- diag(2)
  IGinv <- solve(I2 - Gamma)
  M     <- IGinv %*% Bmat
  Sigma_reduced <- IGinv %*% Sigma_eps %*% t(IGinv)

  eig <- abs(eigen(M, only.values = TRUE)$values)
  if (any(eig >= 1)) {
    warning(sprintf(
      "Reduced-form within AR matrix is non-stationary (max |eigenvalue| = %.3f).",
      max(eig)), call. = FALSE)
  }
  if (any(eigen(Sigma_reduced, only.values = TRUE)$values <= 0)) {
    stop("Reduced-form within innovation covariance is not positive definite.",
         call. = FALSE)
  }
  if (any(eigen(Sigma_B, only.values = TRUE)$values <= 0)) {
    stop("Random-intercept covariance is not positive definite.",
         call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)

  N <- sample_size

  # Random intercepts
  LB <- chol(Sigma_B)
  trait <- matrix(stats::rnorm(2 * N), N, 2) %*% LB +
    matrix(c(mean_BX, mean_BY), N, 2, byrow = TRUE)

  # Within-person process: initialise at stationary distribution of M
  # vec(Psi) = (I - M %x% M)^-1 vec(Sigma_reduced)
  Mkron <- kronecker(M, M)
  vec_Sig <- as.vector(Sigma_reduced)
  vec_Psi <- solve(diag(4) - Mkron, vec_Sig)
  Psi_stationary <- matrix(vec_Psi, 2, 2)
  # Symmetrise (numerical safety)
  Psi_stationary <- (Psi_stationary + t(Psi_stationary)) / 2
  if (any(eigen(Psi_stationary, only.values = TRUE)$values <= 0)) {
    # Fallback: use Sigma_reduced as initial dispersion
    Psi_stationary <- Sigma_reduced
  }

  L1 <- chol(Psi_stationary)
  Lr <- chol(Sigma_reduced)

  xi <- array(NA_real_, dim = c(N, waves, 2))
  xi[, 1, ] <- matrix(stats::rnorm(2 * N), N, 2) %*% L1
  for (t in 2:waves) {
    eps_t <- matrix(stats::rnorm(2 * N), N, 2) %*% Lr
    xi[, t, ] <- xi[, t - 1, ] %*% t(M) + eps_t
  }

  # Observed = trait + within
  z <- xi
  for (t in seq_len(waves)) z[, t, ] <- z[, t, ] + trait

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
    "# simRICLPM_t structural model\n",
    "# x_it = I_xi + xi_xit,  y_it = I_yi + xi_yit\n",
    "# (I_x, I_y)' ~ N(mu_B, Sigma_B)\n",
    "# xi_t = Gamma * xi_t + B * xi_{t-1} + eps_t (non-recursive within-wave)\n",
    sprintf("# Gamma = [[0, %.4g], [%.4g, 0]]\n", contemp_yx, contemp_xy),
    sprintf("# B     = [[%.4g, %.4g], [%.4g, %.4g]]\n",
            beta_x, omega_yx, omega_xy, beta_y),
    sprintf("# Sigma_eps = [[%.4g, %.4g], [%.4g, %.4g]]\n",
            var_x, cov_xy, cov_xy, var_y),
    sprintf("# Sigma_B   = [[%.4g, %.4g], [%.4g, %.4g]],  mu_B = (%.4g, %.4g)\n",
            var_BX, cov_BXBY, cov_BXBY, var_BY, mean_BX, mean_BY)
  )

  list(
    model = model_string,
    data  = out_data,
    parameters = list(
      waves = waves, beta_x = beta_x, beta_y = beta_y,
      omega_xy = omega_xy, omega_yx = omega_yx,
      contemp_xy = contemp_xy, contemp_yx = contemp_yx,
      var_x = var_x, var_y = var_y, cov_xy = cov_xy,
      var_BX = var_BX, var_BY = var_BY, cov_BXBY = cov_BXBY,
      mean_BX = mean_BX, mean_BY = mean_BY,
      sample_size = sample_size
    ),
    reduced_form = list(M = M,
                        Sigma_reduced = Sigma_reduced,
                        Sigma_B       = Sigma_B,
                        Psi_within_stationary = Psi_stationary)
  )
}
