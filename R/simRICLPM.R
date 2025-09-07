#' @title simRICLPM
#' @description Simulate data from an observed, random intercept cross-lagged regression (RI-CLR) model,
#' including an optional third variable that can interact with X and Y.
#'
#' @param waves The number of waves (time points) in the model.
#' @param stability_p The stability parameter for the X variable (autoregressive effect).
#' @param stability_q The stability parameter for the Y variable (autoregressive effect).
#' @param stability_r The stability parameter for the Z variable (autoregressive effect, if applicable).
#' @param cross_p The cross-lagged effect of Y on X at the next time point.
#' @param cross_q The cross-lagged effect of X on Y at the next time point.
#' @param cross_zx The cross-lagged effect of Z on X at the next time point (if applicable).
#' @param cross_zy The cross-lagged effect of Z on Y at the next time point (if applicable).
#' @param variance_p The variance for the latent X variables.
#' @param variance_q The variance for the latent Y variables.
#' @param variance_r The variance for the latent Z variables (if applicable).
#' @param cov_pq The covariance between X and Y within the same time point.
#' @param cov_pr The covariance between X and Z within the same time point (if applicable).
#' @param cov_qr The covariance between Y and Z within the same time point (if applicable).
#' @param variance_between_x The variance for the random intercept of X.
#' @param variance_between_y The variance for the random intercept of Y.
#' @param variance_between_z The variance for the random intercept of Z (if applicable).
#' @param cov_between_xy The covariance of intercept terms between X and Y.
#' @param cov_between_xz The covariance of intercept terms between X and Z (if applicable).
#' @param cov_between_yz The covariance of intercept terms between Y and Z (if applicable).
#' @param ... Additional arguments to pass to the `lavaan::simulateData` function.
#'
#' @return A list containing two elements:
#'    * `model`: The Lavaan model syntax used for data simulation.
#'    * `data`: The simulated data in a data frame format.
#'
#' @export
#'
simRICLPM <- function(waves = 10,
                      stability_p = 0.2,
                      stability_q = 0.2,
                      stability_r = NULL,
                      cross_p = 0.1,
                      cross_q = 0.1,
                      cross_zx = NULL,
                      cross_zy = NULL,
                      variance_p = 1,
                      variance_q = 1,
                      variance_r = NULL,
                      cov_pq = 0.1,
                      cov_pr = NULL,
                      cov_qr = NULL,
                      variance_between_x = 1,
                      variance_between_y = 1,
                      variance_between_z = NULL,
                      cov_between_xy = 0.5,
                      cov_between_xz = NULL,
                      cov_between_yz = NULL,
                      ...) {

  # Input validation
  if (!is.numeric(waves) || waves <= 0 || waves != as.integer(waves)) {
    stop("❌ Error: 'waves' must be a positive integer.")
  }

  check_numeric <- function(param, name) {
    if (!is.null(param) && !is.numeric(param)) {
      stop(paste0("❌ Error: '", name, "' must be numeric."))
    }
  }

  # Validate numeric parameters
  params <- list(stability_p = stability_p, stability_q = stability_q, stability_r = stability_r,
                 cross_p = cross_p, cross_q = cross_q, cross_zx = cross_zx, cross_zy = cross_zy,
                 variance_p = variance_p, variance_q = variance_q, variance_r = variance_r,
                 cov_pq = cov_pq, cov_pr = cov_pr, cov_qr = cov_qr,
                 variance_between_x = variance_between_x, variance_between_y = variance_between_y, variance_between_z = variance_between_z,
                 cov_between_xy = cov_between_xy, cov_between_xz = cov_between_xz, cov_between_yz = cov_between_yz)

  lapply(names(params), function(name) {
    check_numeric(params[[name]], name)
  })

  model_string <- ""
  model_string <- paste0(model_string, "\nBX =~ 1*x1")
  for (w in 2:waves) {
    model_string <- paste0(model_string, " + 1*x", w)
  }
  model_string <- paste0(model_string, "\nBY =~ 1*y1")
  for (w in 2:waves) {
    model_string <- paste0(model_string, " + 1*y", w)
  }

  if (!is.null(stability_r)) {
    model_string <- paste0(model_string, "\nBZ =~ 1*z1")
    for (w in 2:waves) {
      model_string <- paste0(model_string, " + 1*z", w)
    }
  }

  # Indicator means
  for (w in 1:waves) {
    model_string <- paste0(model_string, "\nx", w, "~ 1\n")
    model_string <- paste0(model_string, "\ny", w, "~ 1\n")
    if (!is.null(stability_r)) {
      model_string <- paste0(model_string, "\nz", w, "~ 1\n")
    }
  }

  # Latent variable measurement model
  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\np", w, " =~ 1*x", w,
      "\nq", w, " =~ 1*y", w
    )
    if (!is.null(stability_r)) {
      model_string <- paste0(
        model_string, "\nr", w, " =~ 1*z", w
      )
    }
  }

  # Stability and cross-lagged paths
  for (w in 2:waves) {
    model_string <- paste0(
      model_string, "\n p", w, " ~ ", stability_p, " * p", w - 1, " + ", cross_q, " * q", w - 1
    )
    model_string <- paste0(
      model_string, "\n q", w, " ~ ", stability_q, " * q", w - 1, " + ", cross_p, " * p", w - 1
    )
    if (!is.null(stability_r)) {
      model_string <- paste0(
        model_string, "\n p", w, " ~ ", cross_zx, " * r", w - 1,
        "\n q", w, " ~ ", cross_zy, " * r", w - 1,
        "\n r", w, " ~ ", stability_r, " * r", w - 1,
        " + ", cross_zx, " * p", w - 1, " + ", cross_zy, " * q", w - 1
      )
    }
  }

  # Covariances
  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\n p", w, " ~~ ", variance_p, " * p", w,
      "\n q", w, " ~~ ", variance_q, " * q", w,
      "\n p", w, " ~~ ", cov_pq, " * q", w
    )
    if (!is.null(stability_r)) {
      model_string <- paste0(
        model_string, "\n r", w, " ~~ ", variance_r, " * r", w,
        "\n p", w, " ~~ ", cov_pr, " * r", w,
        "\n q", w, " ~~ ", cov_qr, " * r", w
      )
    }
  }

  # Fix observed residuals to zero
  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\nx", w, " ~~ 0*x", w,
      "\ny", w, " ~~ 0*y", w
    )
    if (!is.null(stability_r)) {
      model_string <- paste0(
        model_string, "\nz", w, " ~~ 0*z", w
      )
    }
  }

  # Random intercept variances and covariances
  model_string <- paste0(
    model_string,
    "\nBX ~~ ", variance_between_x, " * BX",
    "\nBY ~~ ", variance_between_y, " * BY",
    "\nBX ~~ ", cov_between_xy, " * BY"
  )

  if (!is.null(stability_r)) {
    model_string <- paste0(
      model_string,
      "\nBZ ~~ ", variance_between_z, " * BZ",
      "\nBX ~~ ", cov_between_xz, " * BZ",
      "\nBY ~~ ", cov_between_yz, " * BZ"
    )
  }

  dat <- lavaan::simulateData(model = model_string, int.ov.free = TRUE, ...)

  return(list(model = model_string, data = dat))
}
