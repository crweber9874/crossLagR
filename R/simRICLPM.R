#' @title simRICLPM
#' @description Simulate data from an observed, random intercept cross-lagged regression (RI-CLR) model,
#' including an optional third variable that can interact with X and Y.
#'
#' @param waves The number of waves (time points) in the model.
#' @param stability.p The stability parameter for the X variable (autoregressive effect).
#' @param stability.q The stability parameter for the Y variable (autoregressive effect).
#' @param stability.r The stability parameter for the Z variable (autoregressive effect, if applicable).
#' @param cross.p The cross-lagged effect of Y on X at the next time point.
#' @param cross.q The cross-lagged effect of X on Y at the next time point.
#' @param cross.zx The cross-lagged effect of Z on X at the next time point (if applicable).
#' @param cross.zy The cross-lagged effect of Z on Y at the next time point (if applicable).
#' @param variance.p The variance for the latent X variables.
#' @param variance.q The variance for the latent Y variables.
#' @param variance.r The variance for the latent Z variables (if applicable).
#' @param cov.pq The covariance between X and Y within the same time point.
#' @param cov.pr The covariance between X and Z within the same time point (if applicable).
#' @param cov.qr The covariance between Y and Z within the same time point (if applicable).
#' @param variance.between.x The variance for the random intercept of X.
#' @param variance.between.y The variance for the random intercept of Y.
#' @param variance.between.z The variance for the random intercept of Z (if applicable).
#' @param cov.between.xy The covariance of intercept terms between X and Y.
#' @param cov.between.xz The covariance of intercept terms between X and Z (if applicable).
#' @param cov.between.yz The covariance of intercept terms between Y and Z (if applicable).
#' @param ... Additional arguments to pass to the `lavaan::simulateData` function.
#'
#' @return A list containing two elements:
#'    * `model`: The Lavaan model syntax used for data simulation.
#'    * `data`: The simulated data in a data frame format.
#'
#' @export
#'
simRICLPM <- function(waves = 10,
                      stability.p = 0.2,
                      stability.q = 0.2,
                      stability.r = NULL,
                      cross.p = 0.1,
                      cross.q = 0.1,
                      cross.zx = NULL,
                      cross.zy = NULL,
                      variance.p = 1,
                      variance.q = 1,
                      variance.r = NULL,
                      cov.pq = 0.1,
                      cov.pr = NULL,
                      cov.qr = NULL,
                      variance.between.x = 1,
                      variance.between.y = 1,
                      variance.between.z = NULL,
                      cov.between.xy = 0.5,
                      cov.between.xz = NULL,
                      cov.between.yz = NULL,
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
  params <- list(stability.p = stability.p, stability.q = stability.q, stability.r = stability.r,
                 cross.p = cross.p, cross.q = cross.q, cross.zx = cross.zx, cross.zy = cross.zy,
                 variance.p = variance.p, variance.q = variance.q, variance.r = variance.r,
                 cov.pq = cov.pq, cov.pr = cov.pr, cov.qr = cov.qr,
                 variance.between.x = variance.between.x, variance.between.y = variance.between.y, variance.between.z = variance.between.z,
                 cov.between.xy = cov.between.xy, cov.between.xz = cov.between.xz, cov.between.yz = cov.between.yz)

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

  if (!is.null(stability.r)) {
    model_string <- paste0(model_string, "\nBZ =~ 1*z1")
    for (w in 2:waves) {
      model_string <- paste0(model_string, " + 1*z", w)
    }
  }

  # Indicator means
  for (w in 1:waves) {
    model_string <- paste0(model_string, "\nx", w, "~ 1\n")
    model_string <- paste0(model_string, "\ny", w, "~ 1\n")
    if (!is.null(stability.r)) {
      model_string <- paste0(model_string, "\nz", w, "~ 1\n")
    }
  }

  # Latent variable measurement model
  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\np", w, " =~ 1*x", w,
      "\nq", w, " =~ 1*y", w
    )
    if (!is.null(stability.r)) {
      model_string <- paste0(
        model_string, "\nr", w, " =~ 1*z", w
      )
    }
  }

  # Stability and cross-lagged paths
  for (w in 2:waves) {
    model_string <- paste0(
      model_string, "\n p", w, " ~ ", stability.p, " * p", w - 1, " + ", cross.q, " * q", w - 1
    )
    model_string <- paste0(
      model_string, "\n q", w, " ~ ", stability.q, " * q", w - 1, " + ", cross.p, " * p", w - 1
    )
    if (!is.null(stability.r)) {
      model_string <- paste0(
        model_string, "\n p", w, " ~ ", cross.zx, " * r", w - 1,
        "\n q", w, " ~ ", cross.zy, " * r", w - 1,
        "\n r", w, " ~ ", stability.r, " * r", w - 1,
        " + ", cross.zx, " * p", w - 1, " + ", cross.zy, " * q", w - 1
      )
    }
  }

  # Covariances
  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\n p", w, " ~~ ", variance.p, " * p", w,
      "\n q", w, " ~~ ", variance.q, " * q", w,
      "\n p", w, " ~~ ", cov.pq, " * q", w
    )
    if (!is.null(stability.r)) {
      model_string <- paste0(
        model_string, "\n r", w, " ~~ ", variance.r, " * r", w,
        "\n p", w, " ~~ ", cov.pr, " * r", w,
        "\n q", w, " ~~ ", cov.qr, " * r", w
      )
    }
  }

  # Fix observed residuals to zero
  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\nx", w, " ~~ 0*x", w,
      "\ny", w, " ~~ 0*y", w
    )
    if (!is.null(stability.r)) {
      model_string <- paste0(
        model_string, "\nz", w, " ~~ 0*z", w
      )
    }
  }

  # Random intercept variances and covariances
  model_string <- paste0(
    model_string,
    "\nBX ~~ ", variance.between.x, " * BX",
    "\nBY ~~ ", variance.between.y, " * BY",
    "\nBX ~~ ", cov.between.xy, " * BY"
  )

  if (!is.null(stability.r)) {
    model_string <- paste0(
      model_string,
      "\nBZ ~~ ", variance.between.z, " * BZ",
      "\nBX ~~ ", cov.between.xz, " * BZ",
      "\nBY ~~ ", cov.between.yz, " * BZ"
    )
  }

  dat <- lavaan::simulateData(model = model_string, int.ov.free = TRUE, ...)

  return(list(model = model_string, data = dat))
}
