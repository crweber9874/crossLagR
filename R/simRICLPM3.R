#' @title simRICLPM
#' @description Simulate data from an observed, random intercept cross-lagged regression (RI-CLR) model
#' In particular, this creates a synthetic dataset of a cross-lagged model with a specified number of waves and structural parameters.
#'
#' @param waves The number of waves (time points) in the model.
#' @param stability_p The stability parameter for the x variable (autoregressive effect).
#' @param stability_q The stability parameter for the y variable (autoregressive effect).
#' @param stability_r The stability parameter for the z variable (autoregressive effect).
#' @param cov_pq The covariance between x and y within the same time point.
#' @param cov_pr The covariance between x and z within the same time point.
#' @param cov_qr The covariance between y and z within the same time point.
#' @param cross_q The cross-lagged effect of x on y at the next time point.
#' @param cross_p The cross-lagged effect of y on x at the next time point.
#' @param cross_r The cross-lagged effect of z on x and y at the next time point.
#' @param ... Additional arguments to pass to the `lavaan::simulateData` function.
#'
#' @return A list containing two elements:
#'    * `model`: The Lavaan model syntax used for data simulation.
#'    * `data`:  The simulated data in a data frame format.
#'
#' @export
#'

simRICLPM3 <- function(waves = 10,
                       stability_p = 0.2,
                       stability_q = 0.2,
                       stability_r = 0.2,
                       cross_p = 0.1,
                       cross_q = 0.1,
                       cross_r = 0.1,
                       variance_p = 1,
                       variance_q = 1,
                       variance_r = 1,
                       cov_pq = 0.1,
                       cov_pr = 0.1,
                       cov_qr = 0.1,
                       variance_between_x = 1, # random intercepts, x
                       variance_between_y = 1, # random intercepts, y
                       variance_between_z = 1, # random intercepts, z
                       cov_between_xy = 0.5, # covariance of intercept terms x and y
                       cov_between_xz = 0.5, # covariance of intercept terms x and z
                       cov_between_yz = 0.5, # covariance of intercept terms y and z
                       ...) {
  model_string <- ""
  model_string <- paste0(model_string, "\n BX =~  1* x1")
  for (w in 2:waves) {
    model_string <- paste0(model_string, " + 1 *x", w, "")
  }
  model_string <- paste0(model_string, "\n BY =~   1* y1")
  for (w in 2:waves) {
    model_string <- paste0(model_string, " + 1 * y", w)
  }
  model_string <- paste0(model_string, "\n BZ =~   1* z1")
  for (w in 2:waves) {
    model_string <- paste0(model_string, " + 1 * z", w)
  }
  model_string <- paste0(model_string, "\n")

  for (w in 1:waves) {
    model_string <- paste0(model_string, "x", w, "~ 1", "\n")
  }
  for (w in 1:waves) {
    model_string <- paste0(model_string, "y", w, "~ 1", "\n")
  }
  for (w in 1:waves) {
    model_string <- paste0(model_string, "z", w, "~ 1", "\n")
  }
  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\np", w, " =~ 1*x", w,
      "\nq", w, " =~ 1*y", w,
      "\nr", w, " =~ 1*z", w
    )
  }

  # Stability
  for (w in 2:waves) {
    model_string <- paste0(
      model_string, "\n p", w, " ~ ", stability_p, " * p", w - 1, " + ", cross_q, " * q", w - 1, " + ", cross_r, " * r", w - 1,
      "\n q", w, " ~  ", stability_q, " * q", w - 1, " + ", cross_p, " * p", w - 1, " + ", cross_r, " * r", w - 1,
      "\n r", w, " ~  ", stability_r, " * r", w - 1, " + ", cross_p, " * p", w - 1, " + ", cross_q, " * q", w - 1
    )
  }

  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\n p", w, " ~~ ", variance_p, " * p", w,
      "\n q", w, " ~~ ", variance_q, " * q", w,
      "\n r", w, " ~~ ", variance_r, " * r", w,
      "\n p", w, " ~~ ", cov_pq, " * q", w,
      "\n p", w, " ~~ ", cov_pr, " * r", w,
      "\n q", w, " ~~ ", cov_qr, " * r", w
    )
  }

  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\nx", w, "~~",
      "0*", "x", w
    )
  }
  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\ny", w, "~~",
      "0*", "y", w
    )
  }
  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\nz", w, "~~",
      "0*", "z", w
    )
  }

  model_string <- paste0(
    model_string,
    "\n BX ~~", variance_between_x, "* BX",
    "\n BY ~~", variance_between_y, "* BY",
    "\n BZ ~~", variance_between_z, "* BZ",
    "\n BX ~~", cov_between_xy, "* BY",
    "\n BX ~~", cov_between_xz, "* BZ",
    "\n BY ~~", cov_between_yz, "* BZ"
  )

  dat <- lavaan::simulateData(
    model = model_string,
    int.ov.free = FALSE
  )

  return(list(model = model_string, data = dat))
}
