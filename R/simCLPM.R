#' @title simCLPM
#' @description Simulate data from a cross-lagged panel model (CLPM)
#' In particular, this creates a synthetic dataset of a cross-lagged model with a specified number of waves and structural parameters.
#'
#' @param waves The number of waves (time points) in the model.
#' @param stability_p The stability parameter for the x variable (autoregressive effect).
#' @param stability_q The stability parameter for the y variable (autoregressive effect).
#' @param cov_pq The covariance between x and y within the same time point.
#' @param cross_q The cross-lagged effect of x on y at the next time point.
#' @param cross_p The cross-lagged effect of y on x at the next time point.
#' @param variance_p The variance parameter for the x variable.
#' @param variance_q The variance parameter for the y variable.
#' @param ... Additional arguments to pass to the `lavaan::simulateData` function.
#'
#' @return A list containing two elements:
#'    * `model`: The Lavaan model syntax used for data simulation.
#'    * `data`:  The simulated data in a data frame format.
#'
#'
#' @export
#'
#'
simCLPM <- function(waves = 10,
                    stability_p = 0.2,
                    stability_q = 0.2,
                    cross_p = 0.1,
                    cross_q = 0.1,
                    variance_p = 1,
                    variance_q = 1,
                    cov_pq = 0.1,
                    ...) {
  model_string <- ""
  for (w in 1:waves) {
    model_string <- paste0(model_string, "x", w, "~ 1", "\n")
  }
  for (w in 1:waves) {
    model_string <- paste0(model_string, "y", w, "~ 1", "\n")
  }
  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\np", w, " =~ 1*x", w,
      "\nq", w, " =~ 1*y", w
    )
  }
  # Stability
  for (w in 2:waves) {
    model_string <- paste0(
      model_string, "\n p", w, " ~ ", stability_p, " * p", w - 1, " + ", cross_q, " * q", w - 1,
      "\n q", w, " ~  ", stability_q, " * q", w - 1, " + ", cross_p, " * p", w - 1
    )
  }
  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\n p", w, " ~~ ", variance_p, " * p", w,
      "\n q", w, " ~~ ", variance_q, " * q", w,
      "\n p", w, " ~~ ", cov_pq, " * q", w
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
  dat <- lavaan::simulateData(model = model_string,
                              int.ov.free = TRUE)
  return(list(model = model_string, data = dat))
}
