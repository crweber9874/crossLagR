#' @title simCLPM
#' @description Simulate data from a cross-lagged panel model (CLPM)
#' In particular, this creates a synthetic dataset of a cross-lagged model with a specified number of waves and structural parameters.
#'
#' @param waves The number of waves (time points) in the model.
#' @param stability.p The stability parameter for the x variable (autoregressive effect).
#' @param stability.q The stability parameter for the y variable (autoregressive effect).
#' @param cov.pq The covariance between x and y within the same time point.
#' @param cross.q The cross-lagged effect of x on y at the next time point.
#' @param cross.p The cross-lagged effect of y on x at the next time point.
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
                            stability.p = 0.2,
                            stability.q = 0.2,
                            cross.p = 0.1,
                            cross.q = 0.1,
                            variance.p = 1,
                            variance.q = 1,
                            cov.pq = 0.1,
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
      model_string, "\n p", w, " ~ ", stability.p, " * p", w - 1, " + ", cross.q, " * q", w - 1,
      "\n q", w, " ~  ", stability.q, " * q", w - 1, " + ", cross.p, " * p", w - 1
    )
  }

  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\n p", w, " ~~ ", variance.p, " * p", w,
      "\n q", w, " ~~ ", variance.q, " * q", w,
      "\n p", w, " ~~ ", cov.pq, " * q", w
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
