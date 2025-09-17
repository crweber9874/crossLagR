
#' @title simCLPM_timeInvariantU
#' @description Simulate data from a cross-lagged panel model (CLPM) with a time-invariant confounder.
#'
#' @param waves The number of waves (time points) in the model.
#' @param stability_p The stability parameter for the x variable (autoregressive effect).
#' @param stability_q The stability parameter for the y variable (autoregressive effect).
#' @param cov_pq The covariance between x and y within the same time point.
#' @param cross_q The cross-lagged effect of x on y at the next time point.
#' @param cross_p The cross-lagged effect of y on x at the next time point.
#' @param variance_p The variance of the p latent variable.
#' @param variance_q The variance of the q latent variable.
#' @param confounder_p The effect of the time-invariant confounder on x variables.
#' @param confounder_q The effect of the time-invariant confounder on y variables.
#' @param confounder_variance The variance of the time-invariant confounder.
#' @param ... Additional arguments to pass to the `lavaan::simulateData` function.
#'
#' @return A list containing two elements:
#'    * `model`: The Lavaan model syntax used for data simulation.
#'    * `data`:  The simulated data in a data frame format.
#'
#' @export
simCLPM_timeInvariantU <- function(waves = 10,
                                   stability_p = 0.2,
                                   stability_q = 0.2,
                                   cross_p = 0.1,
                                   cross_q = 0.1,
                                   variance_p = 1,
                                   variance_q = 1,
                                   cov_pq = 0.1,
                                   confounder_p = 0.3,
                                   confounder_q = 0.3,
                                   confounder_variance = 1,
                                   ...) {

  model_string <- ""

  # Define latent time-invariant confounder
  model_string <- paste0(model_string, "\nU =~ ")
  for (w in 1:waves) {
    model_string <- paste0(model_string, "1*x", w, " + 1*y", w, " + ")
  }
  # Remove trailing " + " and add new line
  model_string <- substr(model_string, 1, nchar(model_string) - 3)
  model_string <- paste0(model_string, "\n")


  # Intercepts for observed variables
  for (w in 1:waves) {
    model_string <- paste0(model_string, "x", w, "~ 1", "\n")
  }
  for (w in 1:waves) {
    model_string <- paste0(model_string, "y", w, "~ 1", "\n")
  }

  # Define latent variables
  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\np", w, " =~ 1*x", w,
      "\nq", w, " =~ 1*y", w
    )
  }

  # Stability and cross-lagged effects
  for (w in 2:waves) {
    model_string <- paste0(
      model_string,
      "\n p", w, " ~ ", stability_p, " * p", w - 1, " + ", cross_q, " * q", w - 1,
      "\n q", w, " ~ ", stability_q, " * q", w - 1, " + ", cross_p, " * p", w - 1
    )
  }

  # Confounder effects on p and q
  model_string <- paste0(
    model_string,
    "\np1 ~ ", confounder_p, " * U",
    "\nq1 ~ ", confounder_q, " * U"
  )
  for (w in 2:waves) {
    model_string <- paste0(
      model_string,
      "\n p", w, " ~ ", confounder_p, " * U",
      "\n q", w, " ~ ", confounder_q, " * U"
    )
  }


  # Variances and covariances for p and q
  for (w in 1:waves) {
    model_string <- paste0(
      model_string,
      "\n p", w, " ~~ ", variance_p, " * p", w,
      "\n q", w, " ~~ ", variance_q, " * q", w,
      "\n p", w, " ~~ ", cov_pq, " * q", w
    )
  }

  # Confounder variance
  model_string <- paste0(
    model_string,
    "\n U ~~ ", confounder_variance, " * U"
  )

  # Fix observed variable residual variances to zero
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

  # Generate data
  dat <- lavaan::simulateData(model = model_string,
                              int.ov.free = TRUE,
                              ...)

  return(list(model = model_string, data = dat))
}
