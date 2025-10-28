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
#' @param means_p Vector of means for p across waves. If NULL, defaults to 0 for all waves.
#' @param means_q Vector of means for q across waves. If NULL, defaults to 0 for all waves.
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
                    means_p = NULL,
                    means_q = NULL,
                    ...) {

  # Set default means if not provided
  if (is.null(means_p)) {
    means_p <- rep(0, waves)
  }
  if (is.null(means_q)) {
    means_q <- rep(0, waves)
  }

  # Check that means vectors have correct length
  if (length(means_p) != waves) {
    stop("means_p must have length equal to waves")
  }
  if (length(means_q) != waves) {
    stop("means_q must have length equal to waves")
  }

  model_string <- ""

  # Set observed variable intercepts to 0 (since latents will have the means)
  for (w in 1:waves) {
    model_string <- paste0(model_string, "x", w, "~ 0*1", "\n")
  }
  for (w in 1:waves) {
    model_string <- paste0(model_string, "y", w, "~ 0*1", "\n")
  }

  # Define latent variables
  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\np", w, " =~ 1*y", w,
      "\nq", w, " =~ 1*x", w
    )
  }

  # Set means for latent variables
  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\np", w, " ~ ", means_p[w], "*1",
      "\nq", w, " ~ ", means_q[w], "*1", "\n"
    )
  }

  # Stability and cross-lagged effects
  for (w in 2:waves) {
    model_string <- paste0(
      model_string, "\n p", w, " ~ ", stability_p, " * p", w - 1, " + ", cross_q, " * q", w - 1,
      "\n q", w, " ~  ", stability_q, " * q", w - 1, " + ", cross_p, " * p", w - 1
    )
  }

  # Variances and covariances
  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\n p", w, " ~~ ", variance_p, " * p", w,
      "\n q", w, " ~~ ", variance_q, " * q", w,
      "\n p", w, " ~~ ", cov_pq, " * q", w
    )
  }

  # Set observed variable error variances to 0
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

# Example usage:
# simCLPM(waves = 5,
#         means_p = c(0, 0.5, 1, 1.5, 2),  # increasing means for p
#         means_q = c(2, 1.5, 1, 0.5, 0))  # decreasing means for q
