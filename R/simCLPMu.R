#' @title simCLPMu
#' @description Simulate data from a cross-lagged panel model (CLPM) with confounder
#' In particular, this creates a simulated dataset of a cross-lagged model with a specified number of waves, structural parameters, and an optional common confounder.
#'
#' @param waves The number of waves (time points) in the model.
#' @param stability.p The stability parameter for the x variable (autoregressive effect).
#' @param stability.q The stability parameter for the y variable (autoregressive effect).
#' @param cov.pq The covariance between x and y within the same time point.
#' @param cross.q The cross-lagged effect of x on y at the next time point.
#' @param cross.p The cross-lagged effect of y on x at the next time point.
#' @param variance.p The variance of the p latent variable.
#' @param variance.q The variance of the q latent variable.
#' @param include.confounder Logical. Whether to include a common confounder u.
#' @param confounder.p The effect of the confounder u on x variables.
#' @param confounder.q The effect of the confounder u on y variables.
#' @param confounder.variance The variance of the confounder u.
#' @param confounder.stability The stability parameter for the confounder (autoregressive effect).
#' @param ... Additional arguments to pass to the `lavaan::simulateData` function.
#'
#' @return A list containing two elements:
#'    * `model`: The Lavaan model syntax used for data simulation.
#'    * `data`:  The simulated data in a data frame format.
#'
#' @export

simCLPMu <- function(waves = 10,
                    stability.p = 0.2,
                    stability.q = 0.2,
                    cross.p = 0.1,
                    cross.q = 0.1,
                    variance.p = 1,
                    variance.q = 1,
                    cov.pq = 0.1,
                    include.confounder = TRUE,
                    confounder.p = 0.3,
                    confounder.q = 0.3,
                    confounder.variance = 1,
                    confounder.stability = 0.4,
                    ...) {

  model_string <- ""

  # Intercepts for observed variables
  for (w in 1:waves) {
    model_string <- paste0(model_string, "x", w, "~ 1", "\n")
  }
  for (w in 1:waves) {
    model_string <- paste0(model_string, "y", w, "~ 1", "\n")
  }

  # If including confounder, add intercepts for u variables
  if (include.confounder) {
    for (w in 1:waves) {
      model_string <- paste0(model_string, "u", w, "~ 1", "\n")
    }
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
    if (include.confounder) {
      # Include confounder effects in the structural model
      model_string <- paste0(
        model_string,
        "\n p", w, " ~ ", stability.p, " * p", w - 1, " + ", cross.q, " * q", w - 1, " + ", confounder.p, " * u", w,
        "\n q", w, " ~ ", stability.q, " * q", w - 1, " + ", cross.p, " * p", w - 1, " + ", confounder.q, " * u", w
      )
    } else {
      # Original model without confounder
      model_string <- paste0(
        model_string,
        "\n p", w, " ~ ", stability.p, " * p", w - 1, " + ", cross.q, " * q", w - 1,
        "\n q", w, " ~ ", stability.q, " * q", w - 1, " + ", cross.p, " * p", w - 1
      )
    }
  }

  # Include confounder effects for first wave as well
  if (include.confounder) {
    model_string <- paste0(
      model_string,
      "\n p1 ~ ", confounder.p, " * u1",
      "\n q1 ~ ", confounder.q, " * u1"
    )
  }

  # Variances and covariances for p and q
  for (w in 1:waves) {
    model_string <- paste0(
      model_string,
      "\n p", w, " ~~ ", variance.p, " * p", w,
      "\n q", w, " ~~ ", variance.q, " * q", w,
      "\n p", w, " ~~ ", cov.pq, " * q", w
    )
  }

  # Confounder structure
  if (include.confounder) {
    # Confounder autoregressive effects
    for (w in 2:waves) {
      model_string <- paste0(
        model_string,
        "\n u", w, " ~ ", confounder.stability, " * u", w - 1
      )
    }

    # Confounder variances
    for (w in 1:waves) {
      model_string <- paste0(
        model_string,
        "\n u", w, " ~~ ", confounder.variance, " * u", w
      )
    }
  }

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
