#' @title simCLPMu
#' @description Simulate data from a cross-lagged panel model (CLPM) with an
#'   unmeasured confounder. `confounder_type` selects whether the confounder is
#'   time-variant (an AR(1) process `u1..uT`, one draw per wave) or
#'   time-invariant (a single latent `U` common to all waves).
#'
#' @param waves The number of waves (time points) in the model.
#' @param stability_p The stability parameter for the x variable (autoregressive effect).
#' @param stability_q The stability parameter for the y variable (autoregressive effect).
#' @param cov_pq The covariance between x and y within the same time point.
#' @param cross_q The cross-lagged effect of x on y at the next time point.
#' @param cross_p The cross-lagged effect of y on x at the next time point.
#' @param variance_p The variance of the p latent variable.
#' @param variance_q The variance of the q latent variable.
#' @param confounder_p The effect of the confounder on x variables.
#' @param confounder_q The effect of the confounder on y variables.
#' @param confounder_variance The variance of the confounder.
#' @param confounder_stability Autoregressive effect of the confounder. Used
#'   only when `confounder_type = "time_variant"`.
#' @param confounder_type Either `"time_variant"` (default) or `"time_invariant"`.
#' @param sample.nobs Number of observations to simulate.
#' @param ... Additional arguments passed to `lavaan::simulateData`.
#'
#' @return A list containing two elements:
#'    * `model`: The lavaan model syntax used for data simulation.
#'    * `data`:  The simulated data in a data frame format.
#'
#' @export
simCLPMu <- function(waves = 10,
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
                     confounder_stability = 0.4,
                     confounder_type = c("time_variant", "time_invariant"),
                     sample.nobs = 500,
                     ...) {

  confounder_type <- match.arg(confounder_type)
  time_variant <- confounder_type == "time_variant"

  model_string <- ""

  # Time-invariant confounder is a single latent U common to every wave.
  if (!time_variant) {
    loadings <- paste0("1*x", 1:waves, " + 1*y", 1:waves, collapse = " + ")
    model_string <- paste0(model_string, "\nU =~ ", loadings, "\n")
  }

  # Intercepts for observed variables
  for (w in 1:waves) {
    model_string <- paste0(model_string, "x", w, "~ 1", "\n")
  }
  for (w in 1:waves) {
    model_string <- paste0(model_string, "y", w, "~ 1", "\n")
  }

  if (time_variant) {
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

  # Stability, cross-lagged, and confounder effects
  u_at <- function(w) if (time_variant) paste0("u", w) else "U"
  for (w in 2:waves) {
    model_string <- paste0(
      model_string,
      "\n p", w, " ~ ", stability_p, " * p", w - 1, " + ", cross_q, " * q", w - 1,
      " + ", confounder_p, " * ", u_at(w),
      "\n q", w, " ~ ", stability_q, " * q", w - 1, " + ", cross_p, " * p", w - 1,
      " + ", confounder_q, " * ", u_at(w)
    )
  }
  model_string <- paste0(
    model_string,
    "\n p1 ~ ", confounder_p, " * ", u_at(1),
    "\n q1 ~ ", confounder_q, " * ", u_at(1)
  )

  # Variances and covariances for p and q
  for (w in 1:waves) {
    model_string <- paste0(
      model_string,
      "\n p", w, " ~~ ", variance_p, " * p", w,
      "\n q", w, " ~~ ", variance_q, " * q", w,
      "\n p", w, " ~~ ", cov_pq, " * q", w
    )
  }

  # Confounder structure
  if (time_variant) {
    for (w in 2:waves) {
      model_string <- paste0(
        model_string, "\n u", w, " ~ ", confounder_stability, " * u", w - 1
      )
    }
    for (w in 1:waves) {
      model_string <- paste0(
        model_string, "\n u", w, " ~~ ", confounder_variance, " * u", w
      )
    }
  } else {
    model_string <- paste0(
      model_string, "\n U ~~ ", confounder_variance, " * U"
    )
  }

  # Fix observed variable residual variances to zero
  for (w in 1:waves) {
    model_string <- paste0(model_string, "\nx", w, "~~", "0*", "x", w)
  }
  for (w in 1:waves) {
    model_string <- paste0(model_string, "\ny", w, "~~", "0*", "y", w)
  }

  dat <- lavaan::simulateData(model = model_string,
                              sample.nobs = sample.nobs,
                              int.ov.free = TRUE,
                              ...)

  return(list(model = model_string, data = dat))
}

#' @title simCLPM_timeInvariantU
#' @description Deprecated. Use
#'   `simCLPMu(confounder_type = "time_invariant")` instead.
#' @inheritParams simCLPMu
#' @return See [simCLPMu()].
#' @keywords internal
#' @export
simCLPM_timeInvariantU <- function(waves = 10, ...) {
  .Deprecated("simCLPMu(confounder_type = \"time_invariant\")")
  simCLPMu(waves = waves, confounder_type = "time_invariant", ...)
}
