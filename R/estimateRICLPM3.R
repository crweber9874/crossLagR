#' @title estimateRICLPM3
#' @description Generates model syntax for a three-variable Random Intercept
#'   Cross-Lagged Panel Model (RI-CLPM) with user-specified variable names.
#'
#' @param data A data frame containing the user's data.
#' @param time_varying_x A vector of variable names for X, ordered earliest to latest.
#' @param time_varying_y A vector of variable names for Y, ordered earliest to latest.
#' @param time_varying_z A vector of variable names for Z, ordered earliest to latest.
#' @param time_invariant_vars A vector of variable names for time-invariant predictors.
#' @param waves The number of waves (time points) in the model.
#' @return A character string containing lavaan model syntax for the three-variable RI-CLPM.
#'
#' @details
#' This function generates lavaan syntax for a three-variable RI-CLPM following
#' the unified framework of Usami, Murayama, & Hamaker (2019). See
#' \code{\link{estimateRICLPM}} for the unified parameter naming convention.
#'
#' Note: This variant uses equality constraints across all cross-lagged paths
#' (\code{cl}) and all autoregressive paths (\code{ar}).
#'
#' @references
#' Usami, S., Murayama, K., & Hamaker, E. L. (2019). A unified framework of
#'   longitudinal models to examine reciprocal relations. \emph{Psychological
#'   Methods}, 24(5), 637-657.
#'
#' @examples
#' # Not Run
#' # model_syntax = estimateRICLPM3(data = my_data,
#' #   time_varying_x = paste0("x", 1:5),
#' #   time_varying_y = paste0("y", 1:5),
#' #   time_varying_z = paste0("z", 1:5), waves = 5)
#'
#' @export

estimateRICLPM3 <- function(data, time_varying_x, time_varying_y, time_varying_z, time_invariant_vars = NULL, waves = 10) {
  if (length(time_varying_x) != waves || length(time_varying_y) != waves || length(time_varying_z) != waves) {
    stop("The length of time_varying_x, time_varying_y, and time_varying_z must be equal to the number of waves.")
  }

  model_string <- ""

  # Random Intercepts (Trait Factors I)
  model_string <- paste0(model_string, "\n I_x =~  1* ", time_varying_x[1])
  for (w in 2:waves) {
    model_string <- paste0(model_string, " + 1 * ", time_varying_x[w])
  }
  model_string <- paste0(model_string, "\n I_y =~   1* ", time_varying_y[1])
  for (w in 2:waves) {
    model_string <- paste0(model_string, " + 1 * ", time_varying_y[w])
  }
  model_string <- paste0(model_string, "\n I_z =~   1* ", time_varying_z[1])
  for (w in 2:waves) {
    model_string <- paste0(model_string, " + 1 * ", time_varying_z[w])
  }
  model_string <- paste0(model_string, "\n")

  # Within-Person Latent Variables
  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\np", w, " =~ 1* ", time_varying_x[w],
      "\nq", w, " =~ 1* ", time_varying_y[w],
      "\nr", w, " =~ 1* ", time_varying_z[w]
    )
  }

  # Autoregressive and Cross-lagged Effects
  for (w in 2:waves) {
    model_string <- paste0(
      model_string, "\n p", w, " ~  ", "ar*p", w - 1, " + ", "cl*q", w - 1, " + ", "cl*r", w - 1,
      "\n q", w, " ~ ", "ar*q", w - 1, " + ", "cl*p", w - 1, " + ", "cl*r", w - 1,
      "\n r", w, " ~ ", "ar*r", w - 1, " + ", "cl*p", w - 1, " + ", "cl*q", w - 1
    )
  }

  # Dynamic Residual Variances and Covariances
  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\n p", w, " ~~ ", " p", w,
      "\n q", w, " ~~ ", "q", w,
      "\n r", w, " ~~ ", "r", w,
      "\n p", w, " ~~ ", " q", w,
      "\n p", w, " ~~ ", " r", w,
      "\n q", w, " ~~ ", " r", w
    )
  }

  # Time-invariant variables
  if (!is.null(time_invariant_vars)) {
    for (var in time_invariant_vars) {
      model_string <- paste0(model_string, "\n ", var, " ~~ ", var)
      for (w in 1:waves) {
        model_string <- paste0(
          model_string, "\n p", w, " ~ ", var,
          "\n q", w, " ~ ", var,
          "\n r", w, " ~ ", var
        )
      }
    }
  }

  # Trait Factor (I) Variances and Covariances
  model_string <- paste0(
    model_string,
    "\n I_x ~~", "I_x",
    "\n I_y ~~", "I_y",
    "\n I_z ~~", "I_z",
    "\n I_x ~~", "I_y",
    "\n I_x ~~", "I_z",
    "\n I_y ~~", "I_z"
  )

  return(model_string)
}
