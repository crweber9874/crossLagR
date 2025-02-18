#' @title estimateRICLPM
#' @description Generates model syntax for Random Intercept Cross-Lagged Panel Model (RICLPM)
#'
#' @param data A data frame containing the user's data.
#' @param time_varying_x A vector of variable names for the X variables in the data frame. Ordered in time from earliest to latest.
#' @param time_varying_y A vector of variable names for the Y variables in the data frame. Ordered in time from earliest to latest.
#' @param time_varying_z A vector of variable names for the Z variables in the data frame. Ordered in time from earliest to latest.
#' @param time_invariant_vars A vector of variable names for the time-invariant variables.
#' @param waves The number of waves (time points) in the model.
#' @return A character string containing the Lavaan model syntax for RI-CLPM.
#'    * `model`: The Lavaan model syntax used for data simulation.
#'
#' @details
#' This function generates the model syntax for a Random Intercept Cross-Lagged Panel Model (RICLPM)
#' with a specified number of waves. The model includes latent variables for each wave, stability paths,
#' and covariances between the latent variables.
#'
#'
#' @examples
#' # Not Run
#' # Generate model syntax for a RICLPM with 5 waves
#' #  model_syntax = estimateRICLPM(data = my_data, time_varying_x = c("x1", "x2", "x3", "x4", "x5"), time_varying_y = c("y1", "y2", "y3", "y4", "y5"), time_varying_z = c("z1", "z2", "z3", "z4", "z5"), waves = 5, time_invariant_vars = c("w1", "w2"))
#' #  cat(model_syntax)

#' @export
#'

estimateRICLPM3 <- function(data, time_varying_x, time_varying_y, time_varying_z, time_invariant_vars = NULL, waves = 10) {
  if (length(time_varying_x) != waves || length(time_varying_y) != waves || length(time_varying_z) != waves) {
    stop("The length of time_varying_x, time_varying_y, and time_varying_z must be equal to the number of waves.")
  }

  model_string <- ""
  model_string <- paste0(model_string, "\n BX =~  1* ", time_varying_x[1])
  for (w in 2:waves) {
    model_string <- paste0(model_string, " + 1 * ", time_varying_x[w])
  }
  model_string <- paste0(model_string, "\n BY =~   1* ", time_varying_y[1])
  for (w in 2:waves) {
    model_string <- paste0(model_string, " + 1 * ", time_varying_y[w])
  }
  model_string <- paste0(model_string, "\n BZ =~   1* ", time_varying_z[1])
  for (w in 2:waves) {
    model_string <- paste0(model_string, " + 1 * ", time_varying_z[w])
  }
  model_string <- paste0(model_string, "\n")

  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\np", w, " =~ 1* ", time_varying_x[w],
      "\nq", w, " =~ 1* ", time_varying_y[w],
      "\nr", w, " =~ 1* ", time_varying_z[w]
    )
  }

  # Stability
  for (w in 2:waves) {
    model_string <- paste0(
      model_string, "\n p", w, " ~  ", "ar*p", w - 1, " + ", "cl*q", w - 1, " + ", "cl*r", w - 1,
      "\n q", w, " ~ ", "ar*q", w - 1, " + ", "cl*p", w - 1, " + ", "cl*r", w - 1,
      "\n r", w, " ~ ", "ar*r", w - 1, " + ", "cl*p", w - 1, " + ", "cl*q", w - 1
    )
  }

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

  model_string <- paste0(
    model_string,
    "\n BX ~~", "BX",
    "\n BY ~~", "BY",
    "\n BZ ~~", "BZ",
    "\n BX ~~", "BY",
    "\n BX ~~", "BZ",
    "\n BY ~~", "BZ"
  )

  return(model_string)
}
lavaan:::lavaan(
  estimateRICLPM3(
    data = simRICLPM3(waves = 5)$data,
    time_varying_x = c("x1", "x2", "x3", "x4", "x5"),
    time_varying_y = c("y1", "y2", "y3", "y4", "y5"),
    time_varying_z = c("z1", "z2", "z3", "z4", "z5"),
    waves = 5
  ),
  data = simRICLPM3(waves = 5)$data
) %>% summary()
