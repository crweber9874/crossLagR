#' @title model_syntax_clpm
#' @description Generates model syntax for  Cross-Lagged Panel Model (CLPM)
#'
#' @param waves The number of waves (time points) in the model.
#' @return A character string containing the Lavaan model syntax for RICLPM.
#'    * `model`: The Lavaan model syntax used for data simulation.
#'
#' @details
#' This function generates the model syntax for a Cross-Lagged Panel Model (CLPM)
#' with a specified number of waves.
#'
#'
#' @examples
#' # Not Run
#' # Generate model syntax for a RICLPM with 5 waves
#' #  model_syntax = estimateCLPM(waves = 5)
#' #  cat(model_syntax)

estimateCLPM = function(waves = 10) {
  # Check that 'waves' is a positive integer
  if (!is.numeric(waves) || waves <= 0 || waves != as.integer(waves)) {
    stop("Error: Parameter 'waves' must be a positive integer.")
  }

  model_string <- ""

  model_string <- paste0(model_string, "\n")

  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\np", w, " =~ 1*x", w,
      "\nq", w, " =~ 1*y", w
    )
  }

  # Stability
  for (w in 2:waves) {
    model_string <- paste0(
      model_string, "\n p", w, " ~  ", "ar_yeqn*p", w - 1, " + ", "cl_yeqn*q", w - 1,
      "\n q", w, " ~ ", "ar_xeqn*q", w - 1, " + ", "cl_xeqn*p", w - 1
    )
  }

  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\n p", w, " ~~ ", " p", w,
      "\n q", w, " ~~ ", " q", w,
      "\n p", w, " ~~ ", " q", w
    )
  }

  return(model_string)
}
