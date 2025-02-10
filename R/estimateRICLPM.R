#' @title estimateRICLPM
#' @description Generates model syntax for Random Intercept Cross-Lagged Panel Model (RICLPM)
#'
#' @param waves The number of waves (time points) in the model.
#' @return A character string containing the Lavaan model syntax for RICLPM.
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
#' #  model_syntax = estimateRICLPM(waves = 5)
#' #  cat(model_syntax)

#' @export
#'

estimateRICLPM = function(
    waves = 10) {

  model_string <- ""
  model_string <- paste0(model_string, "\n BX =~  1* x1")
  for (w in 2:waves) {
    model_string <- paste0(model_string, " + 1 *x", w, "")
  }
  model_string <- paste0(model_string, "\n BY =~   1* y1")
  for (w in 2:waves) {
    model_string <- paste0(model_string, " + 1 * y", w)
  }
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
      model_string, "\n p", w, " ~  ", "ar*p", w - 1, " + ","cl*q", w - 1,
      "\n q", w, " ~ ", "ar*q", w - 1, " + ", "cl*p", w - 1
    )
  }

  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\n p", w, " ~~ ",  " p", w,
      "\n q", w, " ~~ ",  "q", w,
      "\n p", w, " ~~ ",  " q", w
    )
  }

  model_string <- paste0(
    model_string,
    "\n BX ~~", "BX",
    "\n BY ~~", "BY",
    "\n BX ~~", "BY"
  )

  return(model_string)
}

