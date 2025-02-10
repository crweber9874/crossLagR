#' @title model_syntax_clpm
#' @description  Generates model syntax for CLPM
#'
#' @param waves The number of waves (time points) in the model.
#' @param ... Additional arguments to pass to the `lavaan` function.
#' @return A character string containing the Lavaan model synthex for CLPM.
#'    * `model`: The Lavaan model syntax used for data simulation.
#'

#' @export
#'
#'

estimateCLPM = function(
    waves = 10) {

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

  return(model_string)
}

