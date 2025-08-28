#' @title estimateStTr
#' @description Generates model syntax for state trait model
#'
#' @param data A data frame containing
#' @param time_varying_x A vector of variable names for the X variables  in data frame. Ordered in time from earliest to latest.
#' @param time_varying_y A vector of variable names for the Y variables in data frame. Ordered in time from earliest to latest.
#' @param y_vars A vector of variable names for the Y variables.
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
#' #  model_syntax = estimateRICLPM(data = my_data, x_vars = c("x1", "x2", "x3", "x4", "x5"), y_vars = c("y1", "y2", "y3", "y4", "y5"), waves = 5)
#' #  cat(model_syntax)

#' @export
#'

estimateStTr <- function(data = data,
                           waves = 10,
                           ...) {

# Trait Effects
  model_string <- ""
  model_string <- paste0(model_string, "\n trait_p =~   1*p1")

  for (w in 2:waves) {
    model_string <- paste0(
      model_string <- paste0(model_string, "+ 1*p", w)
    )
  }

  model_string <- paste0(model_string, "\n trait_q =~   1*q1")
  for (w in 2:waves) {
    model_string <- paste0(
      model_string <- paste0(model_string, "+ 1*q", w)
    )
  }
# state effects
  # p and q states, like RICLPM
  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\np", w, " =~ 1*x", w,
      "\nq", w, " =~ 1*y", w
    )
  }
  # Specify Fixed Variances at 0
  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\np", w, " ~~ 0*p", w,
      "\nq", w, " ~~ 0*q", w
    )
  }
# Specify Occasions
  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\no_x", w, " =~ 1*p", w,
      "\no_y", w, " =~ 1*q", w
    )
  }
  return(model_string)
}
#
# # dat sim
# dat <- simRICLPM(waves = 5)$data
#
# # dat["z1"] <- rnorm(nrow(dat))
#
# library(lavaan)
# lavaan(
#   model = estimateStTr(
#   waves = 5
# ),
# data = dat)
#   data = dat, fixed.x = FALSE
# ) %>% summary()
