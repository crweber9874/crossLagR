#' @title estimateRICLPMl
#' @description Generates model syntax for Random Intercept Cross-Lagged Panel Model (RICLPM) with multiple indicators per construct, where the first indicator of each latent variable has a fixed factor loading of 1, and subsequent indicators have free factor loadings.
#'
#' @param time_varying_x A list of vectors, each containing variable names for the X indicators at each wave. This is clunky, but basically declare the
#' @param time_varying_y A list of vectors, each containing variable names for the Y indicators at each wave.
#' @param time_varying_z An optional list of vectors, each containing variable names for the Z indicators at each wave.
#' @param time_invariant_vars A vector of variable names for time-invariant variables.
#' @param waves The number of waves (time points) in the model.
#' @return A character string containing the Lavaan model syntax for RICLPM with freed factor loadings after the first indicator.
#'
#' @details
#' This function generates the model syntax for a Random Intercept Cross-Lagged Panel Model (RICLPM)
#' with a specified number of waves and multiple indicators for each construct. Importantly, for each
#' latent variable at each wave, the factor loading of the first indicator is fixed to 1 for identification
#' purposes, while the factor loadings of subsequent indicators are freely estimated.
#'
#' @examples
#' \dontrun{
 # model_syntax <- estimateRICLPMl(
 #   time_varying_x = list(c("x1_1", "x2_1", "x3_1"), c("x1_2", "x2_2", "x3_2")),
 #   time_varying_y = list(c("y1_1", "y2_1"), c("y1_2", "y2_2")),
 #    time_varying_z = list(c("z1_1", "z2_1"), c("z1_2", "z2_2")),
 #   waves = 2)
#' cat(model_syntax)
#' }
#'
#' @export
estimateRICLPMl <- function(time_varying_x = list(c("x1_1", "x1_2"), c("x2_1", "x2_2")),
                            time_varying_y = list(c("y1_1", "y1_2"), c("y2_1", "y2_2")) ,
                            time_varying_z = NULL,
                            time_invariant_vars = NULL,
                            waves) {
  # Input validation
  if (!is.list(time_varying_x) || !is.list(time_varying_y) ||
      length(time_varying_x) != waves || length(time_varying_y) != waves) {
    stop("❌ Error: 'time_varying_x' and 'time_varying_y' must be lists with length equal to the number of waves.")
  }

  if (!is.null(time_varying_z) && (!is.list(time_varying_z) || length(time_varying_z) != waves)) {
    stop("❌ Error: 'time_varying_z' must be a list with length equal to the number of waves.")
  }

  if (!is.null(time_invariant_vars) && !is.character(time_invariant_vars)) {
    stop("❌ Error: 'time_invariant_vars' must be a character vector.")
  }

  # Model string generation
  model_string <- ""

  # Random intercepts
  for (i in seq_along(time_varying_x[[1]])) {
    model_string <- paste0(model_string, "\nBX =~ 1*", time_varying_x[[1]][i])
    for (w in 2:waves) {
      model_string <- paste0(model_string, " + 1*", time_varying_x[[w]][i])
    }
  }

  for (i in seq_along(time_varying_y[[1]])) {
    model_string <- paste0(model_string, "\nBY =~ 1*", time_varying_y[[1]][i])
    for (w in 2:waves) {
      model_string <- paste0(model_string, " + 1*", time_varying_y[[w]][i])
    }
  }

  if (!is.null(time_varying_z)) {
    for (i in seq_along(time_varying_z[[1]])) {
      model_string <- paste0(model_string, "\nBZ =~ 1*", time_varying_z[[1]][i])
      for (w in 2:waves) {
        model_string <- paste0(model_string, " + 1*", time_varying_z[[w]][i])
      }
    }
  }

  # Latent variables specified with measurement models
  for (w in 1:waves) {
    model_string <- paste0(model_string, "\n")
    # X latent variable
    x_indicators <- time_varying_x[[w]]
    model_string <- paste0(model_string, "\np", w, " =~ 1*", x_indicators[1])
    if (length(x_indicators) > 1) {
      for (i in 2:length(x_indicators)) {
        model_string <- paste0(model_string, " + lx", w, "_", i, "*", x_indicators[i])
      }
    }

    # Y latent variable
    y_indicators <- time_varying_y[[w]]
    model_string <- paste0(model_string, "\nq", w, " =~ 1*", y_indicators[1])
    if (length(y_indicators) > 1) {
      for (i in 2:length(y_indicators)) {
        model_string <- paste0(model_string, " + ly", w, "_", i, "*", y_indicators[i])
      }
    }

    # Z latent variable (if present)
    if (!is.null(time_varying_z)) {
      z_indicators <- time_varying_z[[w]]
      model_string <- paste0(model_string, "\nr", w, " =~ 1*", z_indicators[1])
      if (length(z_indicators) > 1) {
        for (i in 2:length(z_indicators)) {
          model_string <- paste0(model_string, " + lz", w, "_", i, "*", z_indicators[i])
        }
      }
    }
  }

  # Stability and cross-lagged paths
  for (w in 2:waves) {
    model_string <- paste0(
      model_string, "\n p", w, " ~ ar_yeqn*p", w - 1, " + cl_yeqn*q", w - 1
    )
    model_string <- paste0(
      model_string, "\n q", w, " ~ ar_xeqn*q", w - 1, " + cl_xeqn*p", w - 1
    )
    if (!is.null(time_varying_z)) {
      model_string <- paste0(
        model_string, "\n r", w, " ~ ar_zeqn*r", w - 1,
        " + cl_zx*p", w - 1, " + cl_zy*q", w - 1
      )
      model_string <- paste0(
        model_string, "\n p", w, " ~ cl_xz*r", w - 1
      )
      model_string <- paste0(
        model_string, "\n q", w, " ~ cl_yz*r", w - 1
      )
    }
  }

  # Covariances
  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\n p", w, " ~~ p", w,
      "\n q", w, " ~~ q", w,
      "\n p", w, " ~~ q", w
    )
    if (!is.null(time_varying_z)) {
      model_string <- paste0(
        model_string, "\n r", w, " ~~ r", w,
        "\n p", w, " ~~ r", w,
        "\n q", w, " ~~ r", w
      )
    }
  }

  # Time-invariant variables
  if (!is.null(time_invariant_vars)) {
    for (var in time_invariant_vars) {
      model_string <- paste0(model_string, "\n", var, " ~~ ", var)
      for (w in 1:waves) {
        model_string <- paste0(
          model_string, "\n p", w, " ~ ", var,
          "\n q", w, " ~ ", var
        )
        if (!is.null(time_varying_z)) {
          model_string <- paste0(
            model_string, "\n r", w, " ~ ", var
          )
        }
      }
    }
  }


## Add this appropriately, no correlation between the random
  RIX1 + RIX2 + RIX3 + RIY1 + RIY2 + RIY3 ~~ 0*WFY1 + 0*WFX1

  model_string <- paste0(
    model_string,
    "\nBX ~~ BX",
    "\nBY ~~ BY",
    "\nBX ~~ BY"
  )

  if (!is.null(time_varying_z)) {
    model_string <- paste0(
      model_string,
      "\nBZ ~~ BZ",
      "\nBX ~~ BZ",
      "\nBY ~~ BZ"
    )
  }

  return(model_string)
}
