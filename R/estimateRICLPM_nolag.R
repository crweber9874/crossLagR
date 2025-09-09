#' @title estimateRICLPM_nolag
#' @description Generates model syntax for Random Intercept Cross-Lagged Panel Model (RICLPM), no AR paths.
#'
#' @param time_varying_x A vector of variable names for the X variables, ordered from earliest to latest.
#' @param time_varying_y A vector of variable names for the Y variables, ordered from earliest to latest.
#' @param time_varying_z An optional vector of variable names for an additional Z variable, ordered from earliest to latest.
#' @param time_invariant_vars A vector of variable names for time-invariant variables.
#' @param waves The number of waves (time points) in the model.
#' @return A character string containing the Lavaan model syntax for RICLPM.
#'
#' @details
#' This function generates the model syntax for a Random Intercept Cross-Lagged Panel Model (RICLPM)
#' with a specified number of waves. The model includes latent variables for each wave, stability paths,
#' and covariances between the latent variables.
#'
#' @examples
#' \dontrun{
#'


#'
#' @export
estimateRICLPM_nolag <- function(time_varying_x,
                           time_varying_y,
                           time_varying_z = NULL,
                           time_invariant_vars = NULL,
                           waves) {
  # Input validation
  if (!is.character(time_varying_x) || !is.character(time_varying_y)) {
    stop("❌ Error: 'time_varying_x' and 'time_varying_y' must be character vectors of names in data (e.g., x1).")
  }

  if (length(time_varying_x) != waves || length(time_varying_y) != waves) {
    stop("❌ Error: The number of variables in 'time_varying_x' and 'time_varying_y' must be equal to the number of waves.")
  }

  if (!is.null(time_invariant_vars) && !is.character(time_invariant_vars)) {
    stop("❌ Error: 'time_invariant_vars' must be a character vector of names in data (e.g., x1).")
  }

  if (!is.null(time_varying_z)) {
    if (!is.character(time_varying_z)) {
      stop("❌ Error: 'time_varying_z' must be a character vector of names in data (e.g., z1).")
    }
    if (length(time_varying_z) != waves) {
      stop("❌ Error: The number of variables in 'time_varying_z' must be equal to the number of waves.")
    }
  }

  # Model string
  model_string <- ""
  model_string <- paste0(model_string, "\nBX =~  1* ", time_varying_x[1])
  for (w in 2:waves) {
    model_string <- paste0(model_string, " + 1 * ", time_varying_x[w])
  }
  model_string <- paste0(model_string, "\nBY =~   1* ", time_varying_y[1])
  for (w in 2:waves) {
    model_string <- paste0(model_string, " + 1 * ", time_varying_y[w])
  }

  if (!is.null(time_varying_z)) {
    model_string <- paste0(model_string, "\nBZ =~   1* ", time_varying_z[1])
    for (w in 2:waves) {
      model_string <- paste0(model_string, " + 1 * ", time_varying_z[w])
    }
  }

  model_string <- paste0(model_string, "\n")

  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\np", w, " =~ 1* ", time_varying_x[w],
      "\nq", w, " =~ 1* ", time_varying_y[w]
    )
    # If time varying z is present
    if (!is.null(time_varying_z)) {
      model_string <- paste0(
        model_string, "\nr", w, " =~ 1* ", time_varying_z[w]
      )
    }
  }

  # Stability and cross-lagged paths
  for (w in 2:waves) {
    model_string <- paste0(
      model_string, "\n p", w, " ~  cl_yeqn*q", w - 1
    )
    model_string <- paste0(
      model_string, "\n q", w, " ~  cl_xeqn*p", w - 1
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
