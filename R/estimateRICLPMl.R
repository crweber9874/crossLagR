#' @title estimateRICLPMl
#' @description Generates lavaan model syntax for RI-CLPM with multiple indicators
#'   per construct, where the first indicator has a fixed factor loading of 1 and
#'   subsequent indicators have free factor loadings.
#'
#' @param time_varying_x A list of vectors, each containing variable names for
#'   the X indicators at each wave.
#' @param time_varying_y A list of vectors, each containing variable names for
#'   the Y indicators at each wave.
#' @param time_varying_z An optional list of vectors for Z indicators at each wave.
#' @param time_invariant_vars A vector of variable names for time-invariant variables.
#' @param waves The number of waves (time points) in the model.
#' @return A character string containing lavaan model syntax.
#'
#' @details
#' This function generates lavaan syntax for an RI-CLPM with multiple indicators,
#' following the unified framework of Usami, Murayama, & Hamaker (2019). See
#' \code{\link{estimateRICLPM}} for the unified parameter naming convention.
#'
#' Parameter labels:
#' \itemize{
#'   \item \code{ar_x}, \code{ar_y}, \code{ar_z}: Autoregressive (within-person carry-over).
#'   \item \code{cl_xy}, \code{cl_yx}, etc.: Cross-lagged effects.
#'   \item \code{I_x}, \code{I_y}, \code{I_z}: Trait factors (random intercepts).
#' }
#'
#' @references
#' Usami, S., Murayama, K., & Hamaker, E. L. (2019). A unified framework of
#'   longitudinal models to examine reciprocal relations. \emph{Psychological
#'   Methods}, 24(5), 637-657.
#'
#' @examples
#' \dontrun{
#' model_syntax <- estimateRICLPMl(
#'   time_varying_x = list(c("x1_1", "x1_2"), c("x2_1", "x2_2")),
#'   time_varying_y = list(c("y1_1", "y1_2"), c("y2_1", "y2_2")),
#'   waves = 2)
#' cat(model_syntax)
#' }
#'
#' @export
estimateRICLPMl <- function(time_varying_x = list(c("x1_1", "x1_2"), c("x2_1", "x2_2")),
                            time_varying_y = list(c("y1_1", "y1_2"), c("y2_1", "y2_2")),
                            time_varying_z = NULL,
                            time_invariant_vars = NULL,
                            waves) {
  # Input validation
  if (!is.list(time_varying_x) || !is.list(time_varying_y) ||
    length(time_varying_x) != waves || length(time_varying_y) != waves) {
    stop("Error: 'time_varying_x' and 'time_varying_y' must be lists with length equal to the number of waves.")
  }

  if (!is.null(time_varying_z) && (!is.list(time_varying_z) || length(time_varying_z) != waves)) {
    stop("Error: 'time_varying_z' must be a list with length equal to the number of waves.")
  }

  if (!is.null(time_invariant_vars) && !is.character(time_invariant_vars)) {
    stop("Error: 'time_invariant_vars' must be a character vector.")
  }

  model_string <- ""

  # Random intercepts (Trait Factors I)
  for (i in seq_along(time_varying_x[[1]])) {
    model_string <- paste0(model_string, "\nI_x =~ 1*", time_varying_x[[1]][i])
    for (w in 2:waves) {
      model_string <- paste0(model_string, " + 1*", time_varying_x[[w]][i])
    }
  }

  for (i in seq_along(time_varying_y[[1]])) {
    model_string <- paste0(model_string, "\nI_y =~ 1*", time_varying_y[[1]][i])
    for (w in 2:waves) {
      model_string <- paste0(model_string, " + 1*", time_varying_y[[w]][i])
    }
  }

  if (!is.null(time_varying_z)) {
    for (i in seq_along(time_varying_z[[1]])) {
      model_string <- paste0(model_string, "\nI_z =~ 1*", time_varying_z[[1]][i])
      for (w in 2:waves) {
        model_string <- paste0(model_string, " + 1*", time_varying_z[[w]][i])
      }
    }
  }

  # Within-person latent variables with measurement models
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

  # Autoregressive and cross-lagged paths
  for (w in 2:waves) {
    model_string <- paste0(
      model_string, "\n p", w, " ~ ar_x*p", w - 1, " + cl_yx*q", w - 1
    )
    model_string <- paste0(
      model_string, "\n q", w, " ~ ar_y*q", w - 1, " + cl_xy*p", w - 1
    )
    if (!is.null(time_varying_z)) {
      model_string <- paste0(
        model_string, "\n r", w, " ~ ar_z*r", w - 1,
        " + cl_xz*p", w - 1, " + cl_yz*q", w - 1
      )
      model_string <- paste0(
        model_string, "\n p", w, " ~ cl_zx*r", w - 1
      )
      model_string <- paste0(
        model_string, "\n q", w, " ~ cl_zy*r", w - 1
      )
    }
  }

  # Dynamic residual variances and covariances
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

  # Trait Factor (I) variances and covariances
  model_string <- paste0(
    model_string,
    "\nI_x ~~ I_x",
    "\nI_y ~~ I_y",
    "\nI_x ~~ I_y"
  )

  if (!is.null(time_varying_z)) {
    model_string <- paste0(
      model_string,
      "\nI_z ~~ I_z",
      "\nI_x ~~ I_z",
      "\nI_y ~~ I_z"
    )
  }

  return(model_string)
}
