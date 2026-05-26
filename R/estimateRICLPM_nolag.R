#' @title estimateRICLPM_nolag
#' @description Generates lavaan model syntax for RI-CLPM without autoregressive
#'   paths (cross-lagged only).
#'
#' @param time_varying_x A vector of variable names for X, ordered earliest to latest.
#' @param time_varying_y A vector of variable names for Y, ordered earliest to latest.
#' @param time_varying_z An optional vector of variable names for Z.
#' @param time_invariant_vars A vector of variable names for time-invariant variables.
#' @param waves The number of waves (time points) in the model.
#' @return A character string containing lavaan model syntax.
#'
#' @details
#' This function generates lavaan syntax for an RI-CLPM with no autoregressive
#' paths (cross-lagged only), following the unified framework of Usami, Murayama,
#' & Hamaker (2019). See \code{\link{estimateRICLPM}} for the unified parameter
#' naming convention.
#'
#' Parameter labels:
#' \itemize{
#'   \item \code{cl_yx}: Cross-lagged effect of Y on X (no AR for X/Y).
#'   \item \code{cl_xy}: Cross-lagged effect of X on Y.
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
#' syntax <- estimateRICLPM_nolag(
#'   time_varying_x = paste0("x", 1:5),
#'   time_varying_y = paste0("y", 1:5),
#'   waves = 5
#' )
#' cat(syntax)
#' }
#'
#' @export
estimateRICLPM_nolag <- function(time_varying_x,
                                 time_varying_y,
                                 time_varying_z = NULL,
                                 time_invariant_vars = NULL,
                                 waves) {
  # Input validation
  if (!is.character(time_varying_x) || !is.character(time_varying_y)) {
    stop("Error: 'time_varying_x' and 'time_varying_y' must be character vectors.")
  }

  if (length(time_varying_x) != waves || length(time_varying_y) != waves) {
    stop("Error: The number of variables in 'time_varying_x' and 'time_varying_y' must equal the number of waves.")
  }

  if (!is.null(time_invariant_vars) && !is.character(time_invariant_vars)) {
    stop("Error: 'time_invariant_vars' must be a character vector.")
  }

  if (!is.null(time_varying_z)) {
    if (!is.character(time_varying_z)) {
      stop("Error: 'time_varying_z' must be a character vector.")
    }
    if (length(time_varying_z) != waves) {
      stop("Error: The number of variables in 'time_varying_z' must equal the number of waves.")
    }
  }

  # Random Intercepts (Trait Factors I)
  model_string <- ""
  model_string <- paste0(model_string, "\nI_x =~  1* ", time_varying_x[1])
  for (w in 2:waves) {
    model_string <- paste0(model_string, " + 1 * ", time_varying_x[w])
  }
  model_string <- paste0(model_string, "\nI_y =~   1* ", time_varying_y[1])
  for (w in 2:waves) {
    model_string <- paste0(model_string, " + 1 * ", time_varying_y[w])
  }

  if (!is.null(time_varying_z)) {
    model_string <- paste0(model_string, "\nI_z =~   1* ", time_varying_z[1])
    for (w in 2:waves) {
      model_string <- paste0(model_string, " + 1 * ", time_varying_z[w])
    }
  }

  model_string <- paste0(model_string, "\n")

  # Within-Person Latent Variables
  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\np", w, " =~ 1* ", time_varying_x[w],
      "\nq", w, " =~ 1* ", time_varying_y[w]
    )
    if (!is.null(time_varying_z)) {
      model_string <- paste0(
        model_string, "\nr", w, " =~ 1* ", time_varying_z[w]
      )
    }
  }

  # Cross-lagged paths only (no AR for X and Y)
  for (w in 2:waves) {
    model_string <- paste0(
      model_string, "\n p", w, " ~  cl_yx*q", w - 1
    )
    model_string <- paste0(
      model_string, "\n q", w, " ~  cl_xy*p", w - 1
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

  # Dynamic Residual Variances and Covariances
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

  # Trait Factor (I) Variances and Covariances
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
