#' @title estimateCLPM. This follows the Usami, Murayama, & Hamaker (2019) unified framework for longitudinal models.
#' @description Generates lavaan model syntax for the Cross-Lagged Panel Model (CLPM).
#'
#' @param waves The number of waves (time points) in the model.
#' @param constrain_beta Logical. If TRUE, constrains autoregressive effects to equality across waves. Default is TRUE.
#' @param constrain_omega Logical. If TRUE, constrains cross-lagged effects to equality across waves. Default is TRUE.
#' @param constrain_residual_variances Logical. If TRUE, constrains residual variances to equality across waves. Default is TRUE.
#' @param constrain_residual_covariances Logical. If TRUE, constrains residual covariances to equality across waves. Default is TRUE.
#' @param estimate_means Logical. If TRUE, estimates means for the first wave. Default is TRUE.
#' @param start_values Logical. If TRUE, provides starting values for key parameters. Default is FALSE.
#'
#' @return A character string containing the lavaan model syntax for the CLPM.
#'
#' @details
#' This function generates lavaan syntax for a Cross-Lagged Panel Model following
#' the unified framework of Usami, Murayama, & Hamaker (2019). In this framework,
#' the CLPM uses only the decomposition equation (group-mean centering) and the
#' dynamic equation (lagged regression), with no unique factors, trait factors,
#' or growth factors.
#'
#' Parameter labels follow the unified naming convention:
#' \itemize{
#'   \item \code{ar_y}, \code{ar_x}: Autoregressive effects. In the CLPM these
#'     reflect rank-order stability, conflating between- and within-person processes.
#'   \item \code{cl_xy}, \code{cl_yx}: Cross-lagged effects (Granger-causal paths).
#'     \code{cl_xy} = effect of X on Y; \code{cl_yx} = effect of Y on X.
#'   \item \code{d_var_x}, \code{d_var_y}: Dynamic residual (innovation) variances.
#'   \item \code{d_cov_xy}: Dynamic residual covariance (within-time association).
#' }
#'
#' @references
#' Usami, S., Murayama, K., & Hamaker, E. L. (2019). A unified framework of
#'   longitudinal models to examine reciprocal relations. \emph{Psychological
#'   Methods}, 24(5), 637-657.
#'
#' @examples
#' \dontrun{
#' # Generate CLPM syntax with 5 waves
#' syntax <- estimateCLPM(waves = 5)
#' cat(syntax)
#'
#' # Fit with lavaan
#' library(lavaan)
#' fit <- lavaan(syntax, data = my_data, meanstructure = TRUE)
#' summary(fit, fit.measures = TRUE)
#' }
#'
#' @export

estimateCLPM <- function(waves = 5,
                         constrain_beta = TRUE,
                         constrain_omega = TRUE,
                         constrain_residual_variances = TRUE,
                         constrain_residual_covariances = TRUE,
                         estimate_means = TRUE,
                         start_values = FALSE) {
  # Input validation
  if (!is.numeric(waves) || waves <= 0 || waves != as.integer(waves)) {
    stop("Error: Parameter 'waves' must be a positive integer.")
  }

  if (waves < 2) {
    stop("Error: CLPM requires at least 2 waves. Cannot estimate with fewer than 2 waves.")
  }

  logical_params <- list(
    constrain_beta = constrain_beta,
    constrain_omega = constrain_omega,
    constrain_residual_variances = constrain_residual_variances,
    constrain_residual_covariances = constrain_residual_covariances,
    estimate_means = estimate_means,
    start_values = start_values
  )

  for (param_name in names(logical_params)) {
    if (!is.logical(logical_params[[param_name]])) {
      stop(paste0("Error: Parameter '", param_name, "' must be logical (TRUE/FALSE)."))
    }
  }

  model_string <- ""

  for (w in 1:waves) {
    model_string <- paste0(
      model_string,
      "    p", w, " =~ 1*y", w, "\n",
      "    q", w, " =~ 1*x", w, "\n"
    )
  }

  if (estimate_means) {
    # Estimate means for first wave, fix others to 0
    model_string <- paste0(model_string, "    p1 ~ 1\n")
    model_string <- paste0(model_string, "    q1 ~ 1\n")

    for (w in 2:waves) {
      model_string <- paste0(model_string, "    p", w, " ~ 0*1\n")
      model_string <- paste0(model_string, "    q", w, " ~ 0*1\n")
    }
  } else {
    # Fix all means to 0
    for (w in 1:waves) {
      model_string <- paste0(model_string, "    p", w, " ~ 0*1\n")
      model_string <- paste0(model_string, "    q", w, " ~ 0*1\n")
    }
  }

  # Fix observed variable intercepts to 0
  for (w in 1:waves) {
    model_string <- paste0(model_string, "    y", w, " ~ 0*1\n")
    model_string <- paste0(model_string, "    x", w, "  ~0*1\n")
  }

  # Autoregressive and cross-lagged effects
  start_ar <- if (start_values) "start(0.1)*" else ""
  start_cl <- if (start_values) "start(0.1)*" else ""

  # Constraints to simplify model
  if (constrain_beta && constrain_omega) {
    # Both autoregressive and cross-lagged effects constrained
    for (w in 2:waves) {
      model_string <- paste0(
        model_string,
        "    p", w, " ~ ", start_ar, "ar_y*p", w - 1, " + ", start_cl, "cl_xy*q", w - 1, "\n",
        "    q", w, " ~ ", start_ar, "ar_x*q", w - 1, " + ", start_cl, "cl_yx*p", w - 1, "\n"
      )
    }
  } else if (constrain_beta && !constrain_omega) {
    # Only autoregressive effects constrained
    for (w in 2:waves) {
      model_string <- paste0(
        model_string,
        "    p", w, " ~ ", start_ar, "ar_y*p", w - 1, " + ", start_cl, "cl_xy", w, "*q", w - 1, "\n",
        "    q", w, " ~ ", start_ar, "ar_x*q", w - 1, " + ", start_cl, "cl_yx", w, "*p", w - 1, "\n"
      )
    }
  } else if (!constrain_beta && constrain_omega) {
    for (w in 2:waves) {
      model_string <- paste0(
        model_string,
        "    p", w, " ~ ", start_ar, "ar_y", w, "*p", w - 1, " + ", start_cl, "cl_xy*q", w - 1, "\n",
        "    q", w, " ~ ", start_ar, "ar_x", w, "*q", w - 1, " + ", start_cl, "cl_yx*p", w - 1, "\n"
      )
    }
  } else {
    # No constraints - all effects free
    for (w in 2:waves) {
      model_string <- paste0(
        model_string,
        "    p", w, " ~ ", start_ar, "ar_y", w, "*p", w - 1, " + ", start_cl, "cl_xy", w, "*q", w - 1, "\n",
        "    q", w, " ~ ", start_ar, "ar_x", w, "*q", w - 1, " + ", start_cl, "cl_yx", w, "*p", w - 1, "\n"
      )
    }
  }

  # Residual variances
  if (constrain_residual_variances) {
    # First wave variances (freely estimated)
    model_string <- paste0(model_string, "    p1 ~~ d_var_y1*p1\n")
    model_string <- paste0(model_string, "    q1 ~~ d_var_x1*q1\n")

    for (w in 2:waves) {
      model_string <- paste0(
        model_string,
        "    p", w, " ~~ d_var_y*p", w, "\n",
        "    q", w, " ~~ d_var_x*q", w, "\n"
      )
    }
  } else {
    for (w in 1:waves) {
      model_string <- paste0(
        model_string,
        "    p", w, " ~~ p", w, "\n",
        "    q", w, " ~~ q", w, "\n"
      )
    }
  }

  # Residual covariances (within-time correlations)
  if (constrain_residual_covariances) {
    model_string <- paste0(model_string, "    p1 ~~ d_cov_xy1*q1\n")

    for (w in 2:waves) {
      model_string <- paste0(model_string, "    p", w, " ~~ d_cov_xy*q", w, "\n")
    }
  } else {
    # All covariances freely estimated
    for (w in 1:waves) {
      model_string <- paste0(model_string, "    p", w, " ~~ q", w, "\n")
    }
  }

  for (w in 1:waves) {
    model_string <- paste0(
      model_string,
      "    y", w, " ~~ 0*y", w, "\n",
      "    x", w, " ~~ 0*x", w, "\n"
    )
  }

  return(model_string)
}
