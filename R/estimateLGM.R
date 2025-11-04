#' @title estimateLGM
#' @description Generates model syntax for Linear Latent Growth Model (LGM) with options to constrain parameters and include multiple variables.
#'
#' @param waves The number of waves (time points) in the model.
#' @param variable_type Specify whether estimating "univariate" (single variable) or "bivariate" (dual variable). Default is "univariate".
#' @param time_scores A numeric vector specifying the time scores for each wave. If NULL, uses 0, 1, 2, ..., waves-1. Default is NULL.
#' @param constrain_residual_variances Logical. If TRUE, constrains residual variances to equality across waves. Default is TRUE.
#' @param constrain_residual_covariances Logical. If TRUE, constrains residual covariances to equality across waves (bivariate only). Default is TRUE.
#' @param estimate_quadratic Logical. If TRUE, includes quadratic growth factor. Default is FALSE.
#' @param center_time Logical. If TRUE, centers time scores around their mean for interpretation. Default is FALSE.
#' @param start_values Logical. If TRUE, provides starting values for key parameters. Default is FALSE.
#'
#' @return A character string containing the Lavaan model syntax for LGM.
#'
#' @details
#' This function generates the model syntax for a Linear Latent Growth Model (LGM)
#' with a specified number of waves. The LGM models systematic change over time
#' using latent intercept and slope factors.
#'
#' For univariate models:
#' - Intercept factor (I): Represents initial level/starting point
#' - Slope factor (S): Represents rate of linear change over time
#' - Quadratic factor (Q): Represents acceleration/deceleration (if specified)
#'
#' For bivariate models:
#' - All components from univariate model for both X and Y variables
#' - Covariances between growth factors across variables
#' - Optional constraints on residual parameters
#'
#' Key parameters estimated:
#' - mean_I_x, mean_I_y: Mean intercepts (initial levels)
#' - mean_S_x, mean_S_y: Mean slopes (rates of change)
#' - mean_Q_x, mean_Q_y: Mean quadratic terms (if included)
#' - var_I_x, var_I_y: Intercept variances (individual differences in initial level)
#' - var_S_x, var_S_y: Slope variances (individual differences in rate of change)
#' - var_Q_x, var_Q_y: Quadratic variances (if included)
#' - cov_IS_x, cov_IS_y: Intercept-slope covariances
#' - var_p, var_q: Latent variable unique variances
#' - cov_pq: Latent variable covariances (bivariate only)
#'
#' Time scoring options:
#' - Default: 0, 1, 2, ..., waves-1 (intercept = initial level)
#' - Custom: User-specified time points
#' - Centered: Time scores centered around their mean for easier interpretation
#'
#' The model assumes linear change over time, with optional quadratic component.
#' Individual differences in both initial levels and rates of change are estimated.
#'
#' @examples
#' # Basic univariate LGM with 5 waves
#' model_syntax <- estimateLGM(waves = 5)
#' cat(model_syntax)
#'
#' # Bivariate LGM with custom time scores
#' model_syntax <- estimateLGM(
#'   waves = 4,
#'   variable_type = "bivariate",
#'   time_scores = c(0, 1, 3, 5)
#' )
#'
#' # LGM with quadratic growth
#' model_syntax <- estimateLGM(
#'   waves = 6,
#'   estimate_quadratic = TRUE,
#'   center_time = TRUE
#' )
#'
#' @export

estimateLGM <- function(waves = 5,
                        variable_type = c("univariate", "bivariate"),
                        time_scores = NULL,
                        constrain_residual_variances = TRUE,
                        constrain_residual_covariances = TRUE,
                        estimate_quadratic = FALSE,
                        center_time = FALSE,
                        start_values = FALSE) {

  # Input validation
  if (!is.numeric(waves) || waves <= 0 || waves != as.integer(waves)) {
    stop("Error: Parameter 'waves' must be a positive integer.")
  }

  if (waves < 3) {
    stop("Error: LGM requires at least 3 waves for identification. Cannot estimate with fewer than 3 waves.")
  }

  if (estimate_quadratic && waves < 4) {
    stop("Error: Quadratic LGM requires at least 4 waves for identification.")
  }

  variable_type <- match.arg(variable_type)

  # Validate logical parameters
  logical_params <- list(
    constrain_residual_variances = constrain_residual_variances,
    constrain_residual_covariances = constrain_residual_covariances,
    estimate_quadratic = estimate_quadratic,
    center_time = center_time,
    start_values = start_values
  )

  for (param_name in names(logical_params)) {
    if (!is.logical(logical_params[[param_name]])) {
      stop(paste0("Error: Parameter '", param_name, "' must be logical (TRUE/FALSE)."))
    }
  }

  # Set up time scores
  if (is.null(time_scores)) {
    time_scores <- 0:(waves - 1)
  } else {
    if (!is.numeric(time_scores) || length(time_scores) != waves) {
      stop("Error: 'time_scores' must be a numeric vector of length equal to 'waves'.")
    }
  }

  # Center time scores if requested
  if (center_time) {
    time_scores <- time_scores - mean(time_scores)
  }

  # Calculate quadratic time scores if needed
  if (estimate_quadratic) {
    quad_scores <- time_scores^2
  }

  model_string <- ""

  if (variable_type == "univariate") {
    # ==================== UNIVARIATE LGM ====================

    # Latent variables with single indicators
    model_string <- paste0(model_string, "# Latent Variables\n")
    for (w in 1:waves) {
      model_string <- paste0(model_string, "p", w, " =~ 1*x", w, "\n")
    }

    # Define latent growth factors
    model_string <- paste0(model_string, "\n# Latent Growth Factors\n")
    model_string <- paste0(model_string, "I =~ 1*p1")
    for (w in 2:waves) {
      model_string <- paste0(model_string, " + 1*p", w)
    }
    model_string <- paste0(model_string, "\n")

    # Slope factor with time loadings
    model_string <- paste0(model_string, "S =~ ", time_scores[1], "*p1")
    for (w in 2:waves) {
      model_string <- paste0(model_string, " + ", time_scores[w], "*p", w)
    }
    model_string <- paste0(model_string, "\n")

    # Quadratic factor if requested
    if (estimate_quadratic) {
      model_string <- paste0(model_string, "Q =~ ", quad_scores[1], "*p1")
      for (w in 2:waves) {
        model_string <- paste0(model_string, " + ", quad_scores[w], "*p", w)
      }
      model_string <- paste0(model_string, "\n")
    }

    # Growth factor means
    start_intercept <- if (start_values) "start(0)*" else ""
    start_slope <- if (start_values) "start(0.1)*" else ""
    start_quad <- if (start_values) "start(0.01)*" else ""

    model_string <- paste0(model_string, "\n# Growth Factor Means\n")
    model_string <- paste0(model_string, "I ~ ", start_intercept, "mean_I*1\n")
    model_string <- paste0(model_string, "S ~ ", start_slope, "mean_S*1\n")
    if (estimate_quadratic) {
      model_string <- paste0(model_string, "Q ~ ", start_quad, "mean_Q*1\n")
    }

    # Growth factor variances and covariances
    model_string <- paste0(model_string, "\n# Growth Factor Variances and Covariances\n")
    model_string <- paste0(model_string, "I ~~ var_I*I\n")
    model_string <- paste0(model_string, "S ~~ var_S*S\n")
    model_string <- paste0(model_string, "I ~~ cov_IS*S\n")

    if (estimate_quadratic) {
      model_string <- paste0(model_string, "Q ~~ var_Q*Q\n")
      model_string <- paste0(model_string, "I ~~ cov_IQ*Q\n")
      model_string <- paste0(model_string, "S ~~ cov_SQ*Q\n")
    }

    # Latent variable means (fixed to 0)
    model_string <- paste0(model_string, "\n# Latent Variable Means (Fixed to 0)\n")
    for (w in 1:waves) {
      model_string <- paste0(model_string, "p", w, " ~ 0*1\n")
    }

    # Observed variable intercepts (fixed to 0)
    model_string <- paste0(model_string, "\n# Observed Variable Intercepts (Fixed to 0)\n")
    for (w in 1:waves) {
      model_string <- paste0(model_string, "x", w, " ~ 0*1\n")
    }

    # Latent variable variances
    model_string <- paste0(model_string, "\n# Latent Variable Variances\n")
    if (constrain_residual_variances) {
      for (w in 1:waves) {
        model_string <- paste0(model_string, "p", w, " ~~ var_p*p", w, "\n")
      }
    } else {
      for (w in 1:waves) {
        model_string <- paste0(model_string, "p", w, " ~~ p", w, "\n")
      }
    }

    # Fix observed variable residuals to 0 (perfect indicators)
    model_string <- paste0(model_string, "\n# Fix Observed Variable Residuals to Zero\n")
    for (w in 1:waves) {
      model_string <- paste0(model_string, "x", w, " ~~ 0*x", w, "\n")
    }

  } else if (variable_type == "bivariate") {
    # ==================== BIVARIATE LGM ====================

    # Latent variables with single indicators
    model_string <- paste0(model_string, "# Latent Variables\n")
    for (w in 1:waves) {
      model_string <- paste0(model_string, "p", w, " =~ 1*x", w, "\n")
      model_string <- paste0(model_string, "q", w, " =~ 1*y", w, "\n")
    }

    # Define latent growth factors for X
    model_string <- paste0(model_string, "\n# Latent Growth Factors for X\n")
    model_string <- paste0(model_string, "I_x =~ 1*p1")
    for (w in 2:waves) {
      model_string <- paste0(model_string, " + 1*p", w)
    }
    model_string <- paste0(model_string, "\n")

    model_string <- paste0(model_string, "S_x =~ ", time_scores[1], "*p1")
    for (w in 2:waves) {
      model_string <- paste0(model_string, " + ", time_scores[w], "*p", w)
    }
    model_string <- paste0(model_string, "\n")

    if (estimate_quadratic) {
      model_string <- paste0(model_string, "Q_x =~ ", quad_scores[1], "*p1")
      for (w in 2:waves) {
        model_string <- paste0(model_string, " + ", quad_scores[w], "*p", w)
      }
      model_string <- paste0(model_string, "\n")
    }

    # Define latent growth factors for Y
    model_string <- paste0(model_string, "\n# Latent Growth Factors for Y\n")
    model_string <- paste0(model_string, "I_y =~ 1*q1")
    for (w in 2:waves) {
      model_string <- paste0(model_string, " + 1*q", w)
    }
    model_string <- paste0(model_string, "\n")

    model_string <- paste0(model_string, "S_y =~ ", time_scores[1], "*q1")
    for (w in 2:waves) {
      model_string <- paste0(model_string, " + ", time_scores[w], "*q", w)
    }
    model_string <- paste0(model_string, "\n")

    if (estimate_quadratic) {
      model_string <- paste0(model_string, "Q_y =~ ", quad_scores[1], "*q1")
      for (w in 2:waves) {
        model_string <- paste0(model_string, " + ", quad_scores[w], "*q", w)
      }
      model_string <- paste0(model_string, "\n")
    }

    # Growth factor means
    start_intercept <- if (start_values) "start(0)*" else ""
    start_slope <- if (start_values) "start(0.1)*" else ""
    start_quad <- if (start_values) "start(0.01)*" else ""

    model_string <- paste0(model_string, "\n# Growth Factor Means\n")
    model_string <- paste0(model_string, "I_x ~ ", start_intercept, "mean_I_x*1\n")
    model_string <- paste0(model_string, "S_x ~ ", start_slope, "mean_S_x*1\n")
    model_string <- paste0(model_string, "I_y ~ ", start_intercept, "mean_I_y*1\n")
    model_string <- paste0(model_string, "S_y ~ ", start_slope, "mean_S_y*1\n")

    if (estimate_quadratic) {
      model_string <- paste0(model_string, "Q_x ~ ", start_quad, "mean_Q_x*1\n")
      model_string <- paste0(model_string, "Q_y ~ ", start_quad, "mean_Q_y*1\n")
    }

    # Growth factor variances
    model_string <- paste0(model_string, "\n# Growth Factor Variances\n")
    model_string <- paste0(model_string, "I_x ~~ var_I_x*I_x\n")
    model_string <- paste0(model_string, "S_x ~~ var_S_x*S_x\n")
    model_string <- paste0(model_string, "I_y ~~ var_I_y*I_y\n")
    model_string <- paste0(model_string, "S_y ~~ var_S_y*S_y\n")

    if (estimate_quadratic) {
      model_string <- paste0(model_string, "Q_x ~~ var_Q_x*Q_x\n")
      model_string <- paste0(model_string, "Q_y ~~ var_Q_y*Q_y\n")
    }

    # Within-variable covariances
    model_string <- paste0(model_string, "\n# Within-Variable Growth Factor Covariances\n")
    model_string <- paste0(model_string, "I_x ~~ cov_IS_x*S_x\n")
    model_string <- paste0(model_string, "I_y ~~ cov_IS_y*S_y\n")

    if (estimate_quadratic) {
      model_string <- paste0(model_string, "I_x ~~ cov_IQ_x*Q_x\n")
      model_string <- paste0(model_string, "S_x ~~ cov_SQ_x*Q_x\n")
      model_string <- paste0(model_string, "I_y ~~ cov_IQ_y*Q_y\n")
      model_string <- paste0(model_string, "S_y ~~ cov_SQ_y*Q_y\n")
    }

    # Between-variable covariances
    model_string <- paste0(model_string, "\n# Between-Variable Growth Factor Covariances\n")
    model_string <- paste0(model_string, "I_x ~~ cov_I_xy*I_y\n")
    model_string <- paste0(model_string, "S_x ~~ cov_S_xy*S_y\n")
    model_string <- paste0(model_string, "I_x ~~ cov_IS_xy*S_y\n")
    model_string <- paste0(model_string, "I_y ~~ cov_IS_yx*S_x\n")

    if (estimate_quadratic) {
      model_string <- paste0(model_string, "Q_x ~~ cov_Q_xy*Q_y\n")
      model_string <- paste0(model_string, "I_x ~~ cov_IQ_xy*Q_y\n")
      model_string <- paste0(model_string, "I_y ~~ cov_IQ_yx*Q_x\n")
      model_string <- paste0(model_string, "S_x ~~ cov_SQ_xy*Q_y\n")
      model_string <- paste0(model_string, "S_y ~~ cov_SQ_yx*Q_x\n")
    }

    # Latent variable means (fixed to 0)
    model_string <- paste0(model_string, "\n# Latent Variable Means (Fixed to 0)\n")
    for (w in 1:waves) {
      model_string <- paste0(model_string, "p", w, " ~ 0*1\n")
      model_string <- paste0(model_string, "q", w, " ~ 0*1\n")
    }

    # Observed variable intercepts (fixed to 0)
    model_string <- paste0(model_string, "\n# Observed Variable Intercepts (Fixed to 0)\n")
    for (w in 1:waves) {
      model_string <- paste0(model_string, "x", w, " ~ 0*1\n")
      model_string <- paste0(model_string, "y", w, " ~ 0*1\n")
    }

    # Latent variable variances
    model_string <- paste0(model_string, "\n# Latent Variable Variances\n")
    if (constrain_residual_variances) {
      for (w in 1:waves) {
        model_string <- paste0(model_string, "p", w, " ~~ var_p*p", w, "\n")
        model_string <- paste0(model_string, "q", w, " ~~ var_q*q", w, "\n")
      }
    } else {
      for (w in 1:waves) {
        model_string <- paste0(model_string, "p", w, " ~~ p", w, "\n")
        model_string <- paste0(model_string, "q", w, " ~~ q", w, "\n")
      }
    }

    # Latent variable covariances
    model_string <- paste0(model_string, "\n# Latent Variable Covariances\n")
    if (constrain_residual_covariances) {
      for (w in 1:waves) {
        model_string <- paste0(model_string, "p", w, " ~~ cov_pq*q", w, "\n")
      }
    } else {
      for (w in 1:waves) {
        model_string <- paste0(model_string, "p", w, " ~~ q", w, "\n")
      }
    }

    # Fix observed variable residuals to 0 (perfect indicators)
    model_string <- paste0(model_string, "\n# Fix Observed Variable Residuals to Zero\n")
    for (w in 1:waves) {
      model_string <- paste0(model_string, "x", w, " ~~ 0*x", w, "\n")
      model_string <- paste0(model_string, "y", w, " ~~ 0*y", w, "\n")
    }
  }

  return(model_string)
}
