#' @title estimateTSO
#' @description Generates lavaan model syntax for a bivariate Trait-State-Occasion
#'   (TSO) variance decomposition.
#'
#' @param waves Integer. Number of waves (time points). Must be >= 3.
#' @param constrain_state_variances Logical. If TRUE, constrains state (within-person)
#'   variances to equality across waves. Default is TRUE.
#' @param constrain_state_covariances Logical. If TRUE, constrains within-time
#'   state covariances to equality across waves. Default is TRUE.
#'
#' @return A character string containing lavaan model syntax.
#'
#' @details
#' This is a simple variance-decomposition model for the bivariate panel case
#' with a single indicator per wave per variable. Each observed variable is
#' decomposed into a time-invariant trait component and a time-varying state
#' (occasion) component:
#' \deqn{x_{it} = T_{xi} + s_{xit}}
#' \deqn{y_{it} = T_{yi} + s_{yit}}
#'
#' The trait factors \code{T_x} and \code{T_y} capture between-person variance,
#' while the state residuals capture within-person variance. The variance ratio
#' \code{var(T_x) / var(s_x)} estimates the between-to-within variance ratio.
#'
#' This model is the random-intercept-only special case of the RI-CLPM (no
#' lagged dynamics), which is what is needed for variance decomposition. With
#' a single indicator per occasion, the classic multi-indicator Cole-Martin-Steiger
#' TSO is not identified, so this simplified version is used.
#'
#' Parameter labels:
#' \itemize{
#'   \item \code{T_var_x}, \code{T_var_y}: Trait variances (between-person).
#'   \item \code{T_cov_xy}: Trait covariance.
#'   \item \code{s_var_x}, \code{s_var_y}: State variances (within-person).
#'   \item \code{s_cov_xy}: Within-time state covariance.
#' }
#'
#' @export
estimateTSO <- function(waves,
                        constrain_state_variances = TRUE,
                        constrain_state_covariances = TRUE) {
  if (!is.numeric(waves) || waves <= 0 || waves != as.integer(waves)) {
    stop("Error: 'waves' must be a positive integer.")
  }
  if (waves < 3) {
    stop("Error: TSO requires at least 3 waves.")
  }

  model <- ""

  # Trait factors
  model <- paste0(model, "# Trait factors (between-person)\n")
  model <- paste0(model, "T_x =~ 1*x1")
  for (w in 2:waves) model <- paste0(model, " + 1*x", w)
  model <- paste0(model, "\nT_y =~ 1*y1")
  for (w in 2:waves) model <- paste0(model, " + 1*y", w)
  model <- paste0(model, "\n")

  # Trait variances and covariance
  model <- paste0(model, "\n# Trait variances and covariance\n")
  model <- paste0(model, "T_x ~~ T_var_x*T_x\n")
  model <- paste0(model, "T_y ~~ T_var_y*T_y\n")
  model <- paste0(model, "T_x ~~ T_cov_xy*T_y\n")

  # State (residual) variances on observed
  model <- paste0(model, "\n# State variances (within-person)\n")
  if (constrain_state_variances) {
    for (w in 1:waves) {
      model <- paste0(model, "x", w, " ~~ s_var_x*x", w, "\n")
      model <- paste0(model, "y", w, " ~~ s_var_y*y", w, "\n")
    }
  } else {
    for (w in 1:waves) {
      model <- paste0(model, "x", w, " ~~ x", w, "\n")
      model <- paste0(model, "y", w, " ~~ y", w, "\n")
    }
  }

  # Within-time state covariances
  model <- paste0(model, "\n# Within-time state covariances\n")
  if (constrain_state_covariances) {
    for (w in 1:waves) {
      model <- paste0(model, "x", w, " ~~ s_cov_xy*y", w, "\n")
    }
  } else {
    for (w in 1:waves) {
      model <- paste0(model, "x", w, " ~~ y", w, "\n")
    }
  }

  return(model)
}
