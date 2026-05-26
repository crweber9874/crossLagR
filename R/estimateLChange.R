#' @title estimateLChange
#' @description Generate lavaan model syntax for univariate or bivariate Latent Change
#'   Score Models (LCSM).
#'
#' @param waves Integer. Number of waves (time points). Must be >= 3.
#' @param variable_type Character. Either \code{"univariate"} or \code{"bivariate"}.
#'   Default is \code{"univariate"}.
#' @param constrain_indicator_variances Logical. If TRUE, constrains observed
#'   indicator residual variances (unique factors) to equality across waves.
#'   Default is TRUE.
#' @param constrain_beta Logical. If TRUE, constrains proportional (autoregressive)
#'   effects to equality across waves. Default is TRUE.
#' @param constrain_omega Logical. If TRUE, constrains cross-lagged coupling
#'   parameters to equality across waves (bivariate only). Default is TRUE.
#' @param estimate_change_to_change Logical. If TRUE, includes change-to-change
#'   autoregressive and cross-lagged effects. Requires >= 4 waves. Default is FALSE.
#' @param constrain_change_to_change Logical. If TRUE, constrains change-to-change
#'   autoregressive effects to equality across waves. Default is FALSE.
#' @param constrain_change_cross_lag Logical. If TRUE, constrains change-to-change
#'   cross-lagged effects to equality across waves. Default is FALSE.
#' @param estimate_constant_change Logical. If TRUE, estimates a second-order
#'   constant change (accumulating) factor that loads on all latent change scores.
#'   Default is TRUE.
#'
#' @return A character string containing lavaan model syntax for the specified LCSM.
#'
#' @details
#' This function generates lavaan syntax for a Latent Change Score Model following
#' the unified framework of Usami, Murayama, & Hamaker (2019). The LCS model uses
#' the measurement equation (separating latent true scores from unique factors)
#' and the dynamic equation (with accumulating factors and lagged regression).
#'
#' Parameter labels follow the unified naming convention:
#' \itemize{
#'   \item \code{ar_x}, \code{ar_y}: Proportional change coefficients. In the LCS,
#'     the effective autoregressive parameter is \eqn{1 + ar}, so \code{ar_x} is the
#'     proportional self-feedback (denoted \eqn{\beta_x} in Usami et al.).
#'   \item \code{cl_yx}, \code{cl_xy}: Coupling parameters (cross-lagged effects).
#'     \code{cl_yx} = effect of Y level on X change; \code{cl_xy} = effect of X
#'     level on Y change. Referred to as \eqn{\gamma} in Usami et al.
#'   \item \code{u_var_x}, \code{u_var_y}: Unique factor (measurement error) variances.
#'   \item \code{A_x}, \code{A_y}: Accumulating factors (constant change factors).
#'     Unlike growth factors (I, S) in the LCM-SR which have only direct effects,
#'     these accumulating factors have both direct and indirect effects on scores
#'     through the lagged relations.
#' }
#'
#' Change-to-change parameters (\code{phi_x}, \code{phi_y}, \code{phi_cl_x},
#' \code{phi_cl_y}) are unique to the LCS model and have no direct analogue in
#' other models in the unified framework.
#'
#' @references
#' Usami, S., Murayama, K., & Hamaker, E. L. (2019). A unified framework of
#'   longitudinal models to examine reciprocal relations. \emph{Psychological
#'   Methods}, 24(5), 637-657.
#'
#' McArdle, J. J., & Hamagami, F. (2001). Latent difference score structural
#'   models for linear dynamic analyses with incomplete longitudinal data.
#'
#' @examples
#' \dontrun{
#' # Univariate LCSM with 5 waves
#' syntax <- estimateLChange(waves = 5)
#' cat(syntax)
#'
#' # Bivariate LCSM with coupling parameters
#' syntax <- estimateLChange(waves = 5, variable_type = "bivariate")
#'
#' # Fit with lavaan
#' library(lavaan)
#' fit <- lavaan(syntax, data = my_data, meanstructure = TRUE)
#' summary(fit)
#' }
#'
#' @export

estimateLChange <- function(waves = 10,
                            variable_type = c("univariate", "bivariate"),
                            constrain_indicator_variances = TRUE,
                            constrain_beta = TRUE,
                            constrain_omega = TRUE,
                            estimate_change_to_change = FALSE,
                            constrain_change_to_change = FALSE,
                            constrain_change_cross_lag = FALSE,
                            estimate_constant_change = TRUE) {
  # Validate inputs
  if (!is.numeric(waves) || waves <= 0 || waves != as.integer(waves)) {
    stop("Error: Parameter 'waves' must be a positive integer.")
  }

  if (waves < 3) {
    stop("Error: Latent change models require at least 3 waves. Cannot estimate with less than 3 waves.")
  }

  variable_type <- match.arg(variable_type)

  model_string <- ""
  if (variable_type == "univariate") {
    # ==================== SINGLE VARIABLE MODEL ====================

    # Measurement equation: latent true scores
    for (w in 1:waves) {
      model_string <- paste0(model_string, "    cf", w, " =~ 1*x", w, "\n")
    }

    # Free initial true score mean and variance, others fixed to 0 for identification
    model_string <- paste0(model_string, "    cf1 ~ 1\n")
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf", w, " ~ 0*1\n")
    }

    model_string <- paste0(model_string, "    cf1 ~~ cf1\n")
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf", w, " ~~ 0*cf", w, "\n")
    }

    for (w in 1:waves) {
      model_string <- paste0(model_string, "    x", w, " ~ 0*1\n")
    }

    # Unique factor variances (measurement equation)
    if (constrain_indicator_variances) {
      for (w in 1:waves) {
        model_string <- paste0(model_string, "    x", w, " ~~ u_var_x*x", w, "\n")
      }
    } else {
      for (w in 1:waves) {
        model_string <- paste0(model_string, "    x", w, " ~~ x", w, "\n")
      }
    }

    # Autoregressive structure (fixed = 1)
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf", w, " ~ 1*cf", w - 1, "\n")
    }

    # Latent change scores
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    ld", w, " =~ 1*cf", w, "\n")
    }

    # Latent change score means (constrained to 0)
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    ld", w, " ~ 0*1\n")
    }

    if (estimate_constant_change) {
      # Latent change score variances (constrained to 0) when estimating accumulating factor
      for (w in 2:waves) {
        model_string <- paste0(model_string, "    ld", w, " ~~ 0*ld", w, "\n")
      }

      # Accumulating factor A (loadings = 1)
      model_string <- paste0(model_string, "    A_x =~ 1*ld2\n")
      for (w in 3:waves) {
        model_string <- paste0(model_string, "          + 1*ld", w, "\n")
      }

      # Accumulating factor mean and variance
      model_string <- paste0(model_string, "    A_x ~ 1\n")
      model_string <- paste0(model_string, "    A_x ~~ A_x\n")

      # Accumulating factor covariance with the initial true score
      model_string <- paste0(model_string, "    A_x ~~ cf1\n")
    } else {
      # Free latent change score variances when not estimating accumulating factor
      if (constrain_indicator_variances) {
        for (w in 2:waves) {
          model_string <- paste0(model_string, "    ld", w, " ~~ d_var_x*ld", w, "\n")
        }
      } else {
        for (w in 2:waves) {
          model_string <- paste0(model_string, "    ld", w, " ~~ ld", w, "\n")
        }
      }
    }

    # Proportional effects (autoregressive in unified framework)
    if (constrain_beta) {
      for (w in 2:waves) {
        model_string <- paste0(model_string, "    ld", w, " ~ ar_x*cf", w - 1, "\n")
      }
    } else {
      for (w in 2:waves) {
        model_string <- paste0(model_string, "    ld", w, " ~ cf", w - 1, "\n")
      }
    }
  } else if (variable_type == "bivariate") {
    # ==================== BIVARIATE MODEL ====================
    # -------------------- X VARIABLE --------------------
    # Measurement equation: latent true scores for X
    for (w in 1:waves) {
      model_string <- paste0(model_string, "    cf_x", w, " =~ 1*x", w, "\n")
    }

    # Latent true score means (initial free, others = 0)
    model_string <- paste0(model_string, "    cf_x1 ~ 1\n")
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf_x", w, " ~ 0*1\n")
    }

    # Latent true score variances (initial free, others = 0)
    model_string <- paste0(model_string, "    cf_x1 ~~ start(15)*cf_x1\n")
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf_x", w, " ~~ 0*cf_x", w, "\n")
    }

    # Observed intercepts (fixed to 0)
    for (w in 1:waves) {
      model_string <- paste0(model_string, "    x", w, " ~ 0*1\n")
    }

    # Unique factor variances (measurement equation)
    if (constrain_indicator_variances) {
      for (w in 1:waves) {
        model_string <- paste0(model_string, "    x", w, " ~~ u_var_x*x", w, "\n")
      }
    } else {
      for (w in 1:waves) {
        model_string <- paste0(model_string, "    x", w, " ~~ x", w, "\n")
      }
    }

    # Autoregressive structure (fixed = 1)
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf_x", w, " ~ 1*cf_x", w - 1, "\n")
    }

    # Latent change scores
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    ld_x", w, " =~ 1*cf_x", w, "\n")
    }

    # Latent change score means (constrained to 0)
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    ld_x", w, " ~ 0*1\n")
    }

    if (estimate_constant_change) {
      # Latent change score variances (constrained to 0)
      for (w in 2:waves) {
        model_string <- paste0(model_string, "    ld_x", w, " ~~ 0*ld_x", w, "\n")
      }

      # Accumulating factor A_x (loadings = 1)
      model_string <- paste0(model_string, "    A_x =~ 1*ld_x2\n")
      for (w in 3:waves) {
        model_string <- paste0(model_string, "    + 1*ld_x", w, "\n")
      }

      # Accumulating factor mean and variance
      model_string <- paste0(model_string, "    A_x ~ 1\n")
      model_string <- paste0(model_string, "    A_x ~~ A_var_x*A_x\n")

      # Accumulating factor covariance with the initial true score
      model_string <- paste0(model_string, "    A_x ~~ cf_x1\n")
    } else {
      # Free latent change score variances when not estimating accumulating factor
      if (constrain_indicator_variances) {
        for (w in 2:waves) {
          model_string <- paste0(model_string, "    ld_x", w, " ~~ d_var_x*ld_x", w, "\n")
        }
      } else {
        for (w in 2:waves) {
          model_string <- paste0(model_string, "    ld_x", w, " ~~ ld_x", w, "\n")
        }
      }
    }

    # Proportional effects for X (ar_x in unified framework)
    if (constrain_beta) {
      for (w in 2:waves) {
        model_string <- paste0(model_string, "    ld_x", w, " ~ start(-0.15)*ar_x*cf_x", w - 1, "\n")
      }
    } else {
      for (w in 2:waves) {
        model_string <- paste0(model_string, "    ld_x", w, " ~ start(-0.15)*cf_x", w - 1, "\n")
      }
    }

    # -------------------- Y VARIABLE --------------------
    # Measurement equation: latent true scores for Y
    for (w in 1:waves) {
      model_string <- paste0(model_string, "    cf_y", w, " =~ 1*y", w, "\n")
    }

    # Latent true score means (initial free, others = 0)
    model_string <- paste0(model_string, "    cf_y1 ~ 1\n")
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf_y", w, " ~ 0*1\n")
    }

    # Latent true score variances (initial free, others = 0)
    model_string <- paste0(model_string, "    cf_y1 ~~ start(15)*cf_y1\n")
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf_y", w, " ~~ 0*cf_y", w, "\n")
    }

    # Observed intercepts (fixed to 0)
    for (w in 1:waves) {
      model_string <- paste0(model_string, "    y", w, " ~ 0*1\n")
    }

    # Unique factor variances (measurement equation)
    if (constrain_indicator_variances) {
      for (w in 1:waves) {
        model_string <- paste0(model_string, "    y", w, " ~~ u_var_y*y", w, "\n")
      }
    } else {
      for (w in 1:waves) {
        model_string <- paste0(model_string, "    y", w, " ~~ y", w, "\n")
      }
    }

    # Autoregressive structure (fixed = 1)
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf_y", w, " ~ 1*cf_y", w - 1, "\n")
    }

    # Latent change scores
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    ld_y", w, " =~ 1*cf_y", w, "\n")
    }

    # Latent change score means (constrained to 0)
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    ld_y", w, " ~ 0*1\n")
    }

    if (estimate_constant_change) {
      # Latent change score variances (constrained to 0)
      for (w in 2:waves) {
        model_string <- paste0(model_string, "    ld_y", w, " ~~ 0*ld_y", w, "\n")
      }

      # Accumulating factor A_y (loadings = 1)
      model_string <- paste0(model_string, "    A_y =~ 1*ld_y2\n")
      for (w in 3:waves) {
        model_string <- paste0(model_string, "          + 1*ld_y", w, "\n")
      }

      # Accumulating factor mean and variance
      model_string <- paste0(model_string, "    A_y ~ 1\n")
      model_string <- paste0(model_string, "    A_y ~~ A_var_y*A_y\n")

      # Accumulating factor covariance with the initial true score
      model_string <- paste0(model_string, "    A_y ~~ cf_y1\n")
    } else {
      # Free latent change score variances when not estimating accumulating factor
      if (constrain_indicator_variances) {
        for (w in 2:waves) {
          model_string <- paste0(model_string, "    ld_y", w, " ~~ d_var_y*ld_y", w, "\n")
        }
      } else {
        for (w in 2:waves) {
          model_string <- paste0(model_string, "    ld_y", w, " ~~ ld_y", w, "\n")
        }
      }
    }

    # Proportional effects for Y (ar_y in unified framework)
    if (constrain_beta) {
      for (w in 2:waves) {
        model_string <- paste0(model_string, "    ld_y", w, " ~ ar_y*cf_y", w - 1, "\n")
      }
    } else {
      for (w in 2:waves) {
        model_string <- paste0(model_string, "    ld_y", w, " ~ cf_y", w - 1, "\n")
      }
    }

    # -------------------- COUPLING PARAMETERS (Cross-lagged) --------------------
    # cl_yx = effect of Y level on X change; cl_xy = effect of X level on Y change
    if (constrain_omega) {
      for (w in 2:waves) {
        model_string <- paste0(model_string, "    ld_x", w, " ~ start(-0.2)*cl_yx*cf_y", w - 1, "\n")
        model_string <- paste0(model_string, "    ld_y", w, " ~ start(-0.2)*cl_xy*cf_x", w - 1, "\n")
      }
    } else {
      for (w in 2:waves) {
        model_string <- paste0(model_string, "    ld_x", w, " ~ start(-0.2)*cf_y", w - 1, "\n")
        model_string <- paste0(model_string, "    ld_y", w, " ~ start(-0.2)*cf_x", w - 1, "\n")
      }
    }

    # -------------------- CHANGE-TO-CHANGE EFFECTS --------------------
    # These are unique to the LCS model (no direct analogue in unified framework)
    if (estimate_change_to_change) {
      if (waves < 4) {
        stop("Error: Change-to-change effects require at least 4 waves.")
      }

      # Change-to-change autoregressive effects (phi)
      if (constrain_change_to_change) {
        for (w in 3:waves) {
          model_string <- paste0(model_string, "    ld_y", w, " ~ phi_y*ld_y", w - 1, "\n")
          model_string <- paste0(model_string, "    ld_x", w, " ~ phi_x*ld_x", w - 1, "\n")
        }
      } else {
        for (w in 3:waves) {
          model_string <- paste0(model_string, "    ld_y", w, " ~ ld_y", w - 1, "\n")
          model_string <- paste0(model_string, "    ld_x", w, " ~ ld_x", w - 1, "\n")
        }
      }

      # Change-to-change cross-lagged effects
      if (constrain_change_cross_lag) {
        for (w in 3:waves) {
          model_string <- paste0(model_string, "    ld_x", w, " ~ phi_cl_x*ld_y", w - 1, "\n")
          model_string <- paste0(model_string, "    ld_y", w, " ~ phi_cl_y*ld_x", w - 1, "\n")
        }
      } else {
        for (w in 3:waves) {
          model_string <- paste0(model_string, "    ld_x", w, " ~ ld_y", w - 1, "\n")
          model_string <- paste0(model_string, "    ld_y", w, " ~ ld_x", w - 1, "\n")
        }
      }
    }

    # -------------------- COVARIANCES --------------------
    if (estimate_constant_change) {
      # Covariance between accumulating factors
      model_string <- paste0(model_string, "    A_x ~~ A_y\n")

      # Covariance between initial true scores
      model_string <- paste0(model_string, "    cf_x1 ~~ cf_y1\n")

      # Cross-covariances between initial true scores and accumulating factors
      model_string <- paste0(model_string, "    cf_x1 ~~ A_y\n")
      model_string <- paste0(model_string, "    cf_y1 ~~ A_x\n")
    } else {
      # Covariance between initial true scores only
      model_string <- paste0(model_string, "    cf_x1 ~~ cf_y1\n")

      # Covariances between latent change scores when no accumulating factors
      if (constrain_indicator_variances) {
        for (w in 2:waves) {
          model_string <- paste0(model_string, "    ld_x", w, " ~~ d_cov_xy*ld_y", w, "\n")
        }
      } else {
        for (w in 2:waves) {
          model_string <- paste0(model_string, "    ld_x", w, " ~~ ld_y", w, "\n")
        }
      }
    }
  }

  return(model_string)
}
