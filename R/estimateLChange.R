#' @title estimateLChange
#' @description Syntax to estimate single and two variable latent change models. The function returns the lavaan model syntax, with options to constrain parameters, estimate various change components:
#' like the bivariate change model, the dual change model, and change to change parameters. Three or more waves of data are required for estimation absent change-to-change parameters. Change-to-Change
#' specification requires four or more waves of data.
#'
#' @param waves The number of waves (time points) in the model.
#' @param variable_type Specify whether estimating "univariate" (single variable) or "bivariate" (dual variable). The "bivariate" option allows for cross-lagged effects between two variables.
#' @param constrain_indicator_variances Logical. If TRUE, constrains indicator latent variances (true scores) to equality across waves. Default is TRUE.
#' @param constrain_beta Logical. If TRUE, constrains proportional effects (autoregressive) to equality across waves. Default is TRUE.
#' @param constrain_omega  Logical. If TRUE, constrains cross-lagged (coupling) effects to equality. Default is FALSE for bivariate.
#' @param estimate_change_to_change Logical. Extension to the Latent Change Model. Estimate lagged difference effects -- autoregressive and cross lagged. Requires 4 or more waves of data. Default is FALSE for bivariate.
#' @param constrain_change_to_change  Logical. Extension to the Latent Change Model. Constrain change to change (AR) parameters to equality. Default is FALSE for bivariate.
#' @param estimate_constant_change Logical. If TRUE, estimates constant change factors (intercept and slope components). If FALSE, omits constant change factors (pure proportional and coupling effects only). Default is TRUE.
#'
#' @return A character string containing the Lavaan model syntax for the latent change score model.


#' @details
#' This function generates the model syntax for a Latent Change Score Model with a specified
#' number of waves.
#'
#' For single variable models (variable_type = "univariate"):
#' - Latent true scores at each wave
#' - Latent change scores from wave 2 onwards
#' - A constant change factor (if estimate_constant_change = TRUE)
#' - Proportional effects (level affecting own change)
#'
#' For dual change models (variable_type = "bivariate"):
#' - All components from single variable model for both X and Y
#' - Coupling effects (X level affecting Y change and vice versa)

#' For bivariate change models (variable_type = "bivariate", constrain_change_to_change_estimates = FALSE ):
#' - All components from single variable model for both X and Y
#' - Coupling effects (X level affecting Y change and vice versa)

#' For bivariate dual change models (variable_type = "bivariate", constrain_change_to_change_estimates = TRUE ):
#' - All components from single variable model for both X and Y
#' - Cross-lagged (coupling) effects (X level affecting Y change and vice versa)
#' - Change-to-change effects (autoregressive and cross-lagged).
#'
#' When estimate_constant_change = FALSE:
#' - No constant change factors are estimated
#' - Only proportional and coupling effects drive change
#' - Useful when constant change variances are near zero or causing estimation problems
#'
#' Key parameters estimated:
#' - beta_x, beta_y: Proportional effects (level -> own change)
#' - omega_x, omega_y: Cross-lagged effects (other level -> change)
#' - phi_x, phi_y: Change-to-change effects
#' - cross_change_x, cross_change_y: Change-to-change cross-lagged effects
#' - constant_mean_x, constant_mean_y: Means of constant change factors (if estimated)
#' - latent_variance_x, latent_variance_y: Variances of constant change factors (if estimated)
#' - indicator_variance_x, indicator_variance_y: Measurement error variances

#' @examples
#' # Single variable latent change model
#' model_syntax_single <- estimateLChange(waves = 5, variable_type = "univariate")
#' cat(model_syntax_single)
#'
#' # Dual change model
#' model_syntax_dual <- estimateLChange(waves = 5, variable_type = "bivariate")
#' cat(model_syntax_dual)
#'
#' # Bivariate model without constant change factors
#' model_syntax_no_constant <- estimateLChange(
#'   waves = 5,
#'   variable_type = "bivariate",
#'   estimate_constant_change = FALSE
#' )
#' cat(model_syntax_no_constant)
#'
#' @export

estimateLChange <- function(waves = 10,
                            variable_type = c('univariate', 'bivariate'),
                            constrain_indicator_variances = TRUE,  # Constrain measurement error variances to equality
                            constrain_beta = TRUE,
                            constrain_omega = TRUE,
                            estimate_change_to_change = FALSE,
                            change_to_change = FALSE,
                            constrain_change_cross_lag = FALSE,
                            estimate_constant_change = TRUE
) {

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

    # cf = change factor
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

    if (constrain_indicator_variances) {
      for (w in 1:waves) {
        model_string <- paste0(model_string, "    x", w, " ~~ sigma2_indicator*x", w, "\n")
      }
    } else {
      for (w in 1:waves) {
        model_string <- paste0(model_string, "    x", w, " ~~ x", w, "\n")
      }
    }

    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf", w, " ~ 1*cf", w - 1, "\n")
    }

    # Latent change loadings, fixed to 1 for identification with a single variable indicator model
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    ld", w, " =~ 1*cf", w, "\n")
    }

    # Latent change score means (constrained to 0), for ident
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    ld", w, " ~ 0*1\n")
    }

    if (estimate_constant_change) {
      # Latent change score variances (constrained to 0) when estimating higher order factor, constant change
      for (w in 2:waves) {
        model_string <- paste0(model_string, "    ld", w, " ~~ 0*ld", w, "\n")
      }

      # Constant change factor (loadings = 1)
      model_string <- paste0(model_string, "    general =~ 1*ld2\n")
      for (w in 3:waves) {
        model_string <- paste0(model_string, "          + 1*ld", w, "\n")
      }

      # Constant change factor mean and variance
      model_string <- paste0(model_string, "    general ~ 1\n")
      model_string <- paste0(model_string, "    general ~~ general\n")

      # Constant change factor covariance with the initial true score
      model_string <- paste0(model_string, "    general ~~ cf1\n")
    } else {
      # Free latent change score variances when not estimating constant change
      if (constrain_indicator_variances) {
        for (w in 2:waves) {
          model_string <- paste0(model_string, "    ld", w, " ~~ var_ld*ld", w, "\n")
        }
      } else {
        for (w in 2:waves) {
          model_string <- paste0(model_string, "    ld", w, " ~~ ld", w, "\n")
        }
      }
    }

    # Proportional effects (constrained equal if specified)
    if (constrain_beta) {
      for (w in 2:waves) {
        model_string <- paste0(model_string, "    ld", w, " ~ beta*cf", w - 1, "\n")
      }
    } else {
      for (w in 2:waves) {
        model_string <- paste0(model_string, "    ld", w, " ~ cf", w - 1, "\n")
      }
    }

  } else if (variable_type == "bivariate") {

    for (w in 1:waves) {
      model_string <- paste0(model_string, "    cf_x", w, " =~ 1*x", w, "\n")
    }

    # Latent true score means (initial free, others = 0)
    model_string <- paste0(model_string, "    cf_x1 ~ indicator_mean_x*1\n")
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf_x", w, " ~ 0*1\n")
    }

    # Latent true score variances (initial free, others = 0)
    model_string <- paste0(model_string, "    cf_x1 ~~ start(15)*latent_indicator_mean_x*cf_x1\n")
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf_x", w, " ~~ 0*cf_x", w, "\n")
    }

    # Observed intercepts (fixed to 0)
    for (w in 1:waves) {
      model_string <- paste0(model_string, "    x", w, " ~ 0*1\n")
    }

    # Observed residual variances (constrained to equality)
    if (constrain_indicator_variances) {
      for (w in 1:waves) {
        model_string <- paste0(model_string, "    x", w, " ~~ indicator_variance_x*x", w, "\n")
      }
    } else {
      for (w in 1:waves) {
        model_string <- paste0(model_string, "    x", w, " ~~ x", w, "\n")
      }
    }

    # Autoregressions (fixed = 1)
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf_x", w, " ~ 1*cf_x", w - 1, "\n")
    }

    # Latent change scores (fixed = 1)
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

      # Constant change factor (loadings = 1)
      model_string <- paste0(model_string, "    general_x =~ 1*ld_x2\n")
      for (w in 3:waves) {
        model_string <- paste0(model_string, "    + 1*ld_x", w, "\n")
      }

      # Constant change factor mean
      model_string <- paste0(model_string, "    general_x ~ constant_mean_x*1\n")

      # Constant change factor variance
      model_string <- paste0(model_string, "    general_x ~~ latent_variance_x*general_x\n")

      # Constant change factor covariance with the initial true score
      model_string <- paste0(model_string, "    general_x ~~ cf_x1\n")
    } else {
      # Free latent change score variances when not estimating constant change
      if (constrain_indicator_variances) {
        for (w in 2:waves) {
          model_string <- paste0(model_string, "    ld_x", w, " ~~ var_ld_x*ld_x", w, "\n")
        }
      } else {
        for (w in 2:waves) {
          model_string <- paste0(model_string, "    ld_x", w, " ~~ ld_x", w, "\n")
        }
      }
    }

    # Proportional effects for X (constrained equal)
    if (constrain_beta) {
      for (w in 2:waves) {
        model_string <- paste0(model_string, "    ld_x", w, " ~ start(-0.15)*beta_x*cf_x", w - 1, "\n")
      }
    } else {
      for (w in 2:waves) {
        model_string <- paste0(model_string, "    ld_x", w, " ~ start(-0.15)*cf_x", w - 1, "\n")
      }
    }

    # -------------------- Y VARIABLE --------------------

    # Define the latent true scores for Y
    for (w in 1:waves) {
      model_string <- paste0(model_string, "    cf_y", w, " =~ 1*y", w, "\n")
    }

    # Latent true score means (initial free, others = 0)
    model_string <- paste0(model_string, "    cf_y1 ~ indicator_mean_y*1\n")
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf_y", w, " ~ 0*1\n")
    }

    # Latent true score variances (initial free, others = 0)
    model_string <- paste0(model_string, "    cf_y1 ~~ start(15)*latent_variance_y*cf_y1\n")
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf_y", w, " ~~ 0*cf_y", w, "\n")
    }

    # Observed intercepts (fixed to 0)
    for (w in 1:waves) {
      model_string <- paste0(model_string, "    y", w, " ~ 0*1\n")
    }

    # Observed residual variances (constrained to equality)
    if (constrain_indicator_variances) {
      for (w in 1:waves) {
        model_string <- paste0(model_string, "    y", w, " ~~ indicator_variance_y*y", w, "\n")
      }
    } else {
      for (w in 1:waves) {
        model_string <- paste0(model_string, "    y", w, " ~~ y", w, "\n")
      }
    }

    # Autoregressions (fixed = 1)
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf_y", w, " ~ 1*cf_y", w - 1, "\n")
    }

    # Latent change scores (fixed = 1)
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

      # Constant change factor (loadings = 1)
      model_string <- paste0(model_string, "    general_y =~ 1*ld_y2\n")
      for (w in 3:waves) {
        model_string <- paste0(model_string, "          + 1*ld_y", w, "\n")
      }

      # Constant change factor mean
      model_string <- paste0(model_string, "    general_y ~ constant_mean_y*1\n")

      # Constant change factor variance
      model_string <- paste0(model_string, "    general_y ~~ latent_variance_y*general_y\n")

      # Constant change factor covariance with the initial true score
      model_string <- paste0(model_string, "    general_y ~~ cf_y1\n")
    } else {
      # Free latent change score variances when not estimating constant change
      if (constrain_indicator_variances) {
        for (w in 2:waves) {
          model_string <- paste0(model_string, "    ld_y", w, " ~~ var_ld_y*ld_y", w, "\n")
        }
      } else {
        for (w in 2:waves) {
          model_string <- paste0(model_string, "    ld_y", w, " ~~ ld_y", w, "\n")
        }
      }
    }

    # Proportional effects for Y (constrained equal)
    if (constrain_beta) {
      for (w in 2:waves) {
        model_string <- paste0(model_string, "    ld_y", w, " ~ beta_y*cf_y", w - 1, "\n")
      }
    } else {
      for (w in 2:waves) {
        model_string <- paste0(model_string, "    ld_y", w, " ~ cf_y", w - 1, "\n")
      }
    }

    # -------------------- COUPLING PARAMETERS --------------------
    # This joins the two univariate LCSMs into a bivariate LCSM, with cross lagged effects on latent change
    ################################################################################################################

    if (constrain_omega) {
      for (w in 2:waves) {
        model_string <- paste0(model_string, "    ld_x", w, " ~ start(-0.2)*omega_x*cf_y", w - 1, "\n")
        model_string <- paste0(model_string, "    ld_y", w, " ~ start(-0.2)*omega_y*cf_x", w - 1, "\n")
      }
    } else {
      for (w in 2:waves) {
        model_string <- paste0(model_string, "    ld_x", w, " ~ start(-0.2)*cf_y", w - 1, "\n")
        model_string <- paste0(model_string, "    ld_y", w, " ~ start(-0.2)*cf_x", w - 1, "\n")
      }
    }

    ################################################################################################################
    ### Exensions: Coupling and Change Parameters
    ### Change to Change Additions: This is an extension to the basic bivariate LCSM
    # Lags and Autoregression of Change Parameters
    ################################################################################################################
    if(estimate_change_to_change){
      if (waves < 4) {
        stop("Error: Latent change regressions require at least 4 waves. Cannot estimate with fewer than 4 waves.")
      }
      if (constrain_change_to_change) {
        for (w in 3:waves) {
          model_string <- paste0(model_string, "    ld_y", w, " ~ phi_y*ld_y", w - 1, "\n")
        }
      } else {
        for (w in 3:waves) {
          model_string <- paste0(model_string, "    ld_y", w, " ~ ld_y", w - 1, "\n")
        }
      }

      if (constrain_change_to_change) {
        for (w in 3:waves) {
          model_string <- paste0(model_string, "    ld_x", w, " ~ phi_x*ld_x", w - 1, "\n")
        }
      } else {
        for (w in 3:waves) {
          model_string <- paste0(model_string, "    ld_x", w, " ~ ld_x", w - 1, "\n")
        }
      }
    }

    if (constrain_change_cross_lag) {
      for (w in 3:waves) {
        model_string <- paste0(model_string, "    ld_x", w, " ~ cross_change_x*ld_y", w - 1, "\n")
        model_string <- paste0(model_string, "    ld_y", w, " ~ cross-change_y*ld_x", w - 1, "\n")
      }
    } else {
      for (w in 3:waves) {
        model_string <- paste0(model_string, "    ld_x", w, " ~ ld_y", w - 1, "\n")
        model_string <- paste0(model_string, "    ld_y", w, " ~ ld_x", w - 1, "\n")
      }

    }

    # -------------------- COVARIANCES --------------------

    if (estimate_constant_change) {
      # Covariance between constant change factors
      model_string <- paste0(model_string, "    general_x ~~ general_y\n")

      # Covariance between initial true scores
      model_string <- paste0(model_string, "    cf_x1 ~~ cf_y1\n")

      # Covariance between wave 1 scores and constant change scores
      model_string <- paste0(model_string, "    cf_x1 ~~ general_y\n")
      model_string <- paste0(model_string, "    cf_y1 ~~ general_x\n")
    } else {
      # Covariance between initial true scores only
      model_string <- paste0(model_string, "    cf_x1 ~~ cf_y1\n")

      # Covariances between latent change scores when no constant change factors
      if (constrain_indicator_variances) {
        for (w in 2:waves) {
          model_string <- paste0(model_string, "    ld_x", w, " ~~ cov_ld*ld_y", w, "\n")
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
