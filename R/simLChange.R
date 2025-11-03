#' @title simLChange
#' @description Simulate data from a Latent Change Score Model
#'
#' This function creates synthetic data for either a single-variable or dual-variable
#' latent change score model with specified structural parameters that match the
#' estimateLChange function schema.
#'
#' @param waves The number of waves (time points) in the model.
#' @param variable_type Type of latent change model: "univariate" (single variable) or "bivariate" (dual variable).
#' @param initial_mean_x Mean of the initial true score for X.
#' @param initial_mean_y Mean of the initial true score for Y (for bivariate only).
#' @param initial_var_x Variance of the initial true score for X.
#' @param initial_var_y Variance of the initial true score for Y (for bivariate only).
#' @param constant_mean_x Mean of the constant change factor for X.
#' @param constant_mean_y Mean of the constant change factor for Y (for bivariate only).
#' @param latent_variance_x Variance of the constant change factor for X.
#' @param latent_variance_y Variance of the constant change factor for Y (for bivariate only).
#' @param beta_x Proportional effect for X - effect of X level on X change.
#' @param beta_y Proportional effect for Y - effect of Y level on Y change.
#' @param omega_x Coupling parameter - effect of Y level on X change.
#' @param omega_y Coupling parameter - effect of X level on Y change.
#' @param cross_change_x Change-to-change cross-lagged effect - effect of Y change on X change.
#' @param cross_change_y Change-to-change cross-lagged effect - effect of X change on Y change.
#' @param phi_x Change autoregression - effect of prior X change on current X change.
#' @param phi_y Change autoregression - effect of prior Y change on current Y change.
#' @param indicator_variance_x Measurement error variance for X indicators.
#' @param indicator_variance_y Measurement error variance for Y indicators.
#' @param cov_initial_xy Covariance between initial X and Y true scores.
#' @param cov_constant_change_xy Covariance between constant change factors for X and Y.
#' @param cov_initial_x_constant_y Covariance between initial X and constant change Y.
#' @param cov_initial_y_constant_x Covariance between initial Y and constant change X.
#' @param estimate_change_to_change Logical. Whether to include change-to-change effects.
#' @param ... Additional arguments to pass to the `lavaan::simulateData` function.
#'
#' @return A list containing two elements:
#'    * `model`: The Lavaan model syntax used for data simulation.
#'    * `data`: The simulated data in a data frame format.
#'
#' @details
#' For the single variable latent change model (variable_type = "univariate"),
#' only X-related parameters are used.
#'
#' For the dual change model (variable_type = "bivariate"), the model includes:
#' - Latent true scores (cf_x, cf_y) at each wave
#' - Latent change scores (ld_x, ld_y) from wave 2 onwards
#' - Constant change factors (general_x, general_y)
#' - Proportional effects (levels affecting own changes)
#' - Coupling effects (levels affecting other variable's changes)
#' - Change-to-change effects (lagged changes affecting current changes) if enabled
#'
#' @examples
#' # Single variable latent change model
#' single_change_data <- simLChange(
#'   waves = 5,
#'   variable_type = "univariate",
#'   beta_x = -0.2,
#'   constant_mean_x = 0.5,
#'   sample.nobs = 1000
#' )
#'
#' # Dual change model
#' dual_change_data <- simLChange(
#'   waves = 5,
#'   variable_type = "bivariate",
#'   beta_x = -0.2,
#'   beta_y = -0.3,
#'   omega_x = -0.15,
#'   omega_y = -0.15,
#'   cross_change_x = 0.1,
#'   cross_change_y = 0.1,
#'   sample.nobs = 1000
#' )
#'
#' @export
simLChange <- function(waves = 10,
                       variable_type = c("univariate", "bivariate"),
                       initial_mean_x = 0,
                       initial_mean_y = 0,
                       initial_var_x = 1,
                       initial_var_y = 1,
                       constant_mean_x = 0.5,
                       constant_mean_y = 0.5,
                       latent_variance_x = 1,
                       latent_variance_y = 1,
                       beta_x = -0.2,
                       beta_y = -0.2,
                       omega_x = -0.1,
                       omega_y = -0.1,
                       cross_change_x = 0.1,
                       cross_change_y = 0.1,
                       phi_x = 0.1,
                       phi_y = 0.1,
                       indicator_variance_x = 0.5,
                       indicator_variance_y = 0.5,
                       cov_initial_xy = 0.3,
                       cov_constant_change_xy = 0.3,
                       cov_initial_x_constant_y = 0.2,
                       cov_initial_y_constant_x = 0.2,
                       estimate_change_to_change = FALSE,
                       ...) {

  # Validate inputs
  variable_type <- match.arg(variable_type)

  if (!is.numeric(waves) || waves <= 0 || waves != as.integer(waves)) {
    stop("Error: Parameter 'waves' must be a positive integer.")
  }

  if (waves < 3) {
    stop("Error: Latent change models require at least 3 waves.")
  }

  if (estimate_change_to_change && waves < 4) {
    stop("Error: Change-to-change effects require at least 4 waves.")
  }

  model_string <- ""

  if (variable_type == "univariate") {
    # ==================== SINGLE VARIABLE MODEL ====================

    # Latent true scores (cf = change factor)
    for (w in 1:waves) {
      model_string <- paste0(model_string, "    cf", w, " =~ 1*x", w, "\n")
    }

    # Latent true score means (initial free, others = 0)
    model_string <- paste0(model_string, "    cf1 ~ ", initial_mean_x, "*1\n")
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf", w, " ~ 0*1\n")
    }

    # Latent true score variances (initial free, others = 0)
    model_string <- paste0(model_string, "    cf1 ~~ ", initial_var_x, "*cf1\n")
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf", w, " ~~ 0*cf", w, "\n")
    }

    # Observed intercepts (fixed to 0)
    for (w in 1:waves) {
      model_string <- paste0(model_string, "    x", w, " ~ 0*1\n")
    }

    # Observed residual variances
    for (w in 1:waves) {
      model_string <- paste0(model_string, "    x", w, " ~~ ", indicator_variance_x, "*x", w, "\n")
    }

    # Autoregressions (fixed = 1)
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf", w, " ~ 1*cf", w - 1, "\n")
    }

    # Latent change scores (fixed = 1)
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    ld", w, " =~ 1*cf", w, "\n")
    }

    # Latent change score means (constrained to 0)
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    ld", w, " ~ 0*1\n")
    }

    # Latent change score variances (constrained to 0)
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    ld", w, " ~~ 0*ld", w, "\n")
    }

    # Constant change factor (loadings = 1)
    model_string <- paste0(model_string, "    general =~ 1*ld2\n")
    for (w in 3:waves) {
      model_string <- paste0(model_string, "          + 1*ld", w, "\n")
    }

    # Constant change factor mean and variance
    model_string <- paste0(model_string, "    general ~ ", constant_mean_x, "*1\n")
    model_string <- paste0(model_string, "    general ~~ ", latent_variance_x, "*general\n")

    # Constant change factor covariance with the initial true score
    model_string <- paste0(model_string, "    general ~~ cf1\n")

    # Proportional effects
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    ld", w, " ~ ", beta_x, "*cf", w - 1, "\n")
    }

  } else if (variable_type == "bivariate") {
    # ==================== DUAL VARIABLE MODEL ====================

    # -------------------- X VARIABLE --------------------

    # Latent true scores for X
    for (w in 1:waves) {
      model_string <- paste0(model_string, "    cf_x", w, " =~ 1*x", w, "\n")
    }

    # Latent true score means (initial free, others = 0)
    model_string <- paste0(model_string, "    cf_x1 ~ ", initial_mean_x, "*1\n")
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf_x", w, " ~ 0*1\n")
    }

    # Latent true score variances (initial free, others = 0)
    model_string <- paste0(model_string, "    cf_x1 ~~ ", initial_var_x, "*cf_x1\n")
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf_x", w, " ~~ 0*cf_x", w, "\n")
    }

    # Observed intercepts (fixed to 0)
    for (w in 1:waves) {
      model_string <- paste0(model_string, "    x", w, " ~ 0*1\n")
    }

    # Observed residual variances
    for (w in 1:waves) {
      model_string <- paste0(model_string, "    x", w, " ~~ ", indicator_variance_x, "*x", w, "\n")
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

    # Latent change score variances (constrained to 0)
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    ld_x", w, " ~~ 0*ld_x", w, "\n")
    }

    # Constant change factor (loadings = 1)
    model_string <- paste0(model_string, "    general_x =~ 1*ld_x2\n")
    for (w in 3:waves) {
      model_string <- paste0(model_string, "          + 1*ld_x", w, "\n")
    }

    # Constant change factor mean and variance
    model_string <- paste0(model_string, "    general_x ~ ", constant_mean_x, "*1\n")
    model_string <- paste0(model_string, "    general_x ~~ ", latent_variance_x, "*general_x\n")

    # Constant change factor covariance with the initial true score
    model_string <- paste0(model_string, "    general_x ~~ cf_x1\n")

    # Proportional effects for X
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    ld_x", w, " ~ ", beta_x, "*cf_x", w - 1, "\n")
    }

    # -------------------- Y VARIABLE --------------------

    # Latent true scores for Y
    for (w in 1:waves) {
      model_string <- paste0(model_string, "    cf_y", w, " =~ 1*y", w, "\n")
    }

    # Latent true score means (initial free, others = 0)
    model_string <- paste0(model_string, "    cf_y1 ~ ", initial_mean_y, "*1\n")
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf_y", w, " ~ 0*1\n")
    }

    # Latent true score variances (initial free, others = 0)
    model_string <- paste0(model_string, "    cf_y1 ~~ ", initial_var_y, "*cf_y1\n")
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf_y", w, " ~~ 0*cf_y", w, "\n")
    }

    # Observed intercepts (fixed to 0)
    for (w in 1:waves) {
      model_string <- paste0(model_string, "    y", w, " ~ 0*1\n")
    }

    # Observed residual variances
    for (w in 1:waves) {
      model_string <- paste0(model_string, "    y", w, " ~~ ", indicator_variance_y, "*y", w, "\n")
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

    # Latent change score variances (constrained to 0)
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    ld_y", w, " ~~ 0*ld_y", w, "\n")
    }

    # Constant change factor (loadings = 1)
    model_string <- paste0(model_string, "    general_y =~ 1*ld_y2\n")
    for (w in 3:waves) {
      model_string <- paste0(model_string, "          + 1*ld_y", w, "\n")
    }

    # Constant change factor mean and variance
    model_string <- paste0(model_string, "    general_y ~ ", constant_mean_y, "*1\n")
    model_string <- paste0(model_string, "    general_y ~~ ", latent_variance_y, "*general_y\n")

    # Constant change factor covariance with the initial true score
    model_string <- paste0(model_string, "    general_y ~~ cf_y1\n")

    # Proportional effects for Y
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    ld_y", w, " ~ ", beta_y, "*cf_y", w - 1, "\n")
    }

    # -------------------- COUPLING PARAMETERS --------------------
    # Cross-lagged effects on latent change
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    ld_x", w, " ~ ", omega_y, "*cf_y", w - 1, "\n")
      model_string <- paste0(model_string, "    ld_y", w, " ~ ", omega_x, "*cf_x", w - 1, "\n")
    }

    # -------------------- CHANGE-TO-CHANGE PARAMETERS --------------------
    if (estimate_change_to_change) {
      # Change autoregression
      for (w in 3:waves) {
        model_string <- paste0(model_string, "    ld_x", w, " ~ ", phi_x, "*ld_x", w - 1, "\n")
        model_string <- paste0(model_string, "    ld_y", w, " ~ ", phi_y, "*ld_y", w - 1, "\n")
      }

      # Change-to-change cross-lagged effects
      for (w in 3:waves) {
        model_string <- paste0(model_string, "    ld_x", w, " ~ ", cross_change_x, "*ld_y", w - 1, "\n")
        model_string <- paste0(model_string, "    ld_y", w, " ~ ", cross_change_y, "*ld_x", w - 1, "\n")
      }
    }

    # -------------------- COVARIANCES --------------------
    # Covariance between constant change factors
    model_string <- paste0(model_string, "    general_x ~~ ", cov_constant_change_xy, "*general_y\n")

    # Covariance between initial true scores
    model_string <- paste0(model_string, "    cf_x1 ~~ ", cov_initial_xy, "*cf_y1\n")

    # Covariance between wave 1 scores and constant change scores
    model_string <- paste0(model_string, "    cf_x1 ~~ ", cov_initial_x_constant_y, "*general_y\n")
    model_string <- paste0(model_string, "    cf_y1 ~~ ", cov_initial_y_constant_x, "*general_x\n")
  }

  # Generate data using lavaan
  dat <- lavaan::simulateData(model = model_string,
                              int.ov.free = TRUE,
                              ...)

  return(list(model = model_string, data = dat))
}
