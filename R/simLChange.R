#' @title simLChange
#' @description Simulate data from a Latent Change Score Model
#'
#' This function creates synthetic data for either a single-variable or dual-variable
#' latent change score model with specified structural parameters.
#'
#' @param waves The number of waves (time points) in the model.
#' @param model_type Type of latent change model: "latent_change" (single variable) or "dual_change" (bivariate).
#' @param initial_mean_x Mean of the initial true score for X.
#' @param initial_mean_y Mean of the initial true score for Y (for dual_change only).
#' @param initial_var_x Variance of the initial true score for X.
#' @param initial_var_y Variance of the initial true score for Y (for dual_change only).
#' @param constant_change_mean_x Mean of the constant change factor for X.
#' @param constant_change_mean_y Mean of the constant change factor for Y (for dual_change only).
#' @param constant_change_var_x Variance of the constant change factor for X.
#' @param constant_change_var_y Variance of the constant change factor for Y (for dual_change only).
#' @param ar_x Proportional effect (autoregressive) for X - effect of X level on X change.
#' @param ar_y Proportional effect (autoregressive) for Y - effect of Y level on Y change.
#' @param cl_x Coupling parameter - effect of Y level on X change.
#' @param cl_y Coupling parameter - effect of X level on Y change.
#' @param change_x Change-to-change effect - effect of Y change on X change.
#' @param change_y Change-to-change effect - effect of X change on Y change.
#' @param phi_x Change autoregression - effect of prior X change on current X change.
#' @param phi_y Change autoregression - effect of prior Y change on current Y change.
#' @param residual_variance_x Measurement error variance for X indicators.
#' @param residual_variance_y Measurement error variance for Y indicators.
#' @param cov_initial_xy Covariance between initial X and Y true scores.
#' @param cov_constant_change_xy Covariance between constant change factors for X and Y.
#' @param cov_initial_x_constant_y Covariance between initial X and constant change Y.
#' @param cov_initial_y_constant_x Covariance between initial Y and constant change X.
#' @param ... Additional arguments to pass to the `lavaan::simulateData` function.
#'
#' @return A list containing two elements:
#'    * `model`: The Lavaan model syntax used for data simulation.
#'    * `data`: The simulated data in a data frame format.
#'
#' @details
#' For the single variable latent change model (model_type = "latent_change"),
#' only X-related parameters are used.
#'
#' For the dual change model (model_type = "dual_change"), the model includes:
#' - Latent true scores (cf_x, cf_y) at each wave
#' - Latent change scores (ld_x, ld_y) from wave 2 onwards
#' - Constant change factors (general_x, general_y)
#' - Proportional effects (levels affecting own changes)
#' - Coupling effects (levels affecting other variable's changes)
#' - Change-to-change effects (lagged changes affecting current changes)
#'
#' @examples
#' # Single variable latent change model
#' single_change_data <- simLChange(
#'   waves = 5,
#'   model_type = "latent_change",
#'   ar_x = -0.2,
#'   constant_change_mean_x = 0.5,
#'   sample.nobs = 1000
#' )
#'
#' # Dual change model
#' dual_change_data <- simLChange(
#'   waves = 5,
#'   model_type = "dual_change",
#'   ar_x = -0.2,
#'   ar_y = -0.3,
#'   cl_x = -0.15,
#'   cl_y = -0.15,
#'   change_x = 0.1,
#'   change_y = 0.1,
#'   sample.nobs = 1000
#' )
#'
#' @export
simLChange <- function(waves = 10,
                       model_type = "dual_change",
                       initial_mean_x = 0,
                       initial_mean_y = 0,
                       initial_var_x = 1,
                       initial_var_y = 1,
                       constant_change_mean_x = 0.5,
                       constant_change_mean_y = 0.5,
                       constant_change_var_x = 1,
                       constant_change_var_y = 1,
                       ar_x = -0.2,
                       ar_y = -0.2,
                       cl_x = -0.1,
                       cl_y = -0.1,
                       change_x = 0.1,
                       change_y = 0.1,
                       phi_x = 0.1,
                       phi_y = 0.1,
                       residual_variance_x = 0.5,
                       residual_variance_y = 0.5,
                       cov_initial_xy = 0.3,
                       cov_constant_change_xy = 0.3,
                       cov_initial_x_constant_y = 0.2,
                       cov_initial_y_constant_x = 0.2,
                       ...) {

  # Validate inputs
  valid_model_types <- c("latent_change", "dual_change")
  if (!model_type %in% valid_model_types) {
    stop("model_type must be one of: ", paste(valid_model_types, collapse = ", "))
  }

  if (!is.numeric(waves) || waves <= 0 || waves != as.integer(waves)) {
    stop("Error: Parameter 'waves' must be a positive integer.")
  }

  if (waves < 3) {
    stop("Error: Latent change models require at least 3 waves.")
  }

  model_string <- ""

  if (model_type == "latent_change") {

        for (w in 1:waves) {
      model_string <- paste0(model_string, "    cf", w, " =~ 1*x", w, "\n")
    }

    model_string <- paste0(model_string, "    cf1 ~ ", initial_mean_x, "*1\n")
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf", w, " ~ 0*1\n")
    }

    model_string <- paste0(model_string, "    cf1 ~~ ", initial_var_x, "*cf1\n")
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf", w, " ~~ 0*cf", w, "\n")
    }

    for (w in 1:waves) {
      model_string <- paste0(model_string, "    x", w, " ~ 0*1\n")
    }

    for (w in 1:waves) {
      model_string <- paste0(model_string, "    x", w, " ~~ ", residual_variance_x, "*x", w, "\n")
    }

    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf", w, " ~ 1*cf", w - 1, "\n")
    }

    for (w in 2:waves) {
      model_string <- paste0(model_string, "    ld", w, " =~ 1*cf", w, "\n")
    }

    for (w in 2:waves) {
      model_string <- paste0(model_string, "    ld", w, " ~ 0*1\n")
    }

    for (w in 2:waves) {
      model_string <- paste0(model_string, "    ld", w, " ~~ 0*ld", w, "\n")
    }

    model_string <- paste0(model_string, "    general =~ 1*ld2\n")
    for (w in 3:waves) {
      model_string <- paste0(model_string, "          + 1*ld", w, "\n")
    }

    model_string <- paste0(model_string, "    general ~ ", constant_change_mean_x, "*1\n")

    model_string <- paste0(model_string, "    general ~~ ", constant_change_var_x, "*general\n")

    model_string <- paste0(model_string, "    general ~~ 0*cf1\n")

    for (w in 2:waves) {
      model_string <- paste0(model_string, "    ld", w, " ~ ", ar_x, "*cf", w - 1, "\n")
    }

  } else if (model_type == "dual_change") {

    for (w in 1:waves) {
      model_string <- paste0(model_string, "    cf_x", w, " =~ 1*x", w, "\n")
    }

    model_string <- paste0(model_string, "    cf_x1 ~ ", initial_mean_x, "*1\n")
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf_x", w, " ~ 0*1\n")
    }

    model_string <- paste0(model_string, "    cf_x1 ~~ ", initial_var_x, "*cf_x1\n")
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf_x", w, " ~~ 0*cf_x", w, "\n")
    }

    for (w in 1:waves) {
      model_string <- paste0(model_string, "    x", w, " ~ 0*1\n")
    }

    for (w in 1:waves) {
      model_string <- paste0(model_string, "    x", w, " ~~ ", residual_variance_x, "*x", w, "\n")
    }

    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf_x", w, " ~ 1*cf_x", w - 1, "\n")
    }

    for (w in 2:waves) {
      model_string <- paste0(model_string, "    ld_x", w, " =~ 1*cf_x", w, "\n")
    }

    for (w in 2:waves) {
      model_string <- paste0(model_string, "    ld_x", w, " ~ 0*1\n")
    }

    for (w in 2:waves) {
      model_string <- paste0(model_string, "    ld_x", w, " ~~ 0*ld_x", w, "\n")
    }

    model_string <- paste0(model_string, "    general_x =~ 1*ld_x2\n")
    for (w in 3:waves) {
      model_string <- paste0(model_string, "          + 1*ld_x", w, "\n")
    }

    model_string <- paste0(model_string, "    general_x ~ ", constant_change_mean_x, "*1\n")

    model_string <- paste0(model_string, "    general_x ~~ ", constant_change_var_x, "*general_x\n")

    model_string <- paste0(model_string, "    general_x ~~ 0*cf_x1\n")

    for (w in 2:waves) {
      model_string <- paste0(model_string, "    ld_x", w, " ~ ", ar_x, "*cf_x", w - 1, "\n")
    }

    ####################### Y ###############

    for (w in 1:waves) {
      model_string <- paste0(model_string, "    cf_y", w, " =~ 1*y", w, "\n")
    }
    # Y latent true score means (initial free, others = 0)
    model_string <- paste0(model_string, "    cf_y1 ~ ", initial_mean_y, "*1\n")
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf_y", w, " ~ 0*1\n")
    }

    model_string <- paste0(model_string, "    cf_y1 ~~ ", initial_var_y, "*cf_y1\n")
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf_y", w, " ~~ 0*cf_y", w, "\n")
    }

    for (w in 1:waves) {
      model_string <- paste0(model_string, "    y", w, " ~ 0*1\n")
    }

    for (w in 1:waves) {
      model_string <- paste0(model_string, "    y", w, " ~~ ", residual_variance_y, "*y", w, "\n")
    }

    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf_y", w, " ~ 1*cf_y", w - 1, "\n")
    }
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    ld_y", w, " =~ 1*cf_y", w, "\n")
    }

    for (w in 2:waves) {
      model_string <- paste0(model_string, "    ld_y", w, " ~ 0*1\n")
    }

    for (w in 2:waves) {
      model_string <- paste0(model_string, "    ld_y", w, " ~~ 0*ld_y", w, "\n")
    }

    model_string <- paste0(model_string, "    general_y =~ 1*ld_y2\n")
    for (w in 3:waves) {
      model_string <- paste0(model_string, "          + 1*ld_y", w, "\n")
    }

    model_string <- paste0(model_string, "    general_y ~ ", constant_change_mean_y, "*1\n")

    model_string <- paste0(model_string, "    general_y ~~ ", constant_change_var_y, "*general_y\n")

    model_string <- paste0(model_string, "    general_y ~~ 0*cf_y1\n")

    for (w in 2:waves) {
      model_string <- paste0(model_string, "    ld_y", w, " ~ ", ar_y, "*cf_y", w - 1, "\n")
    }

    for (w in 3:waves) {
      model_string <- paste0(model_string, "    ld_y", w, " ~ ", phi_y, "*ld_y", w - 1, "\n")
    }

    for (w in 3:waves) {
      model_string <- paste0(model_string, "    ld_x", w, " ~ ", phi_x, "*ld_x", w - 1, "\n")
    }

    for (w in 3:waves) {
      model_string <- paste0(model_string, "    ld_x", w, " ~ ", change_x, "*ld_y", w - 1, "\n")
      model_string <- paste0(model_string, "    ld_y", w, " ~ ", change_y, "*ld_x", w - 1, "\n")
    }

    for (w in 2:waves) {
      model_string <- paste0(model_string, "    ld_x", w, " ~ ", cl_y, "*cf_y", w - 1, "\n")
      model_string <- paste0(model_string, "    ld_y", w, " ~ ", cl_x, "*cf_x", w - 1, "\n")
    }

    model_string <- paste0(model_string, "    general_x ~~ ", cov_constant_change_xy, "*general_y\n")
    model_string <- paste0(model_string, "    cf_x1 ~~ ", cov_initial_xy, "*cf_y1\n")

    model_string <- paste0(model_string, "    cf_x1 ~~ ", cov_initial_x_constant_y, "*general_y\n")
    model_string <- paste0(model_string, "    cf_y1 ~~ ", cov_initial_y_constant_x, "*general_x\n")
  }

  dat <- lavaan::simulateData(model = model_string,
                              int.ov.free = TRUE,
                              ...)

  return(list(model = model_string, data = dat))
}
