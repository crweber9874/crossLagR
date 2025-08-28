#' @title simulateRICLPMl
#' @description Simulate data from a multiple indicator Random Intercept Cross-Lagged Regression (RI-CLPM) model,
#' where the first indicator of each latent variable has a fixed factor loading of 1, and subsequent
#' indicators have factor loadings based on the provided parameters. Includes an optional third
#' variable that can interact with X and Y, also following a cross lagged structure.
#'
#' @param waves The number of waves (time points) in the model.
#' @param num_indicators_x The number of indicators for the X variable.
#' @param num_indicators_y The number of indicators for the Y variable.
#' @param num_indicators_z The number of indicators for the Z variable (if specified by user).
#' @param stability.p The stability parameter for the latent X variable. This is the autoregressive parameter for X. This should be specified with care, particularly when comparing a CLPM to RI-CLPM.
#' @param stability.q The stability parameter for the latent Y variable. This is the autoregressive parameter. This should be specified with care, particularly when comparing a CLPM to RI-CLPM.
#' @param stability.r The stability parameter for the latent Z variable (autoregressive effect, if applicable).
#' @param cross.p The cross-lagged effect of latent Y on latent X at the next time point.
#' @param cross.q The cross-lagged effect of latent X on latent Y at the next time point.
#' @param cross.zx The cross-lagged effect of latent Z on latent X at the next time point (if applicable).
#' @param cross.zy The cross-lagged effect of latent Z on latent Y at the next time point (if applicable).
#' @param factor_loading_x The factor loadings for the indicators of X. The first loading is implicitly 1,
#'                         and subsequent loadings can be a single value or a vector of length
#'                         `num_indicators_x - 1`.
#' @param factor_loading_y The factor loadings for the indicators of Y. The first loading is implicitly 1,
#'                         and subsequent loadings can be a single value or a vector of length
#'                         `num_indicators_y - 1`.
#' @param factor_loading_z The factor loadings for the indicators of Z (if applicable). The first loading
#'                         is implicitly 1, and subsequent loadings can be a single value or a vector
#'                         of length `num_indicators_z - 1`.
#' @param residual_variance_x The residual variances for the indicators of X (can be a single value or a
#'                            vector of length `num_indicators_x`).
#' @param residual_variance_y The residual variances for the indicators of Y (can be a single value or a
#'                            vector of length `num_indicators_y`).
#' @param residual_variance_z The residual variances for the indicators of Z (if applicable, can be a
#'                            single value or a vector of length `num_indicators_z`).
#' @param variance.p The variance for the latent X variables.
#' @param variance.q The variance for the latent Y variables.
#' @param variance.r The variance for the latent Z variables (if applicable).
#' @param cov.pq The covariance between latent X and latent Y within the same time point, or wave
#' @param cov.pr The covariance between latent X and latent Z within the same time point (if specified by user).
#' @param cov.qr The covariance between latent Y and latent Z within the same time point (if specified by user).
#' @param variance.between.x The variance for the random intercept of X.
#' @param variance.between.y The variance for the random intercept of Y.
#' @param variance.between.z The variance for the random intercept of Z (if specified by user).
#' @param cov.between.xy The covariance of intercept terms between X and Y.
#' @param cov.between.xz The covariance of intercept terms between X and Z (if applicable).
#' @param cov.between.yz The covariance of intercept terms between Y and Z (if applicable).
#' @param ... Additional arguments to pass to the `lavaan::simulateData` function.
#'
#' @return A list containing two elements:
#'    * `model`: The lavaan model syntax used for the data simulation.
#'    * `data` : The simulated data in a data frame format.
#'
#' @export
simRICLPMl <- function(waves = 10,
                       num_indicators_x = 2,
                       num_indicators_y = 2,
                       num_indicators_z = NULL,
                       stability.p = 0.2,
                       stability.q = 0.2,
                       stability.r = NULL,
                       cross.p = 0.1,
                       cross.q = 0.1,
                       cross.zx = NULL,
                       cross.zy = NULL,
                       factor_loading_x = 1,
                       factor_loading_y = 1,
                       factor_loading_z = 1,
                       residual_variance_x = 0.5,
                       residual_variance_y = 0.5,
                       residual_variance_z = 0.5,
                       variance.p = 1,
                       variance.q = 1,
                       variance.r = NULL,
                       cov.pq = 0.1,
                       cov.pr = NULL,
                       cov.qr = NULL,
                       variance.between.x = 1,
                       variance.between.y = 1,
                       variance.between.z = NULL,
                       cov.between.xy = 0.5,
                       cov.between.xz = NULL,
                       cov.between.yz = NULL,
                       sample.size = 1000,
                       ...) {

  if (!is.numeric(waves) || waves <= 0 || waves != as.integer(waves)) {
    stop("❌ Error: 'waves' must be a positive integer.")
  }
  if (!is.numeric(num_indicators_x) || num_indicators_x <= 0 || num_indicators_x != as.integer(num_indicators_x)) {
    stop("❌ Error: 'num_indicators_x' must be a positive integer.")
  }
  if (!is.numeric(num_indicators_y) || num_indicators_y <= 0 || num_indicators_y != as.integer(num_indicators_y)) {
    stop("❌ Error: 'num_indicators_y' must be a positive integer.")
  }
  if (!is.null(num_indicators_z) && (!is.numeric(num_indicators_z) || num_indicators_z <= 0 || num_indicators_z != as.integer(num_indicators_z))) {
    stop("❌ Error: 'num_indicators_z' must be a positive integer.")
  }

  check_numeric_or_vector <- function(param, name, length = NULL, first_fixed = FALSE) {
    if (!is.null(param)) {
      if (!is.numeric(param)) {
        stop(paste0("❌ Error: '", name, "' must be numeric."))
      }
      if (!is.null(length)) {
        expected_length <- if (first_fixed) length - 1 else length
        if (length(param) != 1 && length(param) != expected_length) {
          stop(paste0("❌ Error: '", name, "' must be of length 1 or ", expected_length, "."))
        }
      }
    }
  }

  # Validate numeric/vector parameters with adjustment for first fixed loading
  check_numeric_or_vector(stability.p, "stability.p")
  check_numeric_or_vector(stability.q, "stability.q")
  check_numeric_or_vector(stability.r, "stability.r")
  check_numeric_or_vector(cross.p, "cross.p")
  check_numeric_or_vector(cross.q, "cross.q")
  check_numeric_or_vector(cross.zx, "cross.zx")
  check_numeric_or_vector(cross.zy, "cross.zy")
  check_numeric_or_vector(factor_loading_x, "factor_loading_x", num_indicators_x, first_fixed = TRUE)
  check_numeric_or_vector(factor_loading_y, "factor_loading_y", num_indicators_y, first_fixed = TRUE)
  check_numeric_or_vector(factor_loading_z, "factor_loading_z", num_indicators_z, first_fixed = TRUE)
  check_numeric_or_vector(residual_variance_x, "residual_variance_x", num_indicators_x)
  check_numeric_or_vector(residual_variance_y, "residual_variance_y", num_indicators_y)
  check_numeric_or_vector(residual_variance_z, "residual_variance_z", num_indicators_z)
  check_numeric_or_vector(variance.p, "variance.p")
  check_numeric_or_vector(variance.q, "variance.q")
  check_numeric_or_vector(variance.r, "variance.r")
  check_numeric_or_vector(cov.pq, "cov.pq")
  check_numeric_or_vector(cov.pr, "cov.pr")
  check_numeric_or_vector(cov.qr, "cov.qr")
  check_numeric_or_vector(variance.between.x, "variance.between.x")
  check_numeric_or_vector(variance.between.y, "variance.between.y")
  check_numeric_or_vector(variance.between.z, "variance.between.z")
  check_numeric_or_vector(cov.between.xy, "cov.between.xy")
  check_numeric_or_vector(cov.between.xz, "cov.between.xz")
  check_numeric_or_vector(cov.between.yz, "cov.between.yz")

  model_string <- ""

  # Random intercepts
  x_indicators_wave1 <- paste0("x", 1:num_indicators_x, "_1")
  model_string <- paste0(model_string, "\nBX =~ 1*", paste(x_indicators_wave1, collapse = " + 1*"))
  for (w in 2:waves) {
    x_indicators_w <- paste0("x", 1:num_indicators_x, "_", w)
    model_string <- paste0(model_string, " + 1*", paste(x_indicators_w, collapse = " + 1*"))
  }

  y_indicators_wave1 <- paste0("y", 1:num_indicators_y, "_1")
  model_string <- paste0(model_string, "\nBY =~ 1*", paste(y_indicators_wave1, collapse = " + 1*"))
  for (w in 2:waves) {
    y_indicators_w <- paste0("y", 1:num_indicators_y, "_", w)
    model_string <- paste0(model_string, " + 1*", paste(y_indicators_w, collapse = " + 1*"))
  }

  if (!is.null(num_indicators_z)) {
    z_indicators_wave1 <- paste0("z", 1:num_indicators_z, "_1")
    model_string <- paste0(model_string, "\nBZ =~ 1*", paste(z_indicators_wave1, collapse = " + 1*"))
    for (w in 2:waves) {
      z_indicators_w <- paste0("z", 1:num_indicators_z, "_", w)
      model_string <- paste0(model_string, " + 1*", paste(z_indicators_w, collapse = " + 1*"))
    }
  }

  # # Indicator means
  # for (w in 1:waves) {
  #   for (i in 1:num_indicators_x) {
  #     model_string <- paste0(model_string, "\nx", i, "_", w, "~ 1\n")
  #   }
  #   for (i in 1:num_indicators_y) {
  #     model_string <- paste0(model_string, "\ny", i, "_", w, "~ 1\n")
  #   }
  #   if (!is.null(num_indicators_z)) {
  #     for (i in 1:num_indicators_z) {
  #       model_string <- paste0(model_string, "\nz", i, "_", w, "~ 1\n")
  #     }
  #   }
  # }


  # Latent variable measurement model with first loading fixed and subsequent loadings specified
  for (w in 1:waves) {
    # X factor for wave w (pw)
    x_indicators_w <- paste0("x", 1:num_indicators_x, "_", w)
    model_string <- paste0(model_string, "\np", w, " =~ 1*", x_indicators_w[1]) # First loading fixed to 1
    if (num_indicators_x > 1) {
      fl_x <- if (length(factor_loading_x) == 1) rep(factor_loading_x, num_indicators_x - 1) else factor_loading_x
      for (i in 2:num_indicators_x) {
        model_string <- paste0(model_string, " + ", fl_x[i - 1], "*", x_indicators_w[i])
      }
    }
    model_string <- paste0(model_string, "\n")

    # Y factor for wave w (qw)
    y_indicators_w <- paste0("y", 1:num_indicators_y, "_", w)
    model_string <- paste0(model_string, "\nq", w, " =~ 1*", y_indicators_w[1]) # First loading fixed to 1
    if (num_indicators_y > 1) {
      fl_y <- if (length(factor_loading_y) == 1) rep(factor_loading_y, num_indicators_y - 1) else factor_loading_y
      for (i in 2:num_indicators_y) {
        model_string <- paste0(model_string, " + ", fl_y[i - 1], "*", y_indicators_w[i])
      }
    }
    model_string <- paste0(model_string, "\n")

    # Z factor for wave w (rw) - if applicable
    if (!is.null(num_indicators_z)) {
      z_indicators_w <- paste0("z", 1:num_indicators_z, "_", w)
      model_string <- paste0(model_string, "\nr", w, " =~ 1*", z_indicators_w[1]) # First loading fixed to 1
      if (num_indicators_z > 1) {
        fl_z <- if (length(factor_loading_z) == 1) rep(factor_loading_z, num_indicators_z - 1) else factor_loading_z
        for (i in 2:num_indicators_z) {
          model_string <- paste0(model_string, " + ", fl_z[i - 1], "*", z_indicators_w[i])
        }
      }
      model_string <- paste0(model_string, "\n")
    }
  }


  # Stability and cross-lagged paths
  for (w in 2:waves) {
    model_string <- paste0(
      model_string, "\n p", w, " ~ ", stability.p, " * p", w - 1, " + ", cross.q, " * q", w - 1
    )
    model_string <- paste0(
      model_string, "\n q", w, " ~ ", stability.q, " * q", w - 1, " + ", cross.p, " * p", w - 1
    )
    if (!is.null(num_indicators_z)) {
      model_string <- paste0(
        model_string, "\n p", w, " ~ ", cross.zx, " * r", w - 1,
        "\n q", w, " ~ ", cross.zy, " * r", w - 1,
        "\n r", w, " ~ ", stability.r, " * r", w - 1,
        " + ", cross.zx, " * p", w - 1, " + ", cross.zy, " * q", w - 1
      )
    }
  }

  # Latent variable variances and covariances
  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\n p", w, " ~~ ", variance.p, " * p", w,
      "\n q", w, " ~~ ", variance.q, " * q", w,
      "\n p", w, " ~~ ", cov.pq, " * q", w
    )
    if (!is.null(num_indicators_z)) {
      model_string <- paste0(
        model_string, "\n r", w, " ~~ ", variance.r, " * r", w,
        "\n p", w, " ~~ ", cov.pr, " * r", w,
        "\n q", w, " ~~ ", cov.qr, " * r", w
      )
    }
  }

  # Residual variances for indicators
  for (w in 1:waves) {
    for (i in 1:num_indicators_x) {
      res_var <- ifelse(length(residual_variance_x) == 1, residual_variance_x, residual_variance_x[i])
      model_string <- paste0(model_string, "\nx", i, "_", w, " ~~ ", res_var, "*x", i, "_", w, "\n")
    }
    for (i in 1:num_indicators_y) {
      res_var <- ifelse(length(residual_variance_y) == 1, residual_variance_y, residual_variance_y[i])
      model_string <- paste0(model_string, "\ny", i, "_", w, " ~~ ", res_var, "*y", i, "_", w, "\n")
    }
    if (!is.null(num_indicators_z)) {
      for (i in 1:num_indicators_z) {
        res_var <- ifelse(length(residual_variance_z) == 1, residual_variance_z, residual_variance_z[i])
        model_string <- paste0(model_string, "\nz", i, "_", w, " ~~ ", res_var, "*z", i, "_", w, "\n")
      }
    }
  }

  # Random intercept variances and covariances
  model_string <- paste0(
    model_string,
    "\nBX ~~ ", variance.between.x, " * BX",
    "\nBY ~~ ", variance.between.y, " * BY",
    "\nBX ~~ ", cov.between.xy, " * BY"
  )

  if (!is.null(num_indicators_z)) {
    model_string <- paste0(
      model_string,
      "\nBZ ~~ ", variance.between.z, " * BZ",
      "\nBX ~~ ", cov.between.xz, " * BZ",
      "\nBY ~~ ", cov.between.yz, " * BZ"
    )
  }

  message("You have simulated data for a Random Intercept Cross-Lagged Panel Model (RI-CLPM) with ",
          waves, " waves and ", num_indicators_x, " indicators X, ",
          num_indicators_y, " indicators Y",
          ifelse(!is.null(num_indicators_z), paste0(", and ", num_indicators_z, " Z indicators."), ","),
          " and a sample size of ", sample.size, "." )
  dat <- lavaan::simulateData(model = model_string,
                              int.ov.free = TRUE,
                              sample.nobs = sample.size,
                              ...)

  return(list(model = model_string, data = dat))
}

