#' @title simCLPM
#' @description Simulate data from a Cross-Lagged Panel Model (CLPM) with options for parameter constraints and starting values.
#'
#' @param waves The number of waves (time points) in the model.
#' @param beta_y The autoregressive effect for the y variable (stability parameter).
#' @param beta_x The autoregressive effect for the x variable (stability parameter).
#' @param omega_xy The cross-lagged effect of x on y at the next time point.
#' @param omega_yx The cross-lagged effect of y on x at the next time point.
#' @param var_y The residual variance for the y latent variables.
#' @param var_x The residual variance for the x latent variables.
#' @param cov_xy The covariance between x and y within the same time point.
#' @param var_y1 The variance for the first wave y variable. If NULL, uses var_y.
#' @param var_x1 The variance for the first wave x variable. If NULL, uses var_x.
#' @param cov_xy1 The covariance for the first wave between x and y. If NULL, uses cov_xy.
#' @param mean_y1 The mean for the first wave y variable. Default is 0.
#' @param mean_x1 The mean for the first wave x variable. Default is 0.
#' @param sample_size The number of observations to simulate. Default is 1000.
#' @param seed Random seed for reproducibility. Default is NULL.
#' @param ... Additional arguments to pass to the `lavaan::simulateData` function.
#'
#' @return A list containing two elements:
#'    * `model`: The Lavaan model syntax used for data simulation.
#'    * `data`: The simulated data in a data frame format.
#'
#' @details
#' This function simulates data from a Cross-Lagged Panel Model (CLPM) where:
#' - p variables represent latent y constructs
#' - q variables represent latent x constructs
#' - beta parameters represent autoregressive (stability) effects
#' - omega parameters represent cross-lagged (reciprocal) effects
#'
#' The CLPM assumes all relationships are due to within-person processes
#' and does not separate stable between-person differences from dynamic effects
#' (unlike the RI-CLPM).
#'
#' Parameter interpretation:
#' - beta_y: How much y at time t predicts y at time t+1
#' - beta_x: How much x at time t predicts x at time t+1
#' - omega_xy: How much x at time t predicts y at time t+1
#' - omega_yx: How much y at time t predicts x at time t+1
#'
#' @examples
#' # Basic CLPM simulation
#' sim_result <- simCLPM(
#'   waves = 5,
#'   beta_y = 0.4,
#'   beta_x = 0.3,
#'   omega_xy = 0.2,
#'   omega_yx = 0.1,
#'   sample_size = 500
#' )
#'
#' # View the generated model
#' cat(sim_result$model)
#'
#' # View the data
#' head(sim_result$data)
#'
#' @export

simCLPM <- function(waves = 5,
                    beta_y = 0.4,
                    beta_x = 0.4,
                    omega_xy = 0.2,
                    omega_yx = 0.2,
                    var_y = 1,
                    var_x = 1,
                    cov_xy = 0.3,
                    var_y1 = NULL,
                    var_x1 = NULL,
                    cov_xy1 = NULL,
                    mean_y1 = 0,
                    mean_x1 = 0,
                    sample_size = 1000,
                    seed = NULL,
                    ...) {

  # Input validation
  if (!is.numeric(waves) || waves <= 0 || waves != as.integer(waves)) {
    stop("Error: Parameter 'waves' must be a positive integer.")
  }

  if (waves < 2) {
    stop("Error: CLPM simulation requires at least 2 waves.")
  }

  # Validate numeric parameters
  numeric_params <- list(
    beta_y = beta_y, beta_x = beta_x, omega_xy = omega_xy, omega_yx = omega_yx,
    var_y = var_y, var_x = var_x, cov_xy = cov_xy,
    mean_y1 = mean_y1, mean_x1 = mean_x1, sample_size = sample_size
  )

  for (param_name in names(numeric_params)) {
    if (!is.numeric(numeric_params[[param_name]]) || length(numeric_params[[param_name]]) != 1) {
      stop(paste0("Error: Parameter '", param_name, "' must be a single numeric value."))
    }
  }

  # Set default values for wave 1 parameters if not specified
  if (is.null(var_y1)) var_y1 <- var_y
  if (is.null(var_x1)) var_x1 <- var_x
  if (is.null(cov_xy1)) cov_xy1 <- cov_xy

  # Validate wave 1 parameters
  if (!is.numeric(var_y1) || !is.numeric(var_x1) || !is.numeric(cov_xy1)) {
    stop("Error: Wave 1 variance and covariance parameters must be numeric.")
  }

  # Check for positive variances
  if (var_y <= 0 || var_x <= 0 || var_y1 <= 0 || var_x1 <= 0) {
    stop("Error: Variance parameters must be positive.")
  }

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # ==================== MODEL GENERATION ====================
  model_string <- ""

  # Measurement model (perfect indicators)
  for (w in 1:waves) {
    model_string <- paste0(
      model_string,
      "    p", w, " =~ 1*y", w, "\n",
      "    q", w, " =~ 1*x", w, "\n"
    )
  }

  # Means for observed variables
  for (w in 1:waves) {
    model_string <- paste0(
      model_string,
      "    y", w, " ~ 1\n",
      "    x", w, " ~ 1\n"
    )
  }

  # Latent variable means
  # First wave means (can be non-zero)
  model_string <- paste0(
    model_string,
    "    p1 ~ ", mean_y1, "*1\n",
    "    q1 ~ ", mean_x1, "*1\n"
  )

  # Subsequent wave means (fixed to 0 due to regressions)
  for (w in 2:waves) {
    model_string <- paste0(
      model_string,
      "    p", w, " ~ 0*1\n",
      "    q", w, " ~ 0*1\n"
    )
  }

  # Autoregressive and cross-lagged effects
  for (w in 2:waves) {
    model_string <- paste0(
      model_string,
      "    p", w, " ~ ", beta_y, "*p", w - 1, " + ", omega_xy, "*q", w - 1, "\n",
      "    q", w, " ~ ", beta_x, "*q", w - 1, " + ", omega_yx, "*p", w - 1, "\n"
    )
  }

  # Variances and covariances
  # First wave (may have different parameters)
  model_string <- paste0(
    model_string,
    "    p1 ~~ ", var_y1, "*p1\n",
    "    q1 ~~ ", var_x1, "*q1\n",
    "    p1 ~~ ", cov_xy1, "*q1\n"
  )

  # Subsequent waves (residual variances and covariances)
  for (w in 2:waves) {
    model_string <- paste0(
      model_string,
      "    p", w, " ~~ ", var_y, "*p", w, "\n",
      "    q", w, " ~~ ", var_x, "*q", w, "\n",
      "    p", w, " ~~ ", cov_xy, "*q", w, "\n"
    )
  }

  # Fix observed variable residual variances to 0 (perfect indicators)
  for (w in 1:waves) {
    model_string <- paste0(
      model_string,
      "    y", w, " ~~ 0*y", w, "\n",
      "    x", w, " ~~ 0*x", w, "\n"
    )
  }

  # Generate data
  tryCatch({
    dat <- lavaan::simulateData(
      model = model_string,
      sample.nobs = sample_size,
      int.ov.free = TRUE,
      ...
    )
  }, error = function(e) {
    stop(paste0("Error in data generation: ", e$message,
                "\nCheck parameter values for model identification issues."))
  })

  # Validation of generated data
  if (nrow(dat) != sample_size) {
    warning("Generated sample size differs from requested sample size.")
  }

  return(list(model = model_string, data = dat))
}
