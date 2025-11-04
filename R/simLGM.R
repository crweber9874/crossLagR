#' @title simLGM
#' @description Simulate data from a Linear Latent Growth Model (LGM) with options for univariate and bivariate models.
#'
#' @param waves The number of waves (time points) in the model.
#' @param variable_type Type of latent growth model: "univariate" (single variable) or "bivariate" (dual variable).
#' @param time_scores A numeric vector specifying the time scores for each wave. If NULL, uses 0, 1, 2, ..., waves-1.
#' @param mean_I_x Mean of the intercept factor for X.
#' @param mean_S_x Mean of the slope factor for X.
#' @param mean_Q_x Mean of the quadratic factor for X (if estimate_quadratic = TRUE).
#' @param mean_I_y Mean of the intercept factor for Y (bivariate only).
#' @param mean_S_y Mean of the slope factor for Y (bivariate only).
#' @param mean_Q_y Mean of the quadratic factor for Y (bivariate only).
#' @param var_I_x Variance of the intercept factor for X.
#' @param var_S_x Variance of the slope factor for X.
#' @param var_Q_x Variance of the quadratic factor for X (if estimate_quadratic = TRUE).
#' @param var_I_y Variance of the intercept factor for Y (bivariate only).
#' @param var_S_y Variance of the slope factor for Y (bivariate only).
#' @param var_Q_y Variance of the quadratic factor for Y (bivariate only).
#' @param cov_IS_x Covariance between intercept and slope factors for X.
#' @param cov_IQ_x Covariance between intercept and quadratic factors for X (if estimate_quadratic = TRUE).
#' @param cov_SQ_x Covariance between slope and quadratic factors for X (if estimate_quadratic = TRUE).
#' @param cov_IS_y Covariance between intercept and slope factors for Y (bivariate only).
#' @param cov_IQ_y Covariance between intercept and quadratic factors for Y (bivariate only).
#' @param cov_SQ_y Covariance between slope and quadratic factors for Y (bivariate only).
#' @param cov_I_xy Covariance between X and Y intercept factors (bivariate only).
#' @param cov_S_xy Covariance between X and Y slope factors (bivariate only).
#' @param cov_Q_xy Covariance between X and Y quadratic factors (bivariate only).
#' @param cov_IS_xy Covariance between X intercept and Y slope factors (bivariate only).
#' @param cov_IS_yx Covariance between Y intercept and X slope factors (bivariate only).
#' @param cov_IQ_xy Covariance between X intercept and Y quadratic factors (bivariate only).
#' @param cov_IQ_yx Covariance between Y intercept and X quadratic factors (bivariate only).
#' @param cov_SQ_xy Covariance between X slope and Y quadratic factors (bivariate only).
#' @param cov_SQ_yx Covariance between Y slope and X quadratic factors (bivariate only).
#' @param sigma2_x Unique variance for X latent variables (p factors).
#' @param sigma2_y Unique variance for Y latent variables (q factors, bivariate only).
#' @param cov_res Covariance between X and Y latent variables at each time point (bivariate only).
#' @param estimate_quadratic Logical. If TRUE, includes quadratic growth factors in simulation.
#' @param center_time Logical. If TRUE, centers time scores around their mean.
#' @param sample_size The number of observations to simulate. Default is 1000.
#' @param seed Random seed for reproducibility. Default is NULL.
#' @param ... Additional arguments to pass to the `lavaan::simulateData` function.
#'
#' @return A list containing two elements:
#'    * `model`: The Lavaan model syntax used for data simulation.
#'    * `data`: The simulated data in a data frame format.
#'
#' @details
#' This function simulates data from a Linear Latent Growth Model (LGM) where:
#' - Growth is modeled using latent intercept, slope, and optional quadratic factors
#' - Individual differences exist in growth parameters
#' - Systematic change occurs over time according to specified growth functions
#'
#' For univariate models:
#' - Single variable measured across time points
#' - Intercept factor determines initial level
#' - Slope factor determines rate of linear change
#' - Optional quadratic factor for acceleration/deceleration
#'
#' For bivariate models:
#' - Two variables measured across time points
#' - Growth factors for both variables with correlations between them
#' - Residual covariances between variables at each time point
#'
#' Parameter interpretation:
#' - mean_I: Average initial level across individuals
#' - mean_S: Average rate of linear change
#' - mean_Q: Average acceleration/deceleration (quadratic models)
#' - var_I: Individual differences in initial levels
#' - var_S: Individual differences in rates of change
#' - var_Q: Individual differences in acceleration
#' - cov_IS: Relationship between initial level and rate of change
#' - sigma2: Unique variance for latent variables (measurement precision)
#' - cov_res: Covariance between latent variables (bivariate models)
#'
#' Time scoring options:
#' - Default: 0, 1, 2, ..., waves-1 (intercept = initial observation)
#' - Custom: User-specified time points for unequal spacing
#' - Centered: Time scores centered around their mean
#'
#' @examples
#' # Basic univariate LGM
#' lgm_data <- simLGM(
#'   waves = 5,
#'   variable_type = "univariate",
#'   mean_I_x = 10,
#'   mean_S_x = 0.5,
#'   var_I_x = 4,
#'   var_S_x = 1,
#'   cov_IS_x = -0.5,
#'   sigma2_x = 2,
#'   sample_size = 500
#' )
#'
#' # Bivariate LGM with correlation between variables
#' bivar_lgm <- simLGM(
#'   waves = 4,
#'   variable_type = "bivariate",
#'   mean_I_x = 10, mean_I_y = 15,
#'   mean_S_x = 0.3, mean_S_y = -0.2,
#'   var_I_x = 2, var_I_y = 3,
#'   var_S_x = 0.5, var_S_y = 0.8,
#'   cov_I_xy = 1.2,
#'   cov_S_xy = 0.3,
#'   sample_size = 800
#' )
#'
#' # Quadratic growth model
#' quad_lgm <- simLGM(
#'   waves = 6,
#'   estimate_quadratic = TRUE,
#'   mean_I_x = 20,
#'   mean_S_x = 2,
#'   mean_Q_x = -0.1,
#'   var_I_x = 9,
#'   var_S_x = 1,
#'   var_Q_x = 0.01,
#'   sample_size = 1200
#' )
#'
#' @export

simLGM <- function(waves = 5,
                   variable_type = c("univariate", "bivariate"),
                   time_scores = NULL,
                   mean_I_x = 10,
                   mean_S_x = 0.5,
                   mean_Q_x = 0,
                   mean_I_y = 10,
                   mean_S_y = 0.5,
                   mean_Q_y = 0,
                   var_I_x = 4,
                   var_S_x = 1,
                   var_Q_x = 0.01,
                   var_I_y = 4,
                   var_S_y = 1,
                   var_Q_y = 0.01,
                   cov_IS_x = 0,
                   cov_IQ_x = 0,
                   cov_SQ_x = 0,
                   cov_IS_y = 0,
                   cov_IQ_y = 0,
                   cov_SQ_y = 0,
                   cov_I_xy = 0,
                   cov_S_xy = 0,
                   cov_Q_xy = 0,
                   cov_IS_xy = 0,
                   cov_IS_yx = 0,
                   cov_IQ_xy = 0,
                   cov_IQ_yx = 0,
                   cov_SQ_xy = 0,
                   cov_SQ_yx = 0,
                   sigma2_x = 1,
                   sigma2_y = 1,
                   cov_res = 0,
                   estimate_quadratic = FALSE,
                   center_time = FALSE,
                   sample_size = 1000,
                   seed = NULL,
                   ...) {

  # Input validation
  if (!is.numeric(waves) || waves <= 0 || waves != as.integer(waves)) {
    stop("Error: Parameter 'waves' must be a positive integer.")
  }

  if (waves < 3) {
    stop("Error: LGM simulation requires at least 3 waves.")
  }

  if (estimate_quadratic && waves < 4) {
    stop("Error: Quadratic LGM simulation requires at least 4 waves.")
  }

  variable_type <- match.arg(variable_type)

  # Validate numeric parameters
  numeric_params <- list(
    mean_I_x = mean_I_x, mean_S_x = mean_S_x, mean_Q_x = mean_Q_x,
    var_I_x = var_I_x, var_S_x = var_S_x, var_Q_x = var_Q_x,
    cov_IS_x = cov_IS_x, cov_IQ_x = cov_IQ_x, cov_SQ_x = cov_SQ_x,
    sigma2_x = sigma2_x, sample_size = sample_size
  )

  if (variable_type == "bivariate") {
    bivariate_params <- list(
      mean_I_y = mean_I_y, mean_S_y = mean_S_y, mean_Q_y = mean_Q_y,
      var_I_y = var_I_y, var_S_y = var_S_y, var_Q_y = var_Q_y,
      cov_IS_y = cov_IS_y, cov_IQ_y = cov_IQ_y, cov_SQ_y = cov_SQ_y,
      cov_I_xy = cov_I_xy, cov_S_xy = cov_S_xy, cov_Q_xy = cov_Q_xy,
      cov_IS_xy = cov_IS_xy, cov_IS_yx = cov_IS_yx,
      cov_IQ_xy = cov_IQ_xy, cov_IQ_yx = cov_IQ_yx,
      cov_SQ_xy = cov_SQ_xy, cov_SQ_yx = cov_SQ_yx,
      sigma2_y = sigma2_y, cov_res = cov_res
    )
    numeric_params <- c(numeric_params, bivariate_params)
  }

  for (param_name in names(numeric_params)) {
    if (!is.numeric(numeric_params[[param_name]]) || length(numeric_params[[param_name]]) != 1) {
      stop(paste0("Error: Parameter '", param_name, "' must be a single numeric value."))
    }
  }

  # Check for positive variances
  if (var_I_x <= 0 || var_S_x <= 0 || sigma2_x <= 0) {
    stop("Error: Variance parameters for X must be positive.")
  }

  if (estimate_quadratic && var_Q_x <= 0) {
    stop("Error: Quadratic variance parameter for X must be positive.")
  }

  if (variable_type == "bivariate") {
    if (var_I_y <= 0 || var_S_y <= 0 || sigma2_y <= 0) {
      stop("Error: Variance parameters for Y must be positive.")
    }
    if (estimate_quadratic && var_Q_y <= 0) {
      stop("Error: Quadratic variance parameter for Y must be positive.")
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

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # ==================== MODEL GENERATION ====================
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
    model_string <- paste0(model_string, "\n# Growth Factor Means\n")
    model_string <- paste0(model_string, "I ~ ", mean_I_x, "*1\n")
    model_string <- paste0(model_string, "S ~ ", mean_S_x, "*1\n")
    if (estimate_quadratic) {
      model_string <- paste0(model_string, "Q ~ ", mean_Q_x, "*1\n")
    }

    # Growth factor variances and covariances
    model_string <- paste0(model_string, "\n# Growth Factor Variances and Covariances\n")
    model_string <- paste0(model_string, "I ~~ ", var_I_x, "*I\n")
    model_string <- paste0(model_string, "S ~~ ", var_S_x, "*S\n")
    model_string <- paste0(model_string, "I ~~ ", cov_IS_x, "*S\n")

    if (estimate_quadratic) {
      model_string <- paste0(model_string, "Q ~~ ", var_Q_x, "*Q\n")
      model_string <- paste0(model_string, "I ~~ ", cov_IQ_x, "*Q\n")
      model_string <- paste0(model_string, "S ~~ ", cov_SQ_x, "*Q\n")
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
    for (w in 1:waves) {
      model_string <- paste0(model_string, "p", w, " ~~ ", sigma2_x, "*p", w, "\n")
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
    model_string <- paste0(model_string, "\n# Growth Factor Means\n")
    model_string <- paste0(model_string, "I_x ~ ", mean_I_x, "*1\n")
    model_string <- paste0(model_string, "S_x ~ ", mean_S_x, "*1\n")
    model_string <- paste0(model_string, "I_y ~ ", mean_I_y, "*1\n")
    model_string <- paste0(model_string, "S_y ~ ", mean_S_y, "*1\n")

    if (estimate_quadratic) {
      model_string <- paste0(model_string, "Q_x ~ ", mean_Q_x, "*1\n")
      model_string <- paste0(model_string, "Q_y ~ ", mean_Q_y, "*1\n")
    }

    # Growth factor variances
    model_string <- paste0(model_string, "\n# Growth Factor Variances\n")
    model_string <- paste0(model_string, "I_x ~~ ", var_I_x, "*I_x\n")
    model_string <- paste0(model_string, "S_x ~~ ", var_S_x, "*S_x\n")
    model_string <- paste0(model_string, "I_y ~~ ", var_I_y, "*I_y\n")
    model_string <- paste0(model_string, "S_y ~~ ", var_S_y, "*S_y\n")

    if (estimate_quadratic) {
      model_string <- paste0(model_string, "Q_x ~~ ", var_Q_x, "*Q_x\n")
      model_string <- paste0(model_string, "Q_y ~~ ", var_Q_y, "*Q_y\n")
    }

    # Within-variable covariances
    model_string <- paste0(model_string, "\n# Within-Variable Growth Factor Covariances\n")
    model_string <- paste0(model_string, "I_x ~~ ", cov_IS_x, "*S_x\n")
    model_string <- paste0(model_string, "I_y ~~ ", cov_IS_y, "*S_y\n")

    if (estimate_quadratic) {
      model_string <- paste0(model_string, "I_x ~~ ", cov_IQ_x, "*Q_x\n")
      model_string <- paste0(model_string, "S_x ~~ ", cov_SQ_x, "*Q_x\n")
      model_string <- paste0(model_string, "I_y ~~ ", cov_IQ_y, "*Q_y\n")
      model_string <- paste0(model_string, "S_y ~~ ", cov_SQ_y, "*Q_y\n")
    }

    # Between-variable covariances
    model_string <- paste0(model_string, "\n# Between-Variable Growth Factor Covariances\n")
    model_string <- paste0(model_string, "I_x ~~ ", cov_I_xy, "*I_y\n")
    model_string <- paste0(model_string, "S_x ~~ ", cov_S_xy, "*S_y\n")
    model_string <- paste0(model_string, "I_x ~~ ", cov_IS_xy, "*S_y\n")
    model_string <- paste0(model_string, "I_y ~~ ", cov_IS_yx, "*S_x\n")

    if (estimate_quadratic) {
      model_string <- paste0(model_string, "Q_x ~~ ", cov_Q_xy, "*Q_y\n")
      model_string <- paste0(model_string, "I_x ~~ ", cov_IQ_xy, "*Q_y\n")
      model_string <- paste0(model_string, "I_y ~~ ", cov_IQ_yx, "*Q_x\n")
      model_string <- paste0(model_string, "S_x ~~ ", cov_SQ_xy, "*Q_y\n")
      model_string <- paste0(model_string, "S_y ~~ ", cov_SQ_yx, "*Q_x\n")
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
    for (w in 1:waves) {
      model_string <- paste0(model_string, "p", w, " ~~ ", sigma2_x, "*p", w, "\n")
      model_string <- paste0(model_string, "q", w, " ~~ ", sigma2_y, "*q", w, "\n")
    }

    # Latent variable covariances
    model_string <- paste0(model_string, "\n# Latent Variable Covariances\n")
    for (w in 1:waves) {
      model_string <- paste0(model_string, "p", w, " ~~ ", cov_res, "*q", w, "\n")
    }

    # Fix observed variable residuals to 0 (perfect indicators)
    model_string <- paste0(model_string, "\n# Fix Observed Variable Residuals to Zero\n")
    for (w in 1:waves) {
      model_string <- paste0(model_string, "x", w, " ~~ 0*x", w, "\n")
      model_string <- paste0(model_string, "y", w, " ~~ 0*y", w, "\n")
    }
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
                "\nCheck parameter values for model identification issues.",
                "\nEnsure variance-covariance matrices are positive definite."))
  })

  # Validation of generated data
  if (nrow(dat) != sample_size) {
    warning("Generated sample size differs from requested sample size.")
  }

  return(list(model = model_string, data = dat))
}
