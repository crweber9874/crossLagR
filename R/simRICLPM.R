#' @title simRICLPM
#' @description Simulate data from a Random Intercept Cross-Lagged Panel Model (RICLPM) with proper intercept structure
#'
#' This function creates synthetic data for a Random Intercept Cross-Lagged Panel Model (RICLPM)
#' that matches the parameter structure and identification constraints of the estimateRICLPM function.
#' The model separates stable between-person differences (random intercepts) from within-person
#' dynamics (autoregressive and cross-lagged effects).
#'
#' @param waves The number of waves (time points) in the model. Must be at least 2.
#' @param sample.nobs The number of observations (participants) to simulate. Default is 1000.
#' @param include_z Logical. If TRUE, includes a third variable Z in the model. Default is FALSE.
#' @param beta_x Autoregressive effect for X (X(t-1) -> X(t)). Default is 0.3.
#' @param beta_y Autoregressive effect for Y (Y(t-1) -> Y(t)). Default is 0.3.
#' @param beta_z Autoregressive effect for Z (Z(t-1) -> Z(t)), if included. Default is 0.3.
#' @param omega_xy Cross-lagged effect of Y on X (Y(t-1) -> X(t)). Default is 0.1.
#' @param omega_yx Cross-lagged effect of X on Y (X(t-1) -> Y(t)). Default is 0.1.
#' @param omega_zx Cross-lagged effect of Z on X (Z(t-1) -> X(t)), if Z included. Default is 0.1.
#' @param omega_zy Cross-lagged effect of Z on Y (Z(t-1) -> Y(t)), if Z included. Default is 0.1.
#' @param omega_xz Cross-lagged effect of X on Z (X(t-1) -> Z(t)), if Z included. Default is 0.1.
#' @param omega_yz Cross-lagged effect of Y on Z (Y(t-1) -> Z(t)), if Z included. Default is 0.1.
#' @param var_p Residual variance for within-person X factors. Default is 1.
#' @param var_q Residual variance for within-person Y factors. Default is 1.
#' @param var_r Residual variance for within-person Z factors, if included. Default is 1.
#' @param cov_pq Residual covariance between within-person X and Y factors. Default is 0.1.
#' @param cov_pr Residual covariance between within-person X and Z factors, if Z included. Default is 0.1.
#' @param cov_qr Residual covariance between within-person Y and Z factors, if Z included. Default is 0.1.
#' @param var_BX Variance of random intercept for X. Default is 1.5.
#' @param var_BY Variance of random intercept for Y. Default is 1.5.
#' @param var_BZ Variance of random intercept for Z, if included. Default is 1.5.
#' @param cov_BXBY Covariance between random intercepts for X and Y. Default is 0.5.
#' @param cov_BXBZ Covariance between random intercepts for X and Z, if Z included. Default is 0.5.
#' @param cov_BYBZ Covariance between random intercepts for Y and Z, if Z included. Default is 0.5.
#' @param mean_BX Mean of random intercept for X. Default is 0.
#' @param mean_BY Mean of random intercept for Y. Default is 0.
#' @param mean_BZ Mean of random intercept for Z, if included. Default is 0.
#' @param constrain_beta Logical. If TRUE, constrains autoregressive effects to equality across waves. Default is TRUE.
#' @param constrain_omega Logical. If TRUE, constrains cross-lagged effects to equality across waves. Default is TRUE.
#' @param constrain_residual_variances Logical. If TRUE, constrains residual variances to equality across waves. Default is TRUE.
#' @param constrain_residual_covariances Logical. If TRUE, constrains residual covariances to equality across waves. Default is TRUE.
#' @param estimate_means Logical. If TRUE, estimates means for within-person factors. Default is FALSE.
#' @param ... Additional arguments to pass to lavaan::simulateData.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{model}: The Lavaan model syntax used for data simulation
#'     \item \code{data}: The simulated data in a data frame format
#'     \item \code{parameters}: A list of the parameters used in simulation
#'   }
#'
#' @details
#' The Random Intercept Cross-Lagged Panel Model (RICLPM) decomposes each observed variable into:
#' - A stable between-person component (random intercept)
#' - A time-varying within-person component (within-person factor)
#' - Measurement error (fixed to 0 for perfect indicators)
#'
#' Variable naming convention:
#' - X variables: x1, x2, x3, ..., x[waves]
#' - Y variables: y1, y2, y3, ..., y[waves]
#' - Z variables (if included): z1, z2, z3, ..., z[waves]
#'
#' Model identification follows standard RICLPM practices:
#' - Random intercept means are estimated (unless estimate_means = FALSE)
#' - Within-person factor means are typically fixed to 0 (representing deviations)
#' - Random intercepts don't correlate with first wave within-person factors
#' - Observed variable residuals are fixed to 0 (perfect indicators)
#'
#' @examples
#' # Basic RICLPM simulation
#' riclpm_data <- simRICLPM(waves = 4, sample.nobs = 500)
#'
#' # RICLPM with stronger cross-lagged effects
#' riclpm_strong <- simRICLPM(
#'   waves = 5,
#'   omega_xy = 0.3,
#'   omega_yx = 0.3,
#'   sample.nobs = 1000
#' )
#'
#' # RICLPM with third variable
#' riclpm_3var <- simRICLPM(
#'   waves = 4,
#'   include_z = TRUE,
#'   sample.nobs = 800
#' )
#'
#' # RICLPM with unconstrained parameters
#' riclpm_free <- simRICLPM(
#'   waves = 4,
#'   constrain_beta = FALSE,
#'   constrain_omega = FALSE,
#'   sample.nobs = 600
#' )
#'
#' @export
simRICLPM <- function(waves = 5,
                      sample.nobs = 1000,
                      include_z = FALSE,
                      beta_x = 0.3,
                      beta_y = 0.3,
                      beta_z = 0.3,
                      omega_xy = 0.1,
                      omega_yx = 0.1,
                      omega_zx = 0.1,
                      omega_zy = 0.1,
                      omega_xz = 0.1,
                      omega_yz = 0.1,
                      var_p = 1,
                      var_q = 1,
                      var_r = 1,
                      cov_pq = 0.1,
                      cov_pr = 0.1,
                      cov_qr = 0.1,
                      var_BX = 1.5,
                      var_BY = 1.5,
                      var_BZ = 1.5,
                      cov_BXBY = 0.5,
                      cov_BXBZ = 0.5,
                      cov_BYBZ = 0.5,
                      mean_BX = 0,
                      mean_BY = 0,
                      mean_BZ = 0,
                      constrain_beta = TRUE,
                      constrain_omega = TRUE,
                      constrain_residual_variances = TRUE,
                      constrain_residual_covariances = TRUE,
                      estimate_means = FALSE,
                      ...) {

  # Input validation
  if (!is.numeric(waves) || waves <= 0 || waves != as.integer(waves)) {
    stop("❌ Error: Parameter 'waves' must be a positive integer.")
  }

  if (waves < 2) {
    stop("❌ Error: RICLPM requires at least 2 waves.")
  }

  if (!is.numeric(sample.nobs) || sample.nobs <= 0 || sample.nobs != as.integer(sample.nobs)) {
    stop("❌ Error: Parameter 'sample.nobs' must be a positive integer.")
  }

  if (!is.logical(include_z)) {
    stop("❌ Error: 'include_z' must be logical (TRUE/FALSE).")
  }

  # Validate logical parameters
  logical_params <- list(
    constrain_beta = constrain_beta,
    constrain_omega = constrain_omega,
    constrain_residual_variances = constrain_residual_variances,
    constrain_residual_covariances = constrain_residual_covariances,
    estimate_means = estimate_means
  )

  for (param_name in names(logical_params)) {
    if (!is.logical(logical_params[[param_name]])) {
      stop(paste0("❌ Error: Parameter '", param_name, "' must be logical (TRUE/FALSE)."))
    }
  }

  # Generate hard-coded variable names
  time_varying_x <- paste0("x", 1:waves)
  time_varying_y <- paste0("y", 1:waves)
  time_varying_z <- if (include_z) paste0("z", 1:waves) else NULL

  # Build model string
  model_string <- ""

  # ==================== RANDOM INTERCEPTS ====================
  model_string <- paste0(model_string, "\n# Random Intercepts\nBX =~ 1*", time_varying_x[1])
  for (w in 2:waves) {
    model_string <- paste0(model_string, " + 1*", time_varying_x[w])
  }

  model_string <- paste0(model_string, "\nBY =~ 1*", time_varying_y[1])
  for (w in 2:waves) {
    model_string <- paste0(model_string, " + 1*", time_varying_y[w])
  }

  if (include_z) {
    model_string <- paste0(model_string, "\nBZ =~ 1*", time_varying_z[1])
    for (w in 2:waves) {
      model_string <- paste0(model_string, " + 1*", time_varying_z[w])
    }
  }

  # ==================== WITHIN-PERSON LATENT VARIABLES ====================
  model_string <- paste0(model_string, "\n\n# Within-Person Latent Variables")
  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\np", w, " =~ 1*", time_varying_x[w],
      "\nq", w, " =~ 1*", time_varying_y[w]
    )
    if (include_z) {
      model_string <- paste0(
        model_string, "\nr", w, " =~ 1*", time_varying_z[w]
      )
    }
  }

  # ==================== MEANS ====================
  if (estimate_means) {
    model_string <- paste0(model_string, "\n\n# Latent Variable Means")
    # Random intercept means (estimated)
    model_string <- paste0(model_string, "\nBX ~ ", mean_BX, "*1")
    model_string <- paste0(model_string, "\nBY ~ ", mean_BY, "*1")
    if (include_z) {
      model_string <- paste0(model_string, "\nBZ ~ ", mean_BZ, "*1")
    }

    # Within-person latent variable means (estimated)
    for (w in 1:waves) {
      model_string <- paste0(model_string, "\np", w, " ~ 1")
      model_string <- paste0(model_string, "\nq", w, " ~ 1")
      if (include_z) {
        model_string <- paste0(model_string, "\nr", w, " ~ 1")
      }
    }
  } else {
    model_string <- paste0(model_string, "\n\n# Latent Variable Means")
    # Random intercept means (estimated even when estimate_means = FALSE)
    model_string <- paste0(model_string, "\nBX ~ ", mean_BX, "*1")
    model_string <- paste0(model_string, "\nBY ~ ", mean_BY, "*1")
    if (include_z) {
      model_string <- paste0(model_string, "\nBZ ~ ", mean_BZ, "*1")
    }

    # Within-person latent variable means (fixed to 0)
    for (w in 1:waves) {
      model_string <- paste0(model_string, "\np", w, " ~ 0*1")
      model_string <- paste0(model_string, "\nq", w, " ~ 0*1")
      if (include_z) {
        model_string <- paste0(model_string, "\nr", w, " ~ 0*1")
      }
    }
  }

  # Observed variable intercepts (always fixed to 0)
  model_string <- paste0(model_string, "\n\n# Observed Variable Intercepts (Fixed to 0)")
  for (w in 1:waves) {
    model_string <- paste0(model_string, "\n", time_varying_x[w], " ~ 0*1")
    model_string <- paste0(model_string, "\n", time_varying_y[w], " ~ 0*1")
    if (include_z) {
      model_string <- paste0(model_string, "\n", time_varying_z[w], " ~ 0*1")
    }
  }

  # ==================== AUTOREGRESSIVE AND CROSS-LAGGED PATHS ====================
  model_string <- paste0(model_string, "\n\n# Autoregressive and Cross-lagged Effects")

  for (w in 2:waves) {
    if (constrain_beta && constrain_omega) {
      # Both autoregressive and cross-lagged effects constrained
      model_string <- paste0(
        model_string, "\np", w, " ~ ", beta_y, "*p", w - 1, " + ", omega_xy, "*q", w - 1
      )
      model_string <- paste0(
        model_string, "\nq", w, " ~ ", beta_x, "*q", w - 1, " + ", omega_yx, "*p", w - 1
      )

      if (include_z) {
        model_string <- paste0(
          model_string, "\nr", w, " ~ ", beta_z, "*r", w - 1,
          " + ", omega_zx, "*p", w - 1, " + ", omega_zy, "*q", w - 1
        )
        model_string <- paste0(
          model_string, "\np", w, " ~ ", omega_xz, "*r", w - 1
        )
        model_string <- paste0(
          model_string, "\nq", w, " ~ ", omega_yz, "*r", w - 1
        )
      }
    } else if (constrain_beta && !constrain_omega) {
      # Only autoregressive effects constrained
      model_string <- paste0(
        model_string, "\np", w, " ~ ", beta_y, "*p", w - 1, " + omega_xy", w, "*q", w - 1
      )
      model_string <- paste0(
        model_string, "\nq", w, " ~ ", beta_x, "*q", w - 1, " + omega_yx", w, "*p", w - 1
      )

      if (include_z) {
        model_string <- paste0(
          model_string, "\nr", w, " ~ ", beta_z, "*r", w - 1,
          " + omega_zx", w, "*p", w - 1, " + omega_zy", w, "*q", w - 1
        )
        model_string <- paste0(
          model_string, "\np", w, " ~ omega_xz", w, "*r", w - 1
        )
        model_string <- paste0(
          model_string, "\nq", w, " ~ omega_yz", w, "*r", w - 1
        )
      }
    } else if (!constrain_beta && constrain_omega) {
      # Only cross-lagged effects constrained
      model_string <- paste0(
        model_string, "\np", w, " ~ beta_y", w, "*p", w - 1, " + ", omega_xy, "*q", w - 1
      )
      model_string <- paste0(
        model_string, "\nq", w, " ~ beta_x", w, "*q", w - 1, " + ", omega_yx, "*p", w - 1
      )

      if (include_z) {
        model_string <- paste0(
          model_string, "\nr", w, " ~ beta_z", w, "*r", w - 1,
          " + ", omega_zx, "*p", w - 1, " + ", omega_zy, "*q", w - 1
        )
        model_string <- paste0(
          model_string, "\np", w, " ~ ", omega_xz, "*r", w - 1
        )
        model_string <- paste0(
          model_string, "\nq", w, " ~ ", omega_yz, "*r", w - 1
        )
      }
    } else {
      # No constraints - all effects have unique labels per wave
      model_string <- paste0(
        model_string, "\np", w, " ~ beta_y", w, "*p", w - 1, " + omega_xy", w, "*q", w - 1
      )
      model_string <- paste0(
        model_string, "\nq", w, " ~ beta_x", w, "*q", w - 1, " + omega_yx", w, "*p", w - 1
      )

      if (include_z) {
        model_string <- paste0(
          model_string, "\nr", w, " ~ beta_z", w, "*r", w - 1,
          " + omega_zx", w, "*p", w - 1, " + omega_zy", w, "*q", w - 1
        )
        model_string <- paste0(
          model_string, "\np", w, " ~ omega_xz", w, "*r", w - 1
        )
        model_string <- paste0(
          model_string, "\nq", w, " ~ omega_yz", w, "*r", w - 1
        )
      }
    }
  }

  # ==================== RESIDUAL VARIANCES ====================
  model_string <- paste0(model_string, "\n\n# Residual Variances")
  if (constrain_residual_variances) {
    for (w in 1:waves) {
      model_string <- paste0(
        model_string, "\np", w, " ~~ ", var_p, "*p", w,
        "\nq", w, " ~~ ", var_q, "*q", w
      )
      if (include_z) {
        model_string <- paste0(
          model_string, "\nr", w, " ~~ ", var_r, "*r", w
        )
      }
    }
  } else {
    # When unconstrained, each wave gets its own variance parameter
    for (w in 1:waves) {
      model_string <- paste0(
        model_string, "\np", w, " ~~ ", var_p, "*p", w,
        "\nq", w, " ~~ ", var_q, "*q", w
      )
      if (include_z) {
        model_string <- paste0(
          model_string, "\nr", w, " ~~ ", var_r, "*r", w
        )
      }
    }
  }

  # ==================== RESIDUAL COVARIANCES ====================
  model_string <- paste0(model_string, "\n\n# Residual Covariances")
  if (constrain_residual_covariances) {
    for (w in 1:waves) {
      model_string <- paste0(
        model_string, "\np", w, " ~~ ", cov_pq, "*q", w
      )
      if (include_z) {
        model_string <- paste0(
          model_string, "\np", w, " ~~ ", cov_pr, "*r", w,
          "\nq", w, " ~~ ", cov_qr, "*r", w
        )
      }
    }
  } else {
    # When unconstrained, each wave gets its own covariance (freely estimated)
    for (w in 1:waves) {
      model_string <- paste0(
        model_string, "\np", w, " ~~ q", w
      )
      if (include_z) {
        model_string <- paste0(
          model_string, "\np", w, " ~~ r", w,
          "\nq", w, " ~~ r", w
        )
      }
    }
  }

  # ==================== FIX OBSERVED VARIABLE RESIDUALS ====================
  model_string <- paste0(model_string, "\n\n# Fix Observed Variable Residuals to Zero")
  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\n", time_varying_x[w], " ~~ 0*", time_varying_x[w],
      "\n", time_varying_y[w], " ~~ 0*", time_varying_y[w]
    )
    if (include_z) {
      model_string <- paste0(
        model_string, "\n", time_varying_z[w], " ~~ 0*", time_varying_z[w]
      )
    }
  }

  # ==================== RANDOM INTERCEPT VARIANCES AND COVARIANCES ====================
  model_string <- paste0(model_string, "\n\n# Random Intercept Variances and Covariances")
  model_string <- paste0(
    model_string,
    "\nBX ~~ ", var_BX, "*BX",
    "\nBY ~~ ", var_BY, "*BY",
    "\nBX ~~ ", cov_BXBY, "*BY"
  )

  if (include_z) {
    model_string <- paste0(
      model_string,
      "\nBZ ~~ ", var_BZ, "*BZ",
      "\nBX ~~ ", cov_BXBZ, "*BZ",
      "\nBY ~~ ", cov_BYBZ, "*BZ"
    )
  }

  # ==================== RANDOM INTERCEPT CONSTRAINTS ====================
  # Constrain covariances between random intercepts and first wave within-person factors to 0
  model_string <- paste0(model_string, "\n\n# Random Intercept Constraints (no correlation with first wave)")
  model_string <- paste0(
    model_string,
    "\np1 ~~ 0*BX",
    "\np1 ~~ 0*BY",
    "\nq1 ~~ 0*BX",
    "\nq1 ~~ 0*BY"
  )

  if (include_z) {
    model_string <- paste0(
      model_string,
      "\nr1 ~~ 0*BX",
      "\nr1 ~~ 0*BY",
      "\nr1 ~~ 0*BZ",
      "\np1 ~~ 0*BZ",
      "\nq1 ~~ 0*BZ"
    )
  }

  # Store parameters for reference
  parameters <- list(
    waves = waves,
    sample.nobs = sample.nobs,
    include_z = include_z,
    beta_x = beta_x,
    beta_y = beta_y,
    beta_z = beta_z,
    omega_xy = omega_xy,
    omega_yx = omega_yx,
    omega_zx = omega_zx,
    omega_zy = omega_zy,
    omega_xz = omega_xz,
    omega_yz = omega_yz,
    var_p = var_p,
    var_q = var_q,
    var_r = var_r,
    cov_pq = cov_pq,
    cov_pr = cov_pr,
    cov_qr = cov_qr,
    var_BX = var_BX,
    var_BY = var_BY,
    var_BZ = var_BZ,
    cov_BXBY = cov_BXBY,
    cov_BXBZ = cov_BXBZ,
    cov_BYBZ = cov_BYBZ,
    mean_BX = mean_BX,
    mean_BY = mean_BY,
    mean_BZ = mean_BZ,
    constrain_beta = constrain_beta,
    constrain_omega = constrain_omega,
    constrain_residual_variances = constrain_residual_variances,
    constrain_residual_covariances = constrain_residual_covariances,
    estimate_means = estimate_means
  )

  # Generate data
  tryCatch({
    dat <- lavaan::simulateData(
      model = model_string,
      sample.nobs = sample.nobs,
      int.ov.free = TRUE,
      ...
    )

    # Validation message
    message(paste0("✅ Successfully simulated RICLPM data with ", waves, " waves",
                   if (include_z) " (including Z variable)" else "",
                   " and ", sample.nobs, " observations."))

    return(list(
      model = model_string,
      data = dat,
      parameters = parameters
    ))

  }, error = function(e) {
    stop(paste0("❌ Error in data simulation: ", e$message,
                "\nCheck your parameter values and model specification."))
  })
}
