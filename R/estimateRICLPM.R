#' @title estimateRICLPM
#' @description Generates model syntax for Random Intercept Cross-Lagged Panel Model (RICLPM) with hard-coded variable names and comprehensive parameter control options
#'
#' @param waves The number of waves (time points) in the model.
#' @param include_z Logical. If TRUE, includes a third variable Z in the model. Default is FALSE.
#' @param time_invariant_vars A vector of variable names for time-invariant variables.
#' @param constrain_beta Logical. If TRUE, constrains autoregressive effects to equality across waves. Default is TRUE.
#' @param constrain_omega Logical. If TRUE, constrains cross-lagged effects to equality across waves. Default is TRUE.
#' @param constrain_residual_variances Logical. If TRUE, constrains residual variances to equality across waves. Default is TRUE.
#' @param constrain_residual_covariances Logical. If TRUE, constrains residual covariances to equality across waves. Default is TRUE.
#' @param estimate_means Logical. If TRUE, estimates means for latent variables. Default is FALSE.
#' @param start_values Logical. If TRUE, provides starting values for key parameters. Default is FALSE.
#' @param label_autoregressive Character vector for autoregressive parameter labels. Default is c("beta_y", "beta_x", "beta_z").
#' @param label_crosslagged Character vector for cross-lagged parameter labels. Default is c("omega_xy", "omega_yx", "omega_zx", "omega_zy", "omega_xz", "omega_yz").
#'
#' @return A character string containing the Lavaan model syntax for RICLPM.
#'
#' @details
#' This function generates the model syntax for a Random Intercept Cross-Lagged Panel Model (RICLPM)
#' with a specified number of waves. The RICLPM separates stable between-person differences
#' (random intercepts) from within-person dynamics (autoregressive and cross-lagged effects).
#'
#' Variable naming convention:
#' - X variables: x1, x2, x3, ..., x[waves]
#' - Y variables: y1, y2, y3, ..., y[waves]
#' - Z variables (if included): z1, z2, z3, ..., z[waves]
#'
#' Your data must contain these exact variable names for the specified number of waves.
#'
#' Parameter labels follow consistent naming:
#' - Autoregressive effects: beta_x (X→X), beta_y (Y→Y), beta_z (Z→Z)
#' - Cross-lagged effects: omega_xy (Y→X), omega_yx (X→Y), omega_zx (Z→X),
#'   omega_zy (Z→Y), omega_xz (X→Z), omega_yz (Y→Z)
#'
#' Key model components:
#' - Random intercepts (BX, BY, BZ) that capture stable between-person differences
#' - Wave-specific latent variables (p, q, r) that capture within-person deviations
#' - Autoregressive effects within variables across time
#' - Cross-lagged effects between variables across time
#' - Within-time covariances between variables
#'
#' The model assumes that observed variables are indicators of latent within-person
#' deviations plus stable individual differences (random intercepts).
#'
#' Constraint options allow for:
#' - Equality constraints on autoregressive effects across waves
#' - Equality constraints on cross-lagged effects across waves
#' - Equality constraints on residual variances across waves
#' - Equality constraints on residual covariances across waves
#'
#' When constraints are applied, parameters are constrained to equality across all waves
#' from wave 2 onwards (since wave 1 has no predictors).
#'
#' @examples
#' \dontrun{
#' # Basic RICLPM with 5 waves (assumes data has x1-x5, y1-y5)
#' model_syntax <- estimateRICLPM(waves = 5)
#'
#' # RICLPM with unconstrained cross-lagged effects
#' model_syntax <- estimateRICLPM(waves = 5, constrain_omega = FALSE)
#'
#' # RICLPM with third variable (assumes data has x1-x5, y1-y5, z1-z5)
#' model_syntax <- estimateRICLPM(waves = 5, include_z = TRUE)
#'
#' # Fit the model
#' library(lavaan)
#' fitted_model <- lavaan(model_syntax, data = your_data)
#' summary(fitted_model, fit.measures = TRUE, standardized = TRUE)
#' }
#'
#' @export
estimateRICLPM <- function(waves,
                           include_z = FALSE,
                           time_invariant_vars = NULL,
                           constrain_beta = TRUE,
                           constrain_omega = TRUE,
                           constrain_residual_variances = TRUE,
                           constrain_residual_covariances = TRUE,
                           estimate_means = FALSE,
                           start_values = FALSE,
                           label_autoregressive = c("beta_y", "beta_x", "beta_z"),
                           label_crosslagged = c("omega_xy", "omega_yx", "omega_zx", "omega_zy", "omega_xz", "omega_yz")) {

  # Input validation
  if (!is.numeric(waves) || waves <= 0 || waves != as.integer(waves)) {
    stop("❌ Error: 'waves' must be a positive integer.")
  }

  if (waves < 2) {
    stop("❌ Error: RICLPM requires at least 2 waves. Cannot estimate with fewer than 2 waves.")
  }

  if (!is.logical(include_z)) {
    stop("❌ Error: 'include_z' must be logical (TRUE/FALSE).")
  }

  if (!is.null(time_invariant_vars) && !is.character(time_invariant_vars)) {
    stop("❌ Error: 'time_invariant_vars' must be a character vector of variable names (e.g., c('gender', 'education')).")
  }

  # Generate hard-coded variable names
  time_varying_x <- paste0("x", 1:waves)
  time_varying_y <- paste0("y", 1:waves)
  time_varying_z <- if (include_z) paste0("z", 1:waves) else NULL

  # Validate logical parameters
  logical_params <- list(
    constrain_beta = constrain_beta,
    constrain_omega = constrain_omega,
    constrain_residual_variances = constrain_residual_variances,
    constrain_residual_covariances = constrain_residual_covariances,
    estimate_means = estimate_means,
    start_values = start_values
  )

  for (param_name in names(logical_params)) {
    if (!is.logical(logical_params[[param_name]])) {
      stop(paste0("❌ Error: Parameter '", param_name, "' must be logical (TRUE/FALSE)."))
    }
  }

  # Initialize model string
  model_string <- ""

  # ==================== RANDOM INTERCEPTS ====================
  # Random intercept for X
  model_string <- paste0(model_string, "\n# Random Intercepts\nBX =~ 1*", time_varying_x[1])
  for (w in 2:waves) {
    model_string <- paste0(model_string, " + 1*", time_varying_x[w])
  }

  # Random intercept for Y
  model_string <- paste0(model_string, "\nBY =~ 1*", time_varying_y[1])
  for (w in 2:waves) {
    model_string <- paste0(model_string, " + 1*", time_varying_y[w])
  }

  # Random intercept for Z (if specified)
  if (!is.null(time_varying_z)) {
    model_string <- paste0(model_string, "\nBZ =~ 1*", time_varying_z[1])
    for (w in 2:waves) {
      model_string <- paste0(model_string, " + 1*", time_varying_z[w])
    }
  }

  model_string <- paste0(model_string, "\n")

  # ==================== WITHIN-PERSON LATENT VARIABLES ====================
  model_string <- paste0(model_string, "\n# Within-Person Latent Variables")
  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\np", w, " =~ 1*", time_varying_x[w],
      "\nq", w, " =~ 1*", time_varying_y[w]
    )
    # Add Z variable if specified
    if (!is.null(time_varying_z)) {
      model_string <- paste0(
        model_string, "\nr", w, " =~ 1*", time_varying_z[w]
      )
    }
  }

  # ==================== MEANS ====================
  if (estimate_means) {
    model_string <- paste0(model_string, "\n\n# Latent Variable Means")
    for (w in 1:waves) {
      model_string <- paste0(model_string, "\np", w, " ~ 1\n")
      model_string <- paste0(model_string, "q", w, " ~ 1\n")
      if (!is.null(time_varying_z)) {
        model_string <- paste0(model_string, "r", w, " ~ 1\n")
      }
    }
  } else {
    model_string <- paste0(model_string, "\n\n# Latent Variable Means (Fixed to 0)")
    for (w in 1:waves) {
      model_string <- paste0(model_string, "\np", w, " ~ 0*1\n")
      model_string <- paste0(model_string, "q", w, " ~ 0*1\n")
      if (!is.null(time_varying_z)) {
        model_string <- paste0(model_string, "r", w, " ~ 0*1\n")
      }
    }
  }

  # ==================== AUTOREGRESSIVE AND CROSS-LAGGED PATHS ====================
  model_string <- paste0(model_string, "\n# Autoregressive and Cross-lagged Effects")

  # Starting values if requested
  start_ar <- if (start_values) "start(0.3)*" else ""
  start_cl <- if (start_values) "start(0.1)*" else ""

  for (w in 2:waves) {
    if (constrain_beta && constrain_omega) {
      # Both autoregressive and cross-lagged effects constrained
      model_string <- paste0(
        model_string, "\np", w, " ~ ", start_ar, label_autoregressive[1], "*p", w - 1,
        " + ", start_cl, label_crosslagged[1], "*q", w - 1
      )
      model_string <- paste0(
        model_string, "\nq", w, " ~ ", start_ar, label_autoregressive[2], "*q", w - 1,
        " + ", start_cl, label_crosslagged[2], "*p", w - 1
      )

      if (!is.null(time_varying_z)) {
        model_string <- paste0(
          model_string, "\nr", w, " ~ ", start_ar, label_autoregressive[3], "*r", w - 1,
          " + ", start_cl, label_crosslagged[3], "*p", w - 1,
          " + ", start_cl, label_crosslagged[4], "*q", w - 1
        )
        model_string <- paste0(
          model_string, "\np", w, " ~ ", start_cl, label_crosslagged[5], "*r", w - 1
        )
        model_string <- paste0(
          model_string, "\nq", w, " ~ ", start_cl, label_crosslagged[6], "*r", w - 1
        )
      }
    } else if (constrain_beta && !constrain_omega) {
      # Only autoregressive effects constrained
      model_string <- paste0(
        model_string, "\np", w, " ~ ", start_ar, label_autoregressive[1], "*p", w - 1,
        " + ", start_cl, "omega_xy", w, "*q", w - 1
      )
      model_string <- paste0(
        model_string, "\nq", w, " ~ ", start_ar, label_autoregressive[2], "*q", w - 1,
        " + ", start_cl, "omega_yx", w, "*p", w - 1
      )

      if (!is.null(time_varying_z)) {
        model_string <- paste0(
          model_string, "\nr", w, " ~ ", start_ar, label_autoregressive[3], "*r", w - 1,
          " + ", start_cl, "omega_zx", w, "*p", w - 1,
          " + ", start_cl, "omega_zy", w, "*q", w - 1
        )
        model_string <- paste0(
          model_string, "\np", w, " ~ ", start_cl, "omega_xz", w, "*r", w - 1
        )
        model_string <- paste0(
          model_string, "\nq", w, " ~ ", start_cl, "omega_yz", w, "*r", w - 1
        )
      }
    } else if (!constrain_beta && constrain_omega) {
      # Only cross-lagged effects constrained
      model_string <- paste0(
        model_string, "\np", w, " ~ ", start_ar, "beta_y", w, "*p", w - 1,
        " + ", start_cl, label_crosslagged[1], "*q", w - 1
      )
      model_string <- paste0(
        model_string, "\nq", w, " ~ ", start_ar, "beta_x", w, "*q", w - 1,
        " + ", start_cl, label_crosslagged[2], "*p", w - 1
      )

      if (!is.null(time_varying_z)) {
        model_string <- paste0(
          model_string, "\nr", w, " ~ ", start_ar, "beta_z", w, "*r", w - 1,
          " + ", start_cl, label_crosslagged[3], "*p", w - 1,
          " + ", start_cl, label_crosslagged[4], "*q", w - 1
        )
        model_string <- paste0(
          model_string, "\np", w, " ~ ", start_cl, label_crosslagged[5], "*r", w - 1
        )
        model_string <- paste0(
          model_string, "\nq", w, " ~ ", start_cl, label_crosslagged[6], "*r", w - 1
        )
      }
    } else {
      # No constraints - all effects free
      model_string <- paste0(
        model_string, "\np", w, " ~ ", start_ar, "beta_y", w, "*p", w - 1,
        " + ", start_cl, "omega_xy", w, "*q", w - 1
      )
      model_string <- paste0(
        model_string, "\nq", w, " ~ ", start_ar, "beta_x", w, "*q", w - 1,
        " + ", start_cl, "omega_yx", w, "*p", w - 1
      )

      if (!is.null(time_varying_z)) {
        model_string <- paste0(
          model_string, "\nr", w, " ~ ", start_ar, "beta_z", w, "*r", w - 1,
          " + ", start_cl, "omega_zx", w, "*p", w - 1,
          " + ", start_cl, "omega_zy", w, "*q", w - 1
        )
        model_string <- paste0(
          model_string, "\np", w, " ~ ", start_cl, "omega_xz", w, "*r", w - 1
        )
        model_string <- paste0(
          model_string, "\nq", w, " ~ ", start_cl, "omega_yz", w, "*r", w - 1
        )
      }
    }
  }

  # ==================== RESIDUAL VARIANCES ====================
  model_string <- paste0(model_string, "\n\n# Residual Variances")
  if (constrain_residual_variances) {
    for (w in 1:waves) {
      model_string <- paste0(
        model_string, "\np", w, " ~~ var_p*p", w,
        "\nq", w, " ~~ var_q*q", w
      )
      if (!is.null(time_varying_z)) {
        model_string <- paste0(
          model_string, "\nr", w, " ~~ var_r*r", w
        )
      }
    }
  } else {
    for (w in 1:waves) {
      model_string <- paste0(
        model_string, "\np", w, " ~~ p", w,
        "\nq", w, " ~~ q", w
      )
      if (!is.null(time_varying_z)) {
        model_string <- paste0(
          model_string, "\nr", w, " ~~ r", w
        )
      }
    }
  }

  # ==================== RESIDUAL COVARIANCES ====================
  model_string <- paste0(model_string, "\n\n# Residual Covariances")
  if (constrain_residual_covariances) {
    for (w in 1:waves) {
      model_string <- paste0(
        model_string, "\np", w, " ~~ cov_pq*q", w
      )
      if (!is.null(time_varying_z)) {
        model_string <- paste0(
          model_string, "\np", w, " ~~ cov_pr*r", w,
          "\nq", w, " ~~ cov_qr*r", w
        )
      }
    }
  } else {
    for (w in 1:waves) {
      model_string <- paste0(
        model_string, "\np", w, " ~~ q", w
      )
      if (!is.null(time_varying_z)) {
        model_string <- paste0(
          model_string, "\np", w, " ~~ r", w,
          "\nq", w, " ~~ r", w
        )
      }
    }
  }

  # ==================== TIME-INVARIANT PREDICTORS ====================
  if (!is.null(time_invariant_vars)) {
    model_string <- paste0(model_string, "\n\n# Time-Invariant Predictors")
    for (var in time_invariant_vars) {
      model_string <- paste0(model_string, "\n", var, " ~~ ", var)
      for (w in 1:waves) {
        model_string <- paste0(
          model_string, "\np", w, " ~ ", var,
          "\nq", w, " ~ ", var
        )
        if (!is.null(time_varying_z)) {
          model_string <- paste0(
            model_string, "\nr", w, " ~ ", var
          )
        }
      }
    }
  }

  # ==================== RANDOM INTERCEPT VARIANCES AND COVARIANCES ====================
  model_string <- paste0(model_string, "\n\n# Random Intercept Variances and Covariances")
  model_string <- paste0(
    model_string,
    "\nBX ~~ BX",
    "\nBY ~~ BY",
    "\nBX ~~ BY"
  )

  if (!is.null(time_varying_z)) {
    model_string <- paste0(
      model_string,
      "\nBZ ~~ BZ",
      "\nBX ~~ BZ",
      "\nBY ~~ BZ"
    )
  }

  # ==================== FIX OBSERVED VARIABLE RESIDUALS ====================
  model_string <- paste0(model_string, "\n\n# Fix Observed Variable Residuals to Zero")
  for (w in 1:waves) {
    model_string <- paste0(
      model_string, "\n", time_varying_x[w], " ~~ 0*", time_varying_x[w],
      "\n", time_varying_y[w], " ~~ 0*", time_varying_y[w]
    )
    if (!is.null(time_varying_z)) {
      model_string <- paste0(
        model_string, "\n", time_varying_z[w], " ~~ 0*", time_varying_z[w]
      )
    }
  }

  return(model_string)
}
