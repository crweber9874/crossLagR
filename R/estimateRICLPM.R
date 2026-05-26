#' @title estimateRICLPM
#' @description Generates lavaan model syntax for the Random Intercept Cross-Lagged
#'   Panel Model (RI-CLPM).
#'
#' @param waves The number of waves (time points) in the model.
#' @param include_z Logical. If TRUE, includes a third variable Z in the model. Default is FALSE.
#' @param time_invariant_vars A vector of variable names for time-invariant variables.
#' @param constrain_beta Logical. If TRUE, constrains autoregressive effects to equality across waves. Default is TRUE.
#' @param constrain_omega Logical. If TRUE, constrains cross-lagged effects to equality across waves. Default is TRUE.
#' @param constrain_residual_variances Logical. If TRUE, constrains residual variances to equality across waves. Default is TRUE.
#' @param constrain_residual_covariances Logical. If TRUE, constrains residual covariances to equality across waves. Default is TRUE.
#' @param start_values Logical. If TRUE, provides starting values for key parameters. Default is FALSE.
#' @param label_autoregressive Character vector for autoregressive parameter labels.
#'   Default is \code{c("ar_x", "ar_y", "ar_z")}.
#' @param label_crosslagged Character vector for cross-lagged parameter labels.
#'   Default is \code{c("cl_yx", "cl_xy", "cl_xz", "cl_yz", "cl_zx", "cl_zy")}.
#'
#' @return A character string containing the lavaan model syntax for the RI-CLPM.
#'
#' @details
#' This function generates lavaan syntax for a Random Intercept Cross-Lagged Panel
#' Model following the unified framework of Usami, Murayama, & Hamaker (2019). The
#' RI-CLPM uses both the decomposition equation (separating stable trait factors
#' from within-person deviations) and the dynamic equation (lagged regression on
#' the within-person deviations).
#'
#' Parameter labels follow the unified naming convention:
#' \itemize{
#'   \item \code{ar_x}, \code{ar_y}, \code{ar_z}: Autoregressive effects. In the
#'     RI-CLPM these represent \emph{within-person carry-over}, not rank-order
#'     stability (which is captured by the trait factor I).
#'   \item \code{cl_xy}, \code{cl_yx}: Cross-lagged effects operating at the
#'     within-person level. \code{cl_xy} = effect of X on Y; \code{cl_yx} = effect
#'     of Y on X.
#'   \item \code{d_var_x}, \code{d_var_y}, \code{d_var_z}: Dynamic residual
#'     (innovation) variances for within-person deviations.
#'   \item \code{d_cov_xy}, \code{d_cov_xz}, \code{d_cov_yz}: Dynamic residual
#'     covariances (within-time, within-person associations).
#'   \item \code{I_x}, \code{I_y}, \code{I_z}: Stable trait factors (random
#'     intercepts) representing time-invariant between-person differences. These
#'     are unique to models that include decomposition of trait and state (no
#'     direct analogue in the CLPM).
#' }
#'
#' Within-person latent deviations are labeled \code{p} (for X), \code{q} (for Y),
#' and \code{r} (for Z). These structural latent variables have no direct analogue
#' across all models and retain model-specific names.
#'
#' @references
#' Usami, S., Murayama, K., & Hamaker, E. L. (2019). A unified framework of
#'   longitudinal models to examine reciprocal relations. \emph{Psychological
#'   Methods}, 24(5), 637-657.
#'
#' Hamaker, E. L., Kuiper, R. M., & Grasman, R. P. P. P. (2015). A critique
#'   of the cross-lagged panel model. \emph{Psychological Methods}, 20, 102-116.
#'
#' @examples
#' \dontrun{
#' # Basic RI-CLPM with 5 waves (assumes data has x1-x5, y1-y5)
#' model_syntax <- estimateRICLPM(waves = 5)
#'
#' # RI-CLPM with unconstrained cross-lagged effects
#' model_syntax <- estimateRICLPM(waves = 5, constrain_omega = FALSE)
#'
#' # RI-CLPM with third variable (assumes data has x1-x5, y1-y5, z1-z5)
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
                           start_values = FALSE,
                           label_autoregressive = c("ar_x", "ar_y", "ar_z"),
                           label_crosslagged = c("cl_yx", "cl_xy", "cl_xz", "cl_yz", "cl_zx", "cl_zy")) {

  # Input validation
  if (!is.numeric(waves) || waves <= 0 || waves != as.integer(waves)) {
    stop("Error: 'waves' must be a positive integer.")
  }

  if (waves < 2) {
    stop("Error: RICLPM requires at least 2 waves. Cannot estimate with fewer than 2 waves.")
  }

  if (!is.logical(include_z)) {
    stop("Error: 'include_z' must be logical (TRUE/FALSE).")
  }

  if (!is.null(time_invariant_vars) && !is.character(time_invariant_vars)) {
    stop("Error: 'time_invariant_vars' must be a character vector of variable names (e.g., c('gender', 'education')).")
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
    start_values = start_values
  )

  for (param_name in names(logical_params)) {
    if (!is.logical(logical_params[[param_name]])) {
      stop(paste0("Error: Parameter '", param_name, "' must be logical (TRUE/FALSE)."))
    }
  }

  # Initialize model string
  model_string <- ""

  # ==================== RANDOM INTERCEPTS (Trait Factors I) ====================
  # Random intercept for X
  model_string <- paste0(model_string, "\n# Random Intercepts (Trait Factors I)\nI_x =~ 1*", time_varying_x[1])
  for (w in 2:waves) {
    model_string <- paste0(model_string, " + 1*", time_varying_x[w])
  }

  # Random intercept for Y
  model_string <- paste0(model_string, "\nI_y =~ 1*", time_varying_y[1])
  for (w in 2:waves) {
    model_string <- paste0(model_string, " + 1*", time_varying_y[w])
  }

  # Random intercept for Z (if specified)
  if (!is.null(time_varying_z)) {
    model_string <- paste0(model_string, "\nI_z =~ 1*", time_varying_z[1])
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

  # ==================== MEAN STRUCTURE ====================
  # Canonical RI-CLPM (Mulder & Hamaker 2021):
  #   - Within-person factors p_w, q_w have zero means (deviations by definition).
  #   - Trait factors I_x, I_y have free means (interpretable as the average
  #     between-person level).
  #   - Observed intercepts are fixed to zero so the trait-factor mean is
  #     identified rather than absorbed into wave-specific intercepts.
  # Requires lavaan() to be called with meanstructure = TRUE (or use sem()/cfa()).
  model_string <- paste0(model_string, "\n\n# Within-person factor means (fixed to 0)")
  for (w in 1:waves) {
    model_string <- paste0(model_string, "\np", w, " ~ 0*1")
    model_string <- paste0(model_string, "\nq", w, " ~ 0*1")
    if (!is.null(time_varying_z)) {
      model_string <- paste0(model_string, "\nr", w, " ~ 0*1")
    }
  }
  model_string <- paste0(model_string, "\n")

  model_string <- paste0(model_string, "\n# Trait factor means (free)")
  model_string <- paste0(model_string, "\nI_x ~ 1")
  model_string <- paste0(model_string, "\nI_y ~ 1")
  if (!is.null(time_varying_z)) {
    model_string <- paste0(model_string, "\nI_z ~ 1")
  }
  model_string <- paste0(model_string, "\n")

  model_string <- paste0(model_string, "\n# Observed intercepts (fixed to 0)")
  for (w in 1:waves) {
    model_string <- paste0(model_string, "\n", time_varying_x[w], " ~ 0*1")
    model_string <- paste0(model_string, "\n", time_varying_y[w], " ~ 0*1")
    if (!is.null(time_varying_z)) {
      model_string <- paste0(model_string, "\n", time_varying_z[w], " ~ 0*1")
    }
  }
  model_string <- paste0(model_string, "\n")

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
        " + ", start_cl, "cl_yx", w, "*q", w - 1
      )
      model_string <- paste0(
        model_string, "\nq", w, " ~ ", start_ar, label_autoregressive[2], "*q", w - 1,
        " + ", start_cl, "cl_xy", w, "*p", w - 1
      )

      if (!is.null(time_varying_z)) {
        model_string <- paste0(
          model_string, "\nr", w, " ~ ", start_ar, label_autoregressive[3], "*r", w - 1,
          " + ", start_cl, "cl_xz", w, "*p", w - 1,
          " + ", start_cl, "cl_yz", w, "*q", w - 1
        )
        model_string <- paste0(
          model_string, "\np", w, " ~ ", start_cl, "cl_zx", w, "*r", w - 1
        )
        model_string <- paste0(
          model_string, "\nq", w, " ~ ", start_cl, "cl_zy", w, "*r", w - 1
        )
      }
    } else if (!constrain_beta && constrain_omega) {
      # Only cross-lagged effects constrained
      model_string <- paste0(
        model_string, "\np", w, " ~ ", start_ar, "ar_x", w, "*p", w - 1,
        " + ", start_cl, label_crosslagged[1], "*q", w - 1
      )
      model_string <- paste0(
        model_string, "\nq", w, " ~ ", start_ar, "ar_y", w, "*q", w - 1,
        " + ", start_cl, label_crosslagged[2], "*p", w - 1
      )

      if (!is.null(time_varying_z)) {
        model_string <- paste0(
          model_string, "\nr", w, " ~ ", start_ar, "ar_z", w, "*r", w - 1,
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
        model_string, "\np", w, " ~ ", start_ar, "ar_x", w, "*p", w - 1,
        " + ", start_cl, "cl_yx", w, "*q", w - 1
      )
      model_string <- paste0(
        model_string, "\nq", w, " ~ ", start_ar, "ar_y", w, "*q", w - 1,
        " + ", start_cl, "cl_xy", w, "*p", w - 1
      )

      if (!is.null(time_varying_z)) {
        model_string <- paste0(
          model_string, "\nr", w, " ~ ", start_ar, "ar_z", w, "*r", w - 1,
          " + ", start_cl, "cl_xz", w, "*p", w - 1,
          " + ", start_cl, "cl_yz", w, "*q", w - 1
        )
        model_string <- paste0(
          model_string, "\np", w, " ~ ", start_cl, "cl_zx", w, "*r", w - 1
        )
        model_string <- paste0(
          model_string, "\nq", w, " ~ ", start_cl, "cl_zy", w, "*r", w - 1
        )
      }
    }
  }

  # ==================== RESIDUAL VARIANCES ====================
  model_string <- paste0(model_string, "\n\n# Dynamic Residual Variances")
  if (constrain_residual_variances) {
    for (w in 1:waves) {
      model_string <- paste0(
        model_string, "\np", w, " ~~ d_var_x*p", w,
        "\nq", w, " ~~ d_var_y*q", w
      )
      if (!is.null(time_varying_z)) {
        model_string <- paste0(
          model_string, "\nr", w, " ~~ d_var_z*r", w
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
  model_string <- paste0(model_string, "\n\n# Dynamic Residual Covariances")
  if (constrain_residual_covariances) {
    for (w in 1:waves) {
      model_string <- paste0(
        model_string, "\np", w, " ~~ d_cov_xy*q", w
      )
      if (!is.null(time_varying_z)) {
        model_string <- paste0(
          model_string, "\np", w, " ~~ d_cov_xz*r", w,
          "\nq", w, " ~~ d_cov_yz*r", w
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

  # ==================== TRAIT FACTOR VARIANCES AND COVARIANCES ====================
  model_string <- paste0(model_string, "\n\n# Trait Factor (I) Variances and Covariances")
  model_string <- paste0(
    model_string,
    "\nI_x ~~ I_x",
    "\nI_y ~~ I_y",
    "\nI_x ~~ I_y"
  )

  if (!is.null(time_varying_z)) {
    model_string <- paste0(
      model_string,
      "\nI_z ~~ I_z",
      "\nI_x ~~ I_z",
      "\nI_y ~~ I_z"
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
