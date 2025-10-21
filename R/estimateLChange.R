#' @title estimateLChange
#' @description Generates model syntax for Latent Change Score Models
#'
#' @param waves The number of waves (time points) in the model.
#' @param model_type Specify whether estimating "latent_change" (single variable) or "dual_change" (bivariate).
#' @param constrain_loadings Logical. If TRUE, constrains factor loadings to equality across waves. Default is TRUE.
#' @param constrain_ar Logical. If TRUE, constrains proportional effects (autoregressive) to equality. Default is TRUE.
#' @param constrain_coupling Logical. If TRUE, constrains coupling effects to equality. Default is TRUE for dual_change.
#' @param constrain_change_to_change Logical. If TRUE, constrains change-to-change effects to equality. Default is TRUE for dual_change.
#' @param constrain_phi Logical. If TRUE, constrains change autoregressions to equality. Default is TRUE for dual_change.
#'
#' @return A character string containing the Lavaan model syntax for the latent change score model.
#'
#' @details
#' This function generates the model syntax for a Latent Change Score Model with a specified
#' number of waves.
#'
#' For single variable models (model_type = "latent_change"):
#' - Latent true scores at each wave
#' - Latent change scores from wave 2 onwards
#' - A constant change factor
#' - Proportional effects (level affecting own change)
#'
#' For dual change models (model_type = "dual_change"):
#' - All components from single variable model for both X and Y
#' - Coupling effects (X level affecting Y change and vice versa)
#' - Change-to-change cross-lagged effects
#' - Change autoregression effects
#'
#' Key parameters estimated:
#' - ar_x, ar_y: Proportional effects (level -> own change)
#' - cl_x, cl_y: Coupling effects (other level -> change)
#' - change_x, change_y: Change-to-change effects
#' - phi_x, phi_y: Change autoregression
#' - alpha_g2: Constant change factor mean for X
#' - sigma2_lx1, sigma2_g2: Variances for initial level and constant change
#' - sigma2_u, sigma2_e: Measurement error variances
#'
#' @examples
#' # Single variable latent change model
#' model_syntax_single <- estimateLChange(waves = 5, model_type = "latent_change")
#' cat(model_syntax_single)
#'
#' # Dual change model
#' model_syntax_dual <- estimateLChange(waves = 5, model_type = "dual_change")
#' cat(model_syntax_dual)
#'
#' # Fit the model
#' library(lavaan)
#' data <- simLChange(waves = 5, model_type = "dual_change", sample.nobs = 1000)$data
#' fit <- lavaan(model_syntax_dual, data = data)
#' summary(fit)
#'
#' @export

estimateLChange <- function(waves = 10,
                            model_type = c('latent_change', 'dual_change'),
                            constrain_loadings = TRUE,
                            constrain_ar = TRUE,
                            constrain_coupling = TRUE,
                            constrain_change_to_change = TRUE,
                            constrain_phi = TRUE) {

  # Validate inputs
  if (!is.numeric(waves) || waves <= 0 || waves != as.integer(waves)) {
    stop("❌ Error: Parameter 'waves' must be a positive integer.")
  }

  if (waves < 3) {
    stop("❌ Error: Latent change models require at least 3 waves.")
  }

  model_type <- match.arg(model_type)

  model_string <- ""

  if (model_type == "latent_change") {
    # ==================== SINGLE VARIABLE MODEL ====================

    # Latent true scores (loadings = 1)
    for (w in 1:waves) {
      model_string <- paste0(model_string, "    cf", w, " =~ 1*x", w, "\n")
    }

    # Set the latent true score means (initial free, others = 0)
    model_string <- paste0(model_string, "    cf1 ~ 1\n")
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf", w, " ~ 0*1\n")
    }

    # Latent true score variances (initial free, others = 0)
    model_string <- paste0(model_string, "    cf1 ~~ cf1\n")
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf", w, " ~~ 0*cf", w, "\n")
    }

    # Observed intercepts (fixed to 0)
    for (w in 1:waves) {
      model_string <- paste0(model_string, "    x", w, " ~ 0*1\n")
    }

    # Observed residual variances (constrained to equality if specified)
    if (constrain_loadings) {
      for (w in 1:waves) {
        model_string <- paste0(model_string, "    x", w, " ~~ sigma2_u*x", w, "\n")
      }
    } else {
      for (w in 1:waves) {
        model_string <- paste0(model_string, "    x", w, " ~~ x", w, "\n")
      }
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

    # Constant change factor mean
    model_string <- paste0(model_string, "    general ~ 1\n")

    # Constant change factor variance
    model_string <- paste0(model_string, "    general ~~ general\n")

    # Constant change factor covariance with the initial true score
    model_string <- paste0(model_string, "    general ~~ cf1\n")

    # Proportional effects (constrained equal if specified)
    if (constrain_ar) {
      for (w in 2:waves) {
        model_string <- paste0(model_string, "    ld", w, " ~ ar_x*cf", w - 1, "\n")
      }
    } else {
      for (w in 2:waves) {
        model_string <- paste0(model_string, "    ld", w, " ~ cf", w - 1, "\n")
      }
    }

  } else if (model_type == "dual_change") {
    # ==================== DUAL CHANGE MODEL ====================

    # -------------------- X VARIABLE --------------------

    # Define the latent true scores for X
    for (w in 1:waves) {
      model_string <- paste0(model_string, "    cf_x", w, " =~ 1*x", w, "\n")
    }

    # Latent true score means (initial free, others = 0)
    model_string <- paste0(model_string, "    cf_x1 ~ gamma_lx1*1\n")
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf_x", w, " ~ 0*1\n")
    }

    # Latent true score variances (initial free, others = 0)
    model_string <- paste0(model_string, "    cf_x1 ~~ start(15)*sigma2_lx1*cf_x1\n")
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    cf_x", w, " ~~ 0*cf_x", w, "\n")
    }

    # Observed intercepts (fixed to 0)
    for (w in 1:waves) {
      model_string <- paste0(model_string, "    x", w, " ~ 0*1\n")
    }

    # Observed residual variances (constrained to equality)
    if (constrain_loadings) {
      for (w in 1:waves) {
        model_string <- paste0(model_string, "    x", w, " ~~ sigma2_u*x", w, "\n")
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

    # Latent change score variances (constrained to 0)
    for (w in 2:waves) {
      model_string <- paste0(model_string, "    ld_x", w, " ~~ 0*ld_x", w, "\n")
    }

    # Constant change factor (loadings = 1)
    model_string <- paste0(model_string, "    general_x =~ 1*ld_x2\n")
    for (w in 3:waves) {
      model_string <- paste0(model_string, "          + 1*ld_x", w, "\n")
    }

    # Constant change factor mean
    model_string <- paste0(model_string, "    general_x ~ alpha_g2*1\n")

    # Constant change factor variance
    model_string <- paste0(model_string, "    general_x ~~ sigma2_g2*general_x\n")

    # Constant change factor covariance with the initial true score
    model_string <- paste0(model_string, "    general_x ~~ cf_x1\n")

    # Proportional effects for X (constrained equal)
    if (constrain_ar) {
      for (w in 2:waves) {
        model_string <- paste0(model_string, "    ld_x", w, " ~ start(-0.15)*ar_x*cf_x", w - 1, "\n")
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

    # Observed residual variances (constrained to equality)
    if (constrain_loadings) {
      for (w in 1:waves) {
        model_string <- paste0(model_string, "    y", w, " ~~ sigma2_e*y", w, "\n")
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
    model_string <- paste0(model_string, "    general_y ~ start(15)*1\n")

    # Constant change factor variance
    model_string <- paste0(model_string, "    general_y ~~ general_y\n")

    # Constant change factor covariance with the initial true score
    model_string <- paste0(model_string, "    general_y ~~ cf_y1\n")

    # Proportional effects for Y (constrained equal)
    if (constrain_ar) {
      for (w in 2:waves) {
        model_string <- paste0(model_string, "    ld_y", w, " ~ ar_y*cf_y", w - 1, "\n")
      }
    } else {
      for (w in 2:waves) {
        model_string <- paste0(model_string, "    ld_y", w, " ~ cf_y", w - 1, "\n")
      }
    }

    # -------------------- CHANGE AUTOREGRESSIONS --------------------

    # Y change autoregression
    if (constrain_phi) {
      for (w in 3:waves) {
        model_string <- paste0(model_string, "    ld_y", w, " ~ phi_y*ld_y", w - 1, "\n")
      }
    } else {
      for (w in 3:waves) {
        model_string <- paste0(model_string, "    ld_y", w, " ~ ld_y", w - 1, "\n")
      }
    }

    # X change autoregression
    if (constrain_phi) {
      for (w in 3:waves) {
        model_string <- paste0(model_string, "    ld_x", w, " ~ phi_x*ld_x", w - 1, "\n")
      }
    } else {
      for (w in 3:waves) {
        model_string <- paste0(model_string, "    ld_x", w, " ~ ld_x", w - 1, "\n")
      }
    }

    # -------------------- CROSS-LAGGED EFFECTS --------------------

    # Change-to-change cross-lagged effects
    if (constrain_change_to_change) {
      for (w in 3:waves) {
        model_string <- paste0(model_string, "    ld_x", w, " ~ change_x*ld_y", w - 1, "\n")
        model_string <- paste0(model_string, "    ld_y", w, " ~ change_y*ld_x", w - 1, "\n")
      }
    } else {
      for (w in 3:waves) {
        model_string <- paste0(model_string, "    ld_x", w, " ~ ld_y", w - 1, "\n")
        model_string <- paste0(model_string, "    ld_y", w, " ~ ld_x", w - 1, "\n")
      }
    }

    # Coupling parameters (level affecting other variable's change)
    if (constrain_coupling) {
      for (w in 2:waves) {
        model_string <- paste0(model_string, "    ld_x", w, " ~ start(-0.2)*cl_y*cf_y", w - 1, "\n")
        model_string <- paste0(model_string, "    ld_y", w, " ~ start(-0.2)*cl_x*cf_x", w - 1, "\n")
      }
    } else {
      for (w in 2:waves) {
        model_string <- paste0(model_string, "    ld_x", w, " ~ start(-0.2)*cf_y", w - 1, "\n")
        model_string <- paste0(model_string, "    ld_y", w, " ~ start(-0.2)*cf_x", w - 1, "\n")
      }
    }

    # -------------------- COVARIANCES --------------------

    # Covariance between constant change factors
    model_string <- paste0(model_string, "    general_x ~~ general_y\n")

    # Covariance between initial true scores
    model_string <- paste0(model_string, "    cf_x1 ~~ cf_y1\n")

    # Covariance between wave 1 scores and constant change scores
    model_string <- paste0(model_string, "    cf_x1 ~~ general_y\n")
    model_string <- paste0(model_string, "    cf_y1 ~~ general_x\n")
  }

  return(model_string)
}
