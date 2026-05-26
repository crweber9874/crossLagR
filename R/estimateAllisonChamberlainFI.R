#' @title estimateAllisonChamberlainFI
#' @description Generate lavaan model syntax for the Allison-Chamberlain fixed-effects
#'   SEM approach to cross-lagged panel models.
#'
#' @param waves Integer. Number of waves (time points). Must be >= 3.
#' @param constrain_coefficients Logical. If TRUE, constrains autoregressive and
#'   cross-lagged effects to equality across waves. Default is TRUE.
#' @param constrain_residual_variances Logical. If TRUE, constrains residual
#'   variances to equality across waves (after wave 1). Default is TRUE.
#' @param correlate_errors_with_future Logical. If TRUE, allows residuals of each
#'   equation to be correlated with future values of the opposite predetermined
#'   variable, following the Chamberlain correction. Default is TRUE.
#' @param correlate_alpha_with_exogenous Logical. If TRUE, allows the latent
#'   fixed-effect variables (I_x, I_y) to be correlated with all observed
#'   x and y values. Default is TRUE.
#'
#' @return A character string containing lavaan model syntax for the
#'   Allison-Chamberlain fixed-effects SEM model.
#'
#' @details
#' This function generates lavaan syntax for the Allison-Chamberlain fixed-effects
#' SEM following the unified framework of Usami, Murayama, & Hamaker (2019). This
#' model includes latent fixed effects (analogous to trait factors I) and relaxes
#' strict exogeneity to predetermination.
#'
#' Parameter labels follow the unified naming convention:
#' \itemize{
#'   \item \code{ar_x}, \code{ar_y}: Autoregressive effects. These operate on
#'     observed scores with individual heterogeneity controlled via latent fixed
#'     effects.
#'   \item \code{cl_xy}: Cross-lagged effect of X on Y.
#'   \item \code{cl_yx}: Cross-lagged effect of Y on X.
#'   \item \code{d_var_x}, \code{d_var_y}: Dynamic residual variances.
#'   \item \code{I_x}, \code{I_y}: Latent fixed effects (individual heterogeneity).
#'     These function as trait factors in the unified framework, with unit loadings
#'     on endogenous outcomes.
#' }
#'
#' @references
#' Allison, P. D. (2005). \emph{Fixed Effects Regression Methods for Longitudinal
#'   Data Using SAS}. SAS Institute.
#'
#' Chamberlain, G. (1982). Multivariate regression models for panel data.
#'   \emph{Journal of Econometrics}, 18(1), 5-46.
#'
#' Usami, S., Murayama, K., & Hamaker, E. L. (2019). A unified framework of
#'   longitudinal models to examine reciprocal relations. \emph{Psychological
#'   Methods}, 24(5), 637-657.
#'
#' @examples
#' \dontrun{
#' # Generate Allison-Chamberlain FI syntax with 5 waves
#' syntax <- estimateAllisonChamberlainFI(waves = 5)
#' cat(syntax)
#'
#' # Fit with lavaan
#' library(lavaan)
#' fit <- lavaan::sem(syntax, data = my_data)
#' summary(fit, fit.measures = TRUE, standardized = TRUE)
#' }
#'
#' @export
estimateAllisonChamberlainFI <- function(waves = 5,
                                         constrain_coefficients = TRUE,
                                         constrain_residual_variances = TRUE,
                                         correlate_errors_with_future = TRUE,
                                         correlate_alpha_with_exogenous = TRUE) {
    # Input validation
    if (!is.numeric(waves) || waves <= 0 || waves != as.integer(waves)) {
        stop("Error: Parameter 'waves' must be a positive integer.")
    }
    if (waves < 3) {
        stop("Error: Allison-Chamberlain FI requires at least 3 waves.")
    }

    logical_params <- list(
        constrain_coefficients = constrain_coefficients,
        constrain_residual_variances = constrain_residual_variances,
        correlate_errors_with_future = correlate_errors_with_future,
        correlate_alpha_with_exogenous = correlate_alpha_with_exogenous
    )
    for (param_name in names(logical_params)) {
        if (!is.logical(logical_params[[param_name]])) {
            stop(paste0("Error: Parameter '", param_name, "' must be logical (TRUE/FALSE)."))
        }
    }

    model_string <- ""

    # ==================== LATENT FIXED EFFECTS (I) ====================
    model_string <- paste0(model_string, "    # Latent fixed effects (trait factors I)\n")
    model_string <- paste0(model_string, "    I_x =~ 1*x2")
    for (w in 3:waves) {
        model_string <- paste0(model_string, " + 1*x", w)
    }
    model_string <- paste0(model_string, "\n")

    model_string <- paste0(model_string, "    I_y =~ 1*y2")
    for (w in 3:waves) {
        model_string <- paste0(model_string, " + 1*y", w)
    }
    model_string <- paste0(model_string, "\n")

    # ==================== REGRESSION EQUATIONS ====================
    model_string <- paste0(model_string, "\n    # Regression equations\n")
    if (constrain_coefficients) {
        for (w in 2:waves) {
            model_string <- paste0(
                model_string,
                "    x", w, " ~ ar_x*x", w - 1, " + cl_yx*y", w - 1, "\n",
                "    y", w, " ~ ar_y*y", w - 1, " + cl_xy*x", w - 1, "\n"
            )
        }
    } else {
        for (w in 2:waves) {
            model_string <- paste0(
                model_string,
                "    x", w, " ~ x", w - 1, " + y", w - 1, "\n",
                "    y", w, " ~ y", w - 1, " + x", w - 1, "\n"
            )
        }
    }

    # ==================== ERROR-FUTURE CORRELATIONS ====================
    if (correlate_errors_with_future) {
        model_string <- paste0(model_string, "\n    # Chamberlain predetermined-variable correction\n")
        model_string <- paste0(model_string, "    # (error-future correlations)\n")
        for (w in 2:(waves - 1)) {
            for (w_future in (w + 1):waves) {
                model_string <- paste0(
                    model_string,
                    "    x", w, " ~~ y", w_future, "\n"
                )
                model_string <- paste0(
                    model_string,
                    "    y", w, " ~~ x", w_future, "\n"
                )
            }
        }
    }

    # ==================== RESIDUAL VARIANCES ====================
    model_string <- paste0(model_string, "\n    # Dynamic residual variances\n")
    if (constrain_residual_variances) {
        for (w in 2:waves) {
            model_string <- paste0(
                model_string,
                "    x", w, " ~~ d_var_x*x", w, "\n",
                "    y", w, " ~~ d_var_y*y", w, "\n"
            )
        }
    } else {
        for (w in 2:waves) {
            model_string <- paste0(
                model_string,
                "    x", w, " ~~ x", w, "\n",
                "    y", w, " ~~ y", w, "\n"
            )
        }
    }

    # ==================== I CORRELATIONS ====================
    if (correlate_alpha_with_exogenous) {
        model_string <- paste0(model_string, "\n    # I correlations with exogenous variables\n")

        model_string <- paste0(model_string, "    I_x ~~ x1\n")
        model_string <- paste0(model_string, "    I_x ~~ y1\n")
        model_string <- paste0(model_string, "    I_y ~~ x1\n")
        model_string <- paste0(model_string, "    I_y ~~ y1\n")

        for (w in 2:waves) {
            model_string <- paste0(model_string, "    I_x ~~ x", w, "\n")
            model_string <- paste0(model_string, "    I_x ~~ y", w, "\n")
            model_string <- paste0(model_string, "    I_y ~~ x", w, "\n")
            model_string <- paste0(model_string, "    I_y ~~ y", w, "\n")
        }
    }

    # ==================== I VARIANCES AND COVARIANCE ====================
    model_string <- paste0(model_string, "\n    # Latent fixed-effect (I) variances and covariance\n")
    model_string <- paste0(model_string, "    I_x ~~ I_x\n")
    model_string <- paste0(model_string, "    I_y ~~ I_y\n")
    model_string <- paste0(model_string, "    I_x ~~ I_y\n")

    return(model_string)
}
