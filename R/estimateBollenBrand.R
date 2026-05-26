#' @title estimateBollen_and_Brand
#' @description Generate lavaan model syntax for the Bollen and Brand (2010) dynamic
#'   panel model, as described in Dishop & DeShon (2018).
#'
#' @param waves Integer. Number of waves (time points). Must be >= 3.
#' @param x_effect Character. How x influences y: \code{"none"} (single variable
#'   autoregressive model), \code{"concurrent"} (x_t predicts y_t), or
#'   \code{"lagged"} (x_{t-1} predicts y_t). Default is \code{"lagged"}.
#' @param x_autoregression Logical. If TRUE, x has autoregressive paths and gets
#'   its own latent fixed effect (I_x). When FALSE, x is treated as fully
#'   exogenous. Default is FALSE.
#' @param y_effect_on_x Logical. If TRUE, includes reciprocal cross-lagged paths
#'   from y to x (y_{t-1} predicts x_t). Requires \code{x_autoregression = TRUE}.
#'   Default is FALSE.
#' @param constrain_coefficients Logical. If TRUE, constrains autoregressive and
#'   cross-lagged effects to equality across waves (stationarity). Default is TRUE.
#' @param constrain_residual_variances Logical. If TRUE, constrains residual
#'   variances to equality across waves. Default is TRUE.
#'
#' @return A character string containing lavaan model syntax.
#'
#' @details
#' This function generates lavaan syntax for the Bollen & Brand dynamic panel model
#' following the unified framework of Usami, Murayama, & Hamaker (2019). This model
#' is closely related to the ALT model, using accumulating factors with lagged
#' regression. It conditions on the first observation and treats it as exogenous.
#'
#' Parameter labels follow the unified naming convention:
#' \itemize{
#'   \item \code{ar_y}, \code{ar_x}: Autoregressive effects. In this model these
#'     operate on observed scores (no within/between decomposition), similar to the
#'     ALT model.
#'   \item \code{cl_xy}: Cross-lagged effect of X on Y.
#'   \item \code{cl_yx}: Cross-lagged effect of Y on X (reciprocal, if enabled).
#'   \item \code{d_var_y}, \code{d_var_x}: Dynamic residual variances.
#'   \item \code{I_y}, \code{I_x}: Latent fixed effects (individual heterogeneity),
#'     functioning as accumulating factors with unit loadings on endogenous outcomes.
#'     In the unified framework these are analogous to the A factor in the ALT/LCS
#'     models.
#' }
#'
#' @references
#' Bollen, K. A., & Brand, J. E. (2010). A general panel model with random and
#'   fixed effects: A structural equations approach. \emph{Social Forces}, 89(1), 1-34.
#'
#' Dishop, C. R., & DeShon, R. P. (2018). A tutorial on Bollen and Brand's
#'   approach to modeling dynamics. \emph{Psychological Methods}, 23(4), 1089-1112.
#'
#' Usami, S., Murayama, K., & Hamaker, E. L. (2019). A unified framework of
#'   longitudinal models to examine reciprocal relations. \emph{Psychological
#'   Methods}, 24(5), 637-657.
#'
#' @examples
#' \dontrun{
#' # Reciprocal cross-lag dynamic panel (Dishop & DeShon Figure 5)
#' syntax <- estimateBollen_and_Brand(waves = 5, x_effect = "lagged",
#'                                     x_autoregression = TRUE, y_effect_on_x = TRUE)
#' cat(syntax)
#'
#' # Fit with lavaan
#' library(lavaan)
#' fit <- lavaan::sem(syntax, data = my_wide_data)
#' summary(fit, fit.measures = TRUE, standardized = TRUE)
#' }
#'
#' @export
estimateBollen_and_Brand <- function(waves = 5,
                                      x_effect = "lagged",
                                      x_autoregression = TRUE,
                                      y_effect_on_x = TRUE,
                                      constrain_coefficients = FALSE,
                                      constrain_residual_variances = TRUE) {

    # ==================== INPUT VALIDATION ====================
    if (!is.numeric(waves) || waves <= 0 || waves != as.integer(waves)) {
        stop("Error: 'waves' must be a positive integer.")
    }
    if (waves < 3) {
        stop("Error: Dynamic panel model requires at least 3 waves.")
    }
    if (!x_effect %in% c("none", "concurrent", "lagged")) {
        stop("Error: 'x_effect' must be 'none', 'concurrent', or 'lagged'.")
    }
    if (y_effect_on_x && !x_autoregression) {
        stop("Error: 'y_effect_on_x' requires 'x_autoregression = TRUE'.")
    }
    if (y_effect_on_x && x_effect == "none") {
        stop("Error: 'y_effect_on_x' requires 'x_effect' to be 'concurrent' or 'lagged'.")
    }
    if (x_autoregression && x_effect == "none") {
        stop("Error: 'x_autoregression' requires 'x_effect' to be 'concurrent' or 'lagged'.")
    }

    has_x <- x_effect != "none"
    model <- ""

    # Determine which x waves serve as predictors of y
    if (has_x) {
        if (x_effect == "concurrent") {
            x_pred_of_y <- 2:waves
        } else {
            x_pred_of_y <- 1:(waves - 1)
        }
    }

    # Determine endogenous x waves (those predicted by AR)
    if (has_x && x_autoregression) {
        if (y_effect_on_x) {
            x_endo_waves <- 2:waves
        } else {
            x_endo_waves <- 2:max(x_pred_of_y)
        }
    }

    # Determine exogenous variables
    exog_vars <- "y1"
    if (has_x) {
        if (x_autoregression) {
            exog_vars <- c(exog_vars, "x1")
        } else {
            exog_vars <- c(exog_vars, paste0("x", x_pred_of_y))
        }
    }

    # ==================== LATENT FIXED EFFECTS (I) ====================
    model <- paste0(model, "# Latent fixed effects (individual heterogeneity, accumulating factors)\n")
    model <- paste0(model, "I_y =~ ", paste0("1*y", 2:waves, collapse = " + "), "\n")

    if (has_x && x_autoregression) {
        model <- paste0(model, "I_x =~ ", paste0("1*x", x_endo_waves, collapse = " + "), "\n")
    }

    # ==================== REGRESSION EQUATIONS ====================
    model <- paste0(model, "\n# Regression equations\n")

    for (w in 2:waves) {
        ar_y <- if (constrain_coefficients) "ar_y" else paste0("ar_y", w)

        y_eq <- paste0("y", w, " ~ ", ar_y, "*y", w - 1)

        # Cross-effect of x on y
        if (has_x) {
            cl_xy <- if (constrain_coefficients) "cl_xy" else paste0("cl_xy", w)
            if (x_effect == "concurrent") {
                y_eq <- paste0(y_eq, " + ", cl_xy, "*x", w)
            } else {
                y_eq <- paste0(y_eq, " + ", cl_xy, "*x", w - 1)
            }
        }
        model <- paste0(model, y_eq, "\n")
    }

    # X autoregressive equations
    if (has_x && x_autoregression) {
        for (w in x_endo_waves) {
            ar_x <- if (constrain_coefficients) "ar_x" else paste0("ar_x", w)
            x_eq <- paste0("x", w, " ~ ", ar_x, "*x", w - 1)

            # Reciprocal: y -> x
            if (y_effect_on_x) {
                cl_yx <- if (constrain_coefficients) "cl_yx" else paste0("cl_yx", w)
                x_eq <- paste0(x_eq, " + ", cl_yx, "*y", w - 1)
            }
            model <- paste0(model, x_eq, "\n")
        }
    }

    # ==================== EXOGENOUS VARIABLE VARIANCES ====================
    model <- paste0(model, "\n# Exogenous variable variances\n")
    for (ev in exog_vars) {
        model <- paste0(model, ev, " ~~ ", ev, "\n")
    }

    # ==================== EXOGENOUS VARIABLE COVARIANCES ====================
    if (length(exog_vars) > 1) {
        model <- paste0(model, "\n# Exogenous variable covariances\n")
        for (i in 1:(length(exog_vars) - 1)) {
            for (j in (i + 1):length(exog_vars)) {
                model <- paste0(model, exog_vars[i], " ~~ ", exog_vars[j], "\n")
            }
        }
    }

    # ==================== I VARIANCES AND COVARIANCES ====================
    model <- paste0(model, "\n# Latent fixed effect (I) variances and covariances\n")
    model <- paste0(model, "I_y ~~ I_y\n")

    if (has_x && x_autoregression) {
        model <- paste0(model, "I_x ~~ I_x\n")
        model <- paste0(model, "I_y ~~ I_x\n")
    }

    # ==================== EXOGENOUS-I COVARIANCES ====================
    model <- paste0(model, "\n# Exogenous-I covariances\n")
    i_vars <- "I_y"
    if (has_x && x_autoregression) {
        i_vars <- c(i_vars, "I_x")
    }

    for (ev in exog_vars) {
        for (iv in i_vars) {
            model <- paste0(model, ev, " ~~ ", iv, "\n")
        }
    }

    # ==================== RESIDUAL VARIANCES ====================
    model <- paste0(model, "\n# Dynamic residual variances\n")
    if (constrain_residual_variances) {
        for (w in 2:waves) {
            model <- paste0(model, "y", w, " ~~ d_var_y*y", w, "\n")
        }
        if (has_x && x_autoregression) {
            for (w in x_endo_waves) {
                model <- paste0(model, "x", w, " ~~ d_var_x*x", w, "\n")
            }
        }
    } else {
        for (w in 2:waves) {
            model <- paste0(model, "y", w, " ~~ y", w, "\n")
        }
        if (has_x && x_autoregression) {
            for (w in x_endo_waves) {
                model <- paste0(model, "x", w, " ~~ x", w, "\n")
            }
        }
    }

    return(model)
}
