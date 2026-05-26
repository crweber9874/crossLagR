#' @title estimateBollen_and_Brand_MI
#' @description Generates lavaan syntax for a multi-indicator Bollen & Brand
#'   (2010) dynamic panel model.
#'
#' @param waves Integer. Number of waves (>= 3).
#' @param n_indicators_x,n_indicators_y Integer. Indicators per wave.
#' @param invariance Character. One of \code{"configural"}, \code{"metric"},
#'   \code{"scalar"}, \code{"strict"}. Default \code{"scalar"}.
#' @param correlated_residuals Logical. Default \code{TRUE}.
#' @param x_effect Character. \code{"none"}, \code{"concurrent"}, or
#'   \code{"lagged"}. Default \code{"lagged"}.
#' @param x_autoregression Logical. If TRUE, x latents have AR paths and a
#'   latent fixed effect I_x. Default \code{TRUE}.
#' @param y_effect_on_x Logical. If TRUE, reciprocal y_{t-1} -> x_t.
#'   Default \code{TRUE}.
#' @param constrain_coefficients Logical. Equality constraints across waves.
#' @param constrain_residual_variances Logical.
#'
#' @return Character string of lavaan syntax.
#'
#' @details
#' Latent wave factors \code{p_w} (X) and \code{q_w} (Y) with multi-indicators
#' replace the manifest x_w/y_w of the standard Bollen-Brand model. Latent
#' fixed effects \code{I_y}, \code{I_x} load with 1 on the endogenous wave
#' latents (\code{q2..qT} and \code{p2..pT} respectively). Same-item residuals
#' are correlated across waves.
#'
#' @references
#' Bollen, K. A., & Brand, J. E. (2010). \emph{Social Forces}, 89(1), 1-34.
#'
#' Dishop, C. R., & DeShon, R. P. (2018). \emph{Psychological Methods}, 23(4),
#'   1089-1112.
#'
#' Usami, S., Murayama, K., & Hamaker, E. L. (2019). \emph{Psychological
#'   Methods}, 24(5), 637-657.
#'
#' @export
estimateBollen_and_Brand_MI <- function(waves = 5,
                                        n_indicators_x = 3,
                                        n_indicators_y = 3,
                                        invariance = c("scalar", "metric", "configural", "strict"),
                                        correlated_residuals = TRUE,
                                        x_effect = "lagged",
                                        x_autoregression = TRUE,
                                        y_effect_on_x = TRUE,
                                        constrain_coefficients = FALSE,
                                        constrain_residual_variances = TRUE) {

  if (!is.numeric(waves) || waves < 3) stop("Error: 'waves' must be >= 3.")
  if (n_indicators_x < 1 || n_indicators_y < 1) stop("Error: indicator counts must be positive.")
  if (!x_effect %in% c("none", "concurrent", "lagged")) {
    stop("Error: 'x_effect' must be 'none', 'concurrent', or 'lagged'.")
  }
  if (y_effect_on_x && !x_autoregression) stop("Error: 'y_effect_on_x' requires 'x_autoregression = TRUE'.")
  if (y_effect_on_x && x_effect == "none") stop("Error: 'y_effect_on_x' requires non-'none' x_effect.")
  if (x_autoregression && x_effect == "none") stop("Error: 'x_autoregression' requires non-'none' x_effect.")
  invariance <- match.arg(invariance)

  has_x <- x_effect != "none"
  if (has_x) {
    x_pred_of_y <- if (x_effect == "concurrent") 2:waves else 1:(waves - 1)
  }
  if (has_x && x_autoregression) {
    x_endo_waves <- if (y_effect_on_x) 2:waves else 2:max(x_pred_of_y)
  }

  m <- ""

  # ==================== MEASUREMENT MODEL ====================
  m <- paste0(m, "# Wave-specific latent factors\n")
  for (w in 1:waves) {
    lhs <- paste0("q", w, " =~ 1*y", w, "_1")
    if (n_indicators_y >= 2) for (j in 2:n_indicators_y) {
      lab <- if (invariance %in% c("metric", "scalar", "strict")) paste0("lam_y", j) else paste0("lam_y", j, "_w", w)
      lhs <- paste0(lhs, " + ", lab, "*y", w, "_", j)
    }
    m <- paste0(m, lhs, "\n")
    if (has_x) {
      lhs <- paste0("p", w, " =~ 1*x", w, "_1")
      if (n_indicators_x >= 2) for (j in 2:n_indicators_x) {
        lab <- if (invariance %in% c("metric", "scalar", "strict")) paste0("lam_x", j) else paste0("lam_x", j, "_w", w)
        lhs <- paste0(lhs, " + ", lab, "*x", w, "_", j)
      }
      m <- paste0(m, lhs, "\n")
    }
  }

  # Indicator intercepts
  m <- paste0(m, "\n# Indicator intercepts\n")
  for (w in 1:waves) {
    m <- paste0(m, "y", w, "_1 ~ 0*1\n")
    if (n_indicators_y >= 2) for (j in 2:n_indicators_y) {
      lab <- if (invariance %in% c("scalar", "strict")) paste0("tau_y", j) else paste0("tau_y", j, "_w", w)
      m <- paste0(m, "y", w, "_", j, " ~ ", lab, "*1\n")
    }
    if (has_x) {
      m <- paste0(m, "x", w, "_1 ~ 0*1\n")
      if (n_indicators_x >= 2) for (j in 2:n_indicators_x) {
        lab <- if (invariance %in% c("scalar", "strict")) paste0("tau_x", j) else paste0("tau_x", j, "_w", w)
        m <- paste0(m, "x", w, "_", j, " ~ ", lab, "*1\n")
      }
    }
  }

  # Indicator residual variances
  m <- paste0(m, "\n# Indicator residual variances\n")
  for (w in 1:waves) {
    for (j in 1:n_indicators_y) {
      lab <- if (invariance == "strict") paste0("u_var_y", j) else paste0("u_var_y", j, "_w", w)
      m <- paste0(m, "y", w, "_", j, " ~~ ", lab, "*y", w, "_", j, "\n")
    }
    if (has_x) for (j in 1:n_indicators_x) {
      lab <- if (invariance == "strict") paste0("u_var_x", j) else paste0("u_var_x", j, "_w", w)
      m <- paste0(m, "x", w, "_", j, " ~~ ", lab, "*x", w, "_", j, "\n")
    }
  }

  # Correlated residuals across waves
  if (correlated_residuals && waves >= 2) {
    m <- paste0(m, "\n# Same-item residual covariances across waves\n")
    for (j in 1:n_indicators_y) for (v in 1:(waves - 1)) for (w in (v + 1):waves) {
      m <- paste0(m, "y", v, "_", j, " ~~ y", w, "_", j, "\n")
    }
    if (has_x) for (j in 1:n_indicators_x) for (v in 1:(waves - 1)) for (w in (v + 1):waves) {
      m <- paste0(m, "x", v, "_", j, " ~~ x", w, "_", j, "\n")
    }
  }

  # ==================== LATENT FIXED EFFECTS (I) ====================
  m <- paste0(m, "\n# Latent fixed effects (accumulating)\n")
  lhs <- paste0("I_y =~ 1*q2")
  if (waves >= 3) for (w in 3:waves) lhs <- paste0(lhs, " + 1*q", w)
  m <- paste0(m, lhs, "\n")
  if (has_x && x_autoregression) {
    lhs <- paste0("I_x =~ ", paste0("1*p", x_endo_waves, collapse = " + "))
    m <- paste0(m, lhs, "\n")
  }

  # ==================== REGRESSION EQUATIONS ====================
  m <- paste0(m, "\n# Regression equations\n")
  for (w in 2:waves) {
    ar_y_lab <- if (constrain_coefficients) "ar_y" else paste0("ar_y", w)
    eq <- paste0("q", w, " ~ ", ar_y_lab, "*q", w - 1)
    if (has_x) {
      cl_xy_lab <- if (constrain_coefficients) "cl_xy" else paste0("cl_xy", w)
      eq <- if (x_effect == "concurrent") paste0(eq, " + ", cl_xy_lab, "*p", w)
            else paste0(eq, " + ", cl_xy_lab, "*p", w - 1)
    }
    m <- paste0(m, eq, "\n")
  }
  if (has_x && x_autoregression) {
    for (w in x_endo_waves) {
      ar_x_lab <- if (constrain_coefficients) "ar_x" else paste0("ar_x", w)
      eq <- paste0("p", w, " ~ ", ar_x_lab, "*p", w - 1)
      if (y_effect_on_x) {
        cl_yx_lab <- if (constrain_coefficients) "cl_yx" else paste0("cl_yx", w)
        eq <- paste0(eq, " + ", cl_yx_lab, "*q", w - 1)
      }
      m <- paste0(m, eq, "\n")
    }
  }

  # ==================== EXOGENOUS LATENT VARIANCES/COVARIANCES ====================
  exog_lat <- "q1"
  if (has_x) {
    if (x_autoregression) {
      exog_lat <- c(exog_lat, "p1")
    } else {
      exog_lat <- c(exog_lat, paste0("p", x_pred_of_y))
    }
  }
  m <- paste0(m, "\n# Exogenous wave variances\n")
  for (ev in exog_lat) m <- paste0(m, ev, " ~~ ", ev, "\n")
  if (length(exog_lat) > 1) {
    m <- paste0(m, "\n# Exogenous wave covariances\n")
    for (i in 1:(length(exog_lat) - 1)) for (j in (i + 1):length(exog_lat)) {
      m <- paste0(m, exog_lat[i], " ~~ ", exog_lat[j], "\n")
    }
  }
  m <- paste0(m, "\n# Exogenous wave means free\n")
  for (ev in exog_lat) m <- paste0(m, ev, " ~ 1\n")

  # ==================== I VARIANCES AND COVARIANCES ====================
  m <- paste0(m, "\n# Latent fixed effect variances and covariances\n")
  m <- paste0(m, "I_y ~~ I_y\n")
  i_vars <- "I_y"
  if (has_x && x_autoregression) {
    m <- paste0(m, "I_x ~~ I_x\nI_y ~~ I_x\n")
    i_vars <- c(i_vars, "I_x")
  }
  m <- paste0(m, "\n# Exogenous-I covariances\n")
  for (ev in exog_lat) for (iv in i_vars) m <- paste0(m, ev, " ~~ ", iv, "\n")

  # Endogenous latent means fixed to 0
  m <- paste0(m, "\n# Endogenous latent means fixed to 0\n")
  for (w in 2:waves) m <- paste0(m, "q", w, " ~ 0*1\n")
  if (has_x && x_autoregression) for (w in x_endo_waves) m <- paste0(m, "p", w, " ~ 0*1\n")

  # ==================== DYNAMIC RESIDUAL VARIANCES ====================
  m <- paste0(m, "\n# Dynamic residual variances\n")
  if (constrain_residual_variances) {
    for (w in 2:waves) m <- paste0(m, "q", w, " ~~ d_var_y*q", w, "\n")
    if (has_x && x_autoregression) for (w in x_endo_waves) m <- paste0(m, "p", w, " ~~ d_var_x*p", w, "\n")
  } else {
    for (w in 2:waves) m <- paste0(m, "q", w, " ~~ q", w, "\n")
    if (has_x && x_autoregression) for (w in x_endo_waves) m <- paste0(m, "p", w, " ~~ p", w, "\n")
  }

  return(m)
}
