#' @title estimateRICLPM_MI
#' @description Generates lavaan syntax for a multi-indicator Random Intercept
#'   Cross-Lagged Panel Model (RI-CLPM).
#'
#' @param waves Integer. Number of waves (>= 2).
#' @param n_indicators_x Integer. Indicators for X per wave (>= 1).
#' @param n_indicators_y Integer. Indicators for Y per wave (>= 1).
#' @param invariance Character. One of \code{"configural"}, \code{"metric"},
#'   \code{"scalar"}, \code{"strict"}. Default \code{"scalar"}.
#' @param correlated_residuals Logical. Correlate same-item residuals across
#'   waves. Default \code{TRUE}.
#' @param constrain_beta,constrain_omega,constrain_residual_variances,constrain_residual_covariances Logical.
#'   Equality constraints across waves.
#' @param start_values Logical.
#'
#' @return Character string of lavaan syntax.
#'
#' @details
#' Second-order structure: wave-specific latent factor \code{eta_x_w} aggregates
#' indicators \code{x{w}_1 ... x{w}_{J}}; it is then perfectly decomposed into a
#' between-person trait \code{I_x} and a within-person deviation \code{p_w}
#' (analogously for Y). All AR/CL dynamics operate on \code{p_w}, \code{q_w}.
#' Same-item residuals are correlated across waves (per Mulder & Hamaker, 2021).
#'
#' @references
#' Hamaker, E. L., Kuiper, R. M., & Grasman, R. P. P. P. (2015). A critique of
#'   the cross-lagged panel model. \emph{Psychological Methods}, 20, 102-116.
#'
#' Mulder, J. D., & Hamaker, E. L. (2021). Three extensions of the random
#'   intercept cross-lagged panel model. \emph{SEM}, 28(4), 638-648.
#'
#' Usami, S., Murayama, K., & Hamaker, E. L. (2019). \emph{Psychological
#'   Methods}, 24(5), 637-657.
#'
#' @examples
#' \dontrun{
#' syntax <- estimateRICLPM_MI(waves = 4, n_indicators_x = 3, n_indicators_y = 3)
#' fit <- lavaan::lavaan(syntax, data = my_data, meanstructure = TRUE)
#' }
#'
#' @export
estimateRICLPM_MI <- function(waves = 5,
                              n_indicators_x = 3,
                              n_indicators_y = 3,
                              invariance = c("scalar", "metric", "configural", "strict"),
                              correlated_residuals = TRUE,
                              constrain_beta = TRUE,
                              constrain_omega = TRUE,
                              constrain_residual_variances = TRUE,
                              constrain_residual_covariances = TRUE,
                              start_values = FALSE) {

  if (!is.numeric(waves) || waves < 2 || waves != as.integer(waves)) {
    stop("Error: 'waves' must be an integer >= 2.")
  }
  if (!is.numeric(n_indicators_x) || n_indicators_x < 1) {
    stop("Error: 'n_indicators_x' must be a positive integer.")
  }
  if (!is.numeric(n_indicators_y) || n_indicators_y < 1) {
    stop("Error: 'n_indicators_y' must be a positive integer.")
  }
  invariance <- match.arg(invariance)

  m <- ""

  # ==================== WAVE-SPECIFIC LATENT FACTORS ====================
  m <- paste0(m, "# Wave-specific latent factors (measurement model)\n")
  for (w in 1:waves) {
    lhs <- paste0("eta_x", w, " =~ 1*x", w, "_1")
    if (n_indicators_x >= 2) {
      for (j in 2:n_indicators_x) {
        lab <- if (invariance %in% c("metric", "scalar", "strict")) {
          paste0("lam_x", j)
        } else {
          paste0("lam_x", j, "_w", w)
        }
        lhs <- paste0(lhs, " + ", lab, "*x", w, "_", j)
      }
    }
    m <- paste0(m, lhs, "\n")

    lhs <- paste0("eta_y", w, " =~ 1*y", w, "_1")
    if (n_indicators_y >= 2) {
      for (j in 2:n_indicators_y) {
        lab <- if (invariance %in% c("metric", "scalar", "strict")) {
          paste0("lam_y", j)
        } else {
          paste0("lam_y", j, "_w", w)
        }
        lhs <- paste0(lhs, " + ", lab, "*y", w, "_", j)
      }
    }
    m <- paste0(m, lhs, "\n")
  }

  # ==================== INDICATOR INTERCEPTS ====================
  m <- paste0(m, "\n# Indicator intercepts (marker fixed to 0)\n")
  for (w in 1:waves) {
    m <- paste0(m, "x", w, "_1 ~ 0*1\n")
    if (n_indicators_x >= 2) {
      for (j in 2:n_indicators_x) {
        lab <- if (invariance %in% c("scalar", "strict")) {
          paste0("tau_x", j)
        } else {
          paste0("tau_x", j, "_w", w)
        }
        m <- paste0(m, "x", w, "_", j, " ~ ", lab, "*1\n")
      }
    }
    m <- paste0(m, "y", w, "_1 ~ 0*1\n")
    if (n_indicators_y >= 2) {
      for (j in 2:n_indicators_y) {
        lab <- if (invariance %in% c("scalar", "strict")) {
          paste0("tau_y", j)
        } else {
          paste0("tau_y", j, "_w", w)
        }
        m <- paste0(m, "y", w, "_", j, " ~ ", lab, "*1\n")
      }
    }
  }

  # ==================== INDICATOR RESIDUAL VARIANCES ====================
  m <- paste0(m, "\n# Indicator residual variances\n")
  for (w in 1:waves) {
    for (j in 1:n_indicators_x) {
      lab <- if (invariance == "strict") paste0("u_var_x", j) else paste0("u_var_x", j, "_w", w)
      m <- paste0(m, "x", w, "_", j, " ~~ ", lab, "*x", w, "_", j, "\n")
    }
    for (j in 1:n_indicators_y) {
      lab <- if (invariance == "strict") paste0("u_var_y", j) else paste0("u_var_y", j, "_w", w)
      m <- paste0(m, "y", w, "_", j, " ~~ ", lab, "*y", w, "_", j, "\n")
    }
  }

  # ==================== CORRELATED INDICATOR RESIDUALS ACROSS WAVES ====================
  if (correlated_residuals && waves >= 2) {
    m <- paste0(m, "\n# Same-item residual covariances across waves\n")
    for (j in 1:n_indicators_x) {
      for (v in 1:(waves - 1)) {
        for (w in (v + 1):waves) {
          m <- paste0(m, "x", v, "_", j, " ~~ x", w, "_", j, "\n")
        }
      }
    }
    for (j in 1:n_indicators_y) {
      for (v in 1:(waves - 1)) {
        for (w in (v + 1):waves) {
          m <- paste0(m, "y", v, "_", j, " ~~ y", w, "_", j, "\n")
        }
      }
    }
  }

  # ==================== WAVE-FACTOR RESIDUAL VARIANCES (ZERO) ====================
  # Decomposition: eta = I + p exactly; no residual variance at wave-factor level.
  m <- paste0(m, "\n# Wave-factor residual variances fixed to 0 (full decomposition)\n")
  for (w in 1:waves) {
    m <- paste0(m, "eta_x", w, " ~~ 0*eta_x", w, "\n")
    m <- paste0(m, "eta_y", w, " ~~ 0*eta_y", w, "\n")
  }

  # Wave-factor intercepts fixed to 0; mean carried by trait factor
  m <- paste0(m, "\n# Wave-factor intercepts fixed to 0\n")
  for (w in 1:waves) {
    m <- paste0(m, "eta_x", w, " ~ 0*1\n")
    m <- paste0(m, "eta_y", w, " ~ 0*1\n")
  }

  # ==================== TRAIT FACTORS ====================
  m <- paste0(m, "\n# Trait (random-intercept) factors\n")
  lhs <- paste0("I_x =~ 1*eta_x1")
  for (w in 2:waves) lhs <- paste0(lhs, " + 1*eta_x", w)
  m <- paste0(m, lhs, "\n")
  lhs <- paste0("I_y =~ 1*eta_y1")
  for (w in 2:waves) lhs <- paste0(lhs, " + 1*eta_y", w)
  m <- paste0(m, lhs, "\n")

  # ==================== WITHIN-PERSON FACTORS ====================
  m <- paste0(m, "\n# Within-person factors\n")
  for (w in 1:waves) {
    m <- paste0(m, "p", w, " =~ 1*eta_x", w, "\n")
    m <- paste0(m, "q", w, " =~ 1*eta_y", w, "\n")
  }

  # ==================== MEAN STRUCTURE ====================
  m <- paste0(m, "\n# Mean structure\n")
  m <- paste0(m, "I_x ~ 1\nI_y ~ 1\n")
  for (w in 1:waves) {
    m <- paste0(m, "p", w, " ~ 0*1\nq", w, " ~ 0*1\n")
  }

  # ==================== AR + CL ON WITHIN-PERSON FACTORS ====================
  m <- paste0(m, "\n# AR and CL on within-person factors\n")
  s_ar <- if (start_values) "start(0.3)*" else ""
  s_cl <- if (start_values) "start(0.1)*" else ""
  for (w in 2:waves) {
    ar_x <- if (constrain_beta) "ar_x" else paste0("ar_x", w)
    ar_y <- if (constrain_beta) "ar_y" else paste0("ar_y", w)
    cl_xy <- if (constrain_omega) "cl_xy" else paste0("cl_xy", w)
    cl_yx <- if (constrain_omega) "cl_yx" else paste0("cl_yx", w)
    m <- paste0(m, "p", w, " ~ ", s_ar, ar_x, "*p", w - 1,
                " + ", s_cl, cl_yx, "*q", w - 1, "\n")
    m <- paste0(m, "q", w, " ~ ", s_ar, ar_y, "*q", w - 1,
                " + ", s_cl, cl_xy, "*p", w - 1, "\n")
  }

  # ==================== DYNAMIC RESIDUAL VARIANCES (WITHIN) ====================
  m <- paste0(m, "\n# Dynamic residual variances on within-person factors\n")
  if (constrain_residual_variances) {
    m <- paste0(m, "p1 ~~ d_var_x1*p1\nq1 ~~ d_var_y1*q1\n")
    for (w in 2:waves) {
      m <- paste0(m, "p", w, " ~~ d_var_x*p", w, "\n")
      m <- paste0(m, "q", w, " ~~ d_var_y*q", w, "\n")
    }
  } else {
    for (w in 1:waves) {
      m <- paste0(m, "p", w, " ~~ p", w, "\nq", w, " ~~ q", w, "\n")
    }
  }

  # ==================== WITHIN-TIME COVARIANCES ====================
  m <- paste0(m, "\n# Within-time covariances\n")
  if (constrain_residual_covariances) {
    m <- paste0(m, "p1 ~~ d_cov_xy1*q1\n")
    for (w in 2:waves) {
      m <- paste0(m, "p", w, " ~~ d_cov_xy*q", w, "\n")
    }
  } else {
    for (w in 1:waves) {
      m <- paste0(m, "p", w, " ~~ q", w, "\n")
    }
  }

  # ==================== TRAIT FACTOR VARIANCES & COVARIANCE ====================
  m <- paste0(m, "\n# Trait factor variances and covariance\n")
  m <- paste0(m, "I_x ~~ I_x\nI_y ~~ I_y\nI_x ~~ I_y\n")

  return(m)
}
