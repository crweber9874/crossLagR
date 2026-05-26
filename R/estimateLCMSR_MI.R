#' @title estimateLCMSR_MI
#' @description Generates lavaan syntax for a multi-indicator bivariate Latent
#'   Curve Model with Structured Residuals (LCM-SR / GCLM).
#'
#' @param waves Integer. Number of waves (>= 3).
#' @param n_indicators_x,n_indicators_y Integer. Indicators per wave.
#' @param invariance Character. One of \code{"configural"}, \code{"metric"},
#'   \code{"scalar"}, \code{"strict"}. Default \code{"scalar"}.
#' @param correlated_residuals Logical. Default \code{TRUE}.
#' @param time_scores Numeric vector length \code{waves}. Default
#'   \code{0:(waves-1)}.
#' @param constrain_beta,constrain_omega,constrain_residual_variances,constrain_residual_covariances Logical.
#' @param estimate_quadratic Logical. Add quadratic growth factor (requires
#'   \code{waves >= 4}). Default \code{FALSE}.
#'
#' @return Character string of lavaan syntax.
#'
#' @details
#' Curve-of-factors layout: indicators -> wave latent \code{eta_x_w}; growth
#' factors \code{I_x}, \code{S_x} load on \code{eta_x_w} (direct effects only),
#' and \code{p_w} captures the structured residual. AR/CL operate on
#' \code{p_w}, \code{q_w}. Same-item residuals are correlated across waves.
#'
#' @references
#' Curran, P. J., Howard, A. L., Bainter, S. A., Lane, S. T., & McGinley, J. S.
#'   (2014). \emph{Annual Review of Psychology}, 65, 637-660.
#'
#' Usami, S., Murayama, K., & Hamaker, E. L. (2019). \emph{Psychological
#'   Methods}, 24(5), 637-657.
#'
#' @export
estimateLCMSR_MI <- function(waves = 5,
                             n_indicators_x = 3,
                             n_indicators_y = 3,
                             invariance = c("scalar", "metric", "configural", "strict"),
                             correlated_residuals = TRUE,
                             time_scores = NULL,
                             constrain_beta = TRUE,
                             constrain_omega = TRUE,
                             constrain_residual_variances = TRUE,
                             constrain_residual_covariances = TRUE,
                             estimate_quadratic = FALSE) {

  if (!is.numeric(waves) || waves < 3) stop("Error: 'waves' must be >= 3.")
  if (n_indicators_x < 1 || n_indicators_y < 1) stop("Error: indicator counts must be positive.")
  if (estimate_quadratic && waves < 4) stop("Error: quadratic growth requires waves >= 4.")
  invariance <- match.arg(invariance)

  if (is.null(time_scores)) {
    time_scores <- 0:(waves - 1)
  } else if (!is.numeric(time_scores) || length(time_scores) != waves) {
    stop("Error: 'time_scores' length must equal 'waves'.")
  }
  quad_scores <- if (estimate_quadratic) time_scores^2 else NULL

  m <- ""

  # ==================== MEASUREMENT MODEL ====================
  m <- paste0(m, "# Wave-specific latent factors (measurement model)\n")
  for (w in 1:waves) {
    lhs <- paste0("eta_x", w, " =~ 1*x", w, "_1")
    if (n_indicators_x >= 2) for (j in 2:n_indicators_x) {
      lab <- if (invariance %in% c("metric", "scalar", "strict")) paste0("lam_x", j)
      else paste0("lam_x", j, "_w", w)
      lhs <- paste0(lhs, " + ", lab, "*x", w, "_", j)
    }
    m <- paste0(m, lhs, "\n")

    lhs <- paste0("eta_y", w, " =~ 1*y", w, "_1")
    if (n_indicators_y >= 2) for (j in 2:n_indicators_y) {
      lab <- if (invariance %in% c("metric", "scalar", "strict")) paste0("lam_y", j)
      else paste0("lam_y", j, "_w", w)
      lhs <- paste0(lhs, " + ", lab, "*y", w, "_", j)
    }
    m <- paste0(m, lhs, "\n")
  }

  # Indicator intercepts
  m <- paste0(m, "\n# Indicator intercepts\n")
  for (w in 1:waves) {
    m <- paste0(m, "x", w, "_1 ~ 0*1\n")
    if (n_indicators_x >= 2) for (j in 2:n_indicators_x) {
      lab <- if (invariance %in% c("scalar", "strict")) paste0("tau_x", j) else paste0("tau_x", j, "_w", w)
      m <- paste0(m, "x", w, "_", j, " ~ ", lab, "*1\n")
    }
    m <- paste0(m, "y", w, "_1 ~ 0*1\n")
    if (n_indicators_y >= 2) for (j in 2:n_indicators_y) {
      lab <- if (invariance %in% c("scalar", "strict")) paste0("tau_y", j) else paste0("tau_y", j, "_w", w)
      m <- paste0(m, "y", w, "_", j, " ~ ", lab, "*1\n")
    }
  }

  # Indicator residual variances
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

  # Correlated residuals across waves
  if (correlated_residuals && waves >= 2) {
    m <- paste0(m, "\n# Same-item residual covariances across waves\n")
    for (j in 1:n_indicators_x) for (v in 1:(waves - 1)) for (w in (v + 1):waves) {
      m <- paste0(m, "x", v, "_", j, " ~~ x", w, "_", j, "\n")
    }
    for (j in 1:n_indicators_y) for (v in 1:(waves - 1)) for (w in (v + 1):waves) {
      m <- paste0(m, "y", v, "_", j, " ~~ y", w, "_", j, "\n")
    }
  }

  # Wave-factor residual variances and intercepts = 0 (full decomposition)
  m <- paste0(m, "\n# Wave-factor residuals fixed to 0 (full decomposition)\n")
  for (w in 1:waves) {
    m <- paste0(m, "eta_x", w, " ~~ 0*eta_x", w, "\n")
    m <- paste0(m, "eta_y", w, " ~~ 0*eta_y", w, "\n")
    m <- paste0(m, "eta_x", w, " ~ 0*1\n")
    m <- paste0(m, "eta_y", w, " ~ 0*1\n")
  }

  # ==================== GROWTH FACTORS (DIRECT EFFECTS ON eta) ====================
  m <- paste0(m, "\n# Growth factors (load on wave latents; direct effects only)\n")
  lhs <- paste0("I_x =~ 1*eta_x1")
  for (w in 2:waves) lhs <- paste0(lhs, " + 1*eta_x", w)
  m <- paste0(m, lhs, "\n")
  lhs <- paste0("S_x =~ ", time_scores[1], "*eta_x1")
  for (w in 2:waves) lhs <- paste0(lhs, " + ", time_scores[w], "*eta_x", w)
  m <- paste0(m, lhs, "\n")
  if (estimate_quadratic) {
    lhs <- paste0("Q_x =~ ", quad_scores[1], "*eta_x1")
    for (w in 2:waves) lhs <- paste0(lhs, " + ", quad_scores[w], "*eta_x", w)
    m <- paste0(m, lhs, "\n")
  }
  lhs <- paste0("I_y =~ 1*eta_y1")
  for (w in 2:waves) lhs <- paste0(lhs, " + 1*eta_y", w)
  m <- paste0(m, lhs, "\n")
  lhs <- paste0("S_y =~ ", time_scores[1], "*eta_y1")
  for (w in 2:waves) lhs <- paste0(lhs, " + ", time_scores[w], "*eta_y", w)
  m <- paste0(m, lhs, "\n")
  if (estimate_quadratic) {
    lhs <- paste0("Q_y =~ ", quad_scores[1], "*eta_y1")
    for (w in 2:waves) lhs <- paste0(lhs, " + ", quad_scores[w], "*eta_y", w)
    m <- paste0(m, lhs, "\n")
  }

  # ==================== STRUCTURED RESIDUAL FACTORS ====================
  m <- paste0(m, "\n# Structured residual (within-person) factors\n")
  for (w in 1:waves) {
    m <- paste0(m, "p", w, " =~ 1*eta_x", w, "\n")
    m <- paste0(m, "q", w, " =~ 1*eta_y", w, "\n")
  }

  # Latent means on p, q fixed to 0
  m <- paste0(m, "\n# Within-person factor means fixed to 0\n")
  for (w in 1:waves) {
    m <- paste0(m, "p", w, " ~ 0*1\nq", w, " ~ 0*1\n")
  }

  # ==================== AR + CL ====================
  m <- paste0(m, "\n# AR and CL on structured residuals\n")
  for (w in 2:waves) {
    ar_x <- if (constrain_beta) "ar_x" else paste0("ar_x", w)
    ar_y <- if (constrain_beta) "ar_y" else paste0("ar_y", w)
    cl_xy <- if (constrain_omega) "cl_xy" else paste0("cl_xy", w)
    cl_yx <- if (constrain_omega) "cl_yx" else paste0("cl_yx", w)
    m <- paste0(m, "p", w, " ~ ", ar_x, "*p", w - 1, " + ", cl_yx, "*q", w - 1, "\n")
    m <- paste0(m, "q", w, " ~ ", ar_y, "*q", w - 1, " + ", cl_xy, "*p", w - 1, "\n")
  }

  # Wave-1 residual (free), endogenous residuals
  m <- paste0(m, "\n# Dynamic residual variances\n")
  m <- paste0(m, "p1 ~~ p1\nq1 ~~ q1\np1 ~~ q1\n")
  if (constrain_residual_variances) {
    for (w in 2:waves) {
      m <- paste0(m, "p", w, " ~~ d_var_x*p", w, "\n")
      m <- paste0(m, "q", w, " ~~ d_var_y*q", w, "\n")
    }
  } else {
    for (w in 2:waves) {
      m <- paste0(m, "p", w, " ~~ p", w, "\nq", w, " ~~ q", w, "\n")
    }
  }
  if (constrain_residual_covariances) {
    for (w in 2:waves) m <- paste0(m, "p", w, " ~~ d_cov_xy*q", w, "\n")
  } else {
    for (w in 2:waves) m <- paste0(m, "p", w, " ~~ q", w, "\n")
  }

  # ==================== GROWTH FACTOR MEANS, VARS, COVS ====================
  m <- paste0(m, "\n# Growth factor means\n")
  m <- paste0(m, "I_x ~ 1\nS_x ~ 1\nI_y ~ 1\nS_y ~ 1\n")
  if (estimate_quadratic) m <- paste0(m, "Q_x ~ 1\nQ_y ~ 1\n")

  m <- paste0(m, "\n# Growth factor variances\n")
  m <- paste0(m, "I_x ~~ I_x\nS_x ~~ S_x\nI_y ~~ I_y\nS_y ~~ S_y\n")
  if (estimate_quadratic) m <- paste0(m, "Q_x ~~ Q_x\nQ_y ~~ Q_y\n")

  m <- paste0(m, "\n# Growth factor covariances (within variable)\n")
  m <- paste0(m, "I_x ~~ S_x\nI_y ~~ S_y\n")
  if (estimate_quadratic) {
    m <- paste0(m, "I_x ~~ Q_x\nS_x ~~ Q_x\nI_y ~~ Q_y\nS_y ~~ Q_y\n")
  }

  m <- paste0(m, "\n# Growth factor covariances (between variables)\n")
  m <- paste0(m, "I_x ~~ I_y\nS_x ~~ S_y\nI_x ~~ S_y\nI_y ~~ S_x\n")
  if (estimate_quadratic) {
    m <- paste0(m, "Q_x ~~ Q_y\nI_x ~~ Q_y\nI_y ~~ Q_x\nS_x ~~ Q_y\nS_y ~~ Q_x\n")
  }

  return(m)
}
