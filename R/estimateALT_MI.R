#' @title estimateALT_MI
#' @description Generates lavaan syntax for a multi-indicator bivariate
#'   Autoregressive Latent Trajectory (ALT) model.
#'
#' @param waves Integer. Number of waves (>= 3).
#' @param n_indicators_x,n_indicators_y Integer. Indicators per wave.
#' @param invariance Character. One of \code{"configural"}, \code{"metric"},
#'   \code{"scalar"}, \code{"strict"}. Default \code{"scalar"}.
#' @param correlated_residuals Logical. Default \code{TRUE}.
#' @param time_scores Numeric vector length \code{waves}. Default
#'   \code{0:(waves-1)}.
#' @param constrain_beta,constrain_omega,constrain_residual_variances,constrain_residual_covariances Logical.
#'
#' @return Character string of lavaan syntax.
#'
#' @details
#' Wave-specific latent factors \code{p_w}, \code{q_w} are defined directly from
#' multiple indicators. Growth factors \code{I_x}, \code{S_x}, \code{I_y},
#' \code{S_y} are accumulating: they load on waves 2..T only, so they have
#' direct effects on endogenous waves AND indirect effects through the AR/CL
#' dynamics. Wave 1 is exogenous. Same-item residuals are correlated across
#' waves.
#'
#' @references
#' Bollen, K. A., & Curran, P. J. (2006). \emph{Latent curve models}. Wiley.
#'
#' Curran, P. J., & Bollen, K. A. (2001). The best of both worlds. In Collins
#'   & Sayer (Eds.), \emph{New methods for the analysis of change} (107-135). APA.
#'
#' Usami, S., Murayama, K., & Hamaker, E. L. (2019). \emph{Psychological
#'   Methods}, 24(5), 637-657.
#'
#' @export
estimateALT_MI <- function(waves = 5,
                           n_indicators_x = 3,
                           n_indicators_y = 3,
                           invariance = c("scalar", "metric", "configural", "strict"),
                           correlated_residuals = TRUE,
                           time_scores = NULL,
                           constrain_beta = TRUE,
                           constrain_omega = TRUE,
                           constrain_residual_variances = TRUE,
                           constrain_residual_covariances = TRUE) {

  if (!is.numeric(waves) || waves < 3) stop("Error: 'waves' must be >= 3.")
  if (n_indicators_x < 1 || n_indicators_y < 1) stop("Error: indicator counts must be positive.")
  invariance <- match.arg(invariance)

  if (is.null(time_scores)) {
    time_scores <- 0:(waves - 1)
  } else if (!is.numeric(time_scores) || length(time_scores) != waves) {
    stop("Error: 'time_scores' length must equal 'waves'.")
  }

  m <- ""

  # ==================== MEASUREMENT MODEL ====================
  m <- paste0(m, "# Wave-specific latent factors (multi-indicator)\n")
  for (w in 1:waves) {
    lhs <- paste0("p", w, " =~ 1*x", w, "_1")
    if (n_indicators_x >= 2) for (j in 2:n_indicators_x) {
      lab <- if (invariance %in% c("metric", "scalar", "strict")) paste0("lam_x", j)
      else paste0("lam_x", j, "_w", w)
      lhs <- paste0(lhs, " + ", lab, "*x", w, "_", j)
    }
    m <- paste0(m, lhs, "\n")

    lhs <- paste0("q", w, " =~ 1*y", w, "_1")
    if (n_indicators_y >= 2) for (j in 2:n_indicators_y) {
      lab <- if (invariance %in% c("metric", "scalar", "strict")) paste0("lam_y", j)
      else paste0("lam_y", j, "_w", w)
      lhs <- paste0(lhs, " + ", lab, "*y", w, "_", j)
    }
    m <- paste0(m, lhs, "\n")
  }

  # Indicator intercepts
  m <- paste0(m, "\n# Indicator intercepts (marker fixed to 0)\n")
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

  # ==================== GROWTH FACTORS (ACCUMULATING; LOAD ON WAVES 2:T) ====================
  m <- paste0(m, "\n# Accumulating growth factors (load on waves 2..T)\n")
  lhs <- paste0("I_x =~ 1*p2")
  if (waves >= 3) for (w in 3:waves) lhs <- paste0(lhs, " + 1*p", w)
  m <- paste0(m, lhs, "\n")
  lhs <- paste0("S_x =~ ", time_scores[2], "*p2")
  if (waves >= 3) for (w in 3:waves) lhs <- paste0(lhs, " + ", time_scores[w], "*p", w)
  m <- paste0(m, lhs, "\n")
  lhs <- paste0("I_y =~ 1*q2")
  if (waves >= 3) for (w in 3:waves) lhs <- paste0(lhs, " + 1*q", w)
  m <- paste0(m, lhs, "\n")
  lhs <- paste0("S_y =~ ", time_scores[2], "*q2")
  if (waves >= 3) for (w in 3:waves) lhs <- paste0(lhs, " + ", time_scores[w], "*q", w)
  m <- paste0(m, lhs, "\n")

  # ==================== AR + CL ====================
  m <- paste0(m, "\n# AR and CL on wave latents\n")
  for (w in 2:waves) {
    ar_x <- if (constrain_beta) "ar_x" else paste0("ar_x", w)
    ar_y <- if (constrain_beta) "ar_y" else paste0("ar_y", w)
    cl_xy <- if (constrain_omega) "cl_xy" else paste0("cl_xy", w)
    cl_yx <- if (constrain_omega) "cl_yx" else paste0("cl_yx", w)
    m <- paste0(m, "p", w, " ~ ", ar_x, "*p", w - 1, " + ", cl_yx, "*q", w - 1, "\n")
    m <- paste0(m, "q", w, " ~ ", ar_y, "*q", w - 1, " + ", cl_xy, "*p", w - 1, "\n")
  }

  # ==================== WAVE 1 EXOGENOUS ====================
  m <- paste0(m, "\n# Wave 1 exogenous\n")
  m <- paste0(m, "p1 ~~ p1\nq1 ~~ q1\np1 ~~ q1\n")
  m <- paste0(m, "p1 ~ 1\nq1 ~ 1\n")

  # ==================== DYNAMIC RESIDUAL VARIANCES ====================
  m <- paste0(m, "\n# Dynamic residual variances/covariances (waves 2..T)\n")
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

  # Latent means w>=2 fixed to 0 (growth factors carry mean structure)
  m <- paste0(m, "\n# Endogenous latent means fixed to 0\n")
  for (w in 2:waves) m <- paste0(m, "p", w, " ~ 0*1\nq", w, " ~ 0*1\n")

  # ==================== GROWTH FACTOR MEANS/VARS/COVS ====================
  m <- paste0(m, "\n# Growth factor means\n")
  m <- paste0(m, "I_x ~ 1\nS_x ~ 1\nI_y ~ 1\nS_y ~ 1\n")

  m <- paste0(m, "\n# Growth factor variances\n")
  m <- paste0(m, "I_x ~~ I_x\nS_x ~~ S_x\nI_y ~~ I_y\nS_y ~~ S_y\n")

  m <- paste0(m, "\n# Growth factor covariances (within variable)\n")
  m <- paste0(m, "I_x ~~ S_x\nI_y ~~ S_y\n")
  m <- paste0(m, "\n# Growth factor covariances (between variables)\n")
  m <- paste0(m, "I_x ~~ I_y\nS_x ~~ S_y\nI_x ~~ S_y\nI_y ~~ S_x\n")

  # Wave-1 covariances with growth factors
  m <- paste0(m, "\n# Wave-1 covariances with growth factors\n")
  for (gf in c("I_x", "S_x", "I_y", "S_y")) {
    m <- paste0(m, "p1 ~~ ", gf, "\n")
    m <- paste0(m, "q1 ~~ ", gf, "\n")
  }

  return(m)
}
