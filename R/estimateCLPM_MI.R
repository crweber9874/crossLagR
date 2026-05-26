#' @title estimateCLPM_MI
#' @description Generates lavaan syntax for a multi-indicator Cross-Lagged Panel
#'   Model (CLPM) following Usami, Murayama, & Hamaker (2019).
#'
#' @param waves Integer. Number of waves (>= 2).
#' @param n_indicators_x Integer. Indicators for X at each wave (>= 1).
#' @param n_indicators_y Integer. Indicators for Y at each wave (>= 1).
#' @param invariance Character. One of \code{"configural"}, \code{"metric"},
#'   \code{"scalar"}, \code{"strict"}. Default \code{"scalar"}.
#' @param correlated_residuals Logical. Correlate same-item residuals across
#'   waves. Default \code{TRUE}.
#' @param constrain_beta,constrain_omega,constrain_residual_variances,constrain_residual_covariances Logical.
#'   Equality constraints across waves on AR, CL, dynamic residual variances
#'   and covariances respectively.
#' @param estimate_means Logical. Estimate latent factor means at wave 1.
#' @param start_values Logical. Provide starting values for AR/CL.
#'
#' @return Character string of lavaan model syntax.
#'
#' @details
#' Variables follow naming \code{x{w}_{j}}, \code{y{w}_{j}}. The first indicator
#' is the marker (loading fixed to 1, intercept fixed to 0). With
#' \code{invariance = "scalar"} the remaining loadings and intercepts are
#' constrained equal across waves. \code{"strict"} additionally constrains
#' indicator residual variances across waves. Same-item residuals are
#' correlated across all wave-pairs when \code{correlated_residuals = TRUE} to
#' absorb stable item-specific variance.
#'
#' @references
#' Usami, S., Murayama, K., & Hamaker, E. L. (2019). \emph{Psychological
#'   Methods}, 24(5), 637-657.
#'
#' @examples
#' \dontrun{
#' syntax <- estimateCLPM_MI(waves = 4, n_indicators_x = 3, n_indicators_y = 3)
#' fit <- lavaan::lavaan(syntax, data = my_data, meanstructure = TRUE,
#'                       int.ov.free = TRUE, int.lv.free = TRUE)
#' }
#'
#' @export
estimateCLPM_MI <- function(waves = 5,
                            n_indicators_x = 3,
                            n_indicators_y = 3,
                            invariance = c("scalar", "metric", "configural", "strict"),
                            correlated_residuals = TRUE,
                            constrain_beta = TRUE,
                            constrain_omega = TRUE,
                            constrain_residual_variances = TRUE,
                            constrain_residual_covariances = TRUE,
                            estimate_means = TRUE,
                            start_values = FALSE) {

  if (!is.numeric(waves) || waves < 2 || waves != as.integer(waves)) {
    stop("Error: 'waves' must be an integer >= 2.")
  }
  if (!is.numeric(n_indicators_x) || n_indicators_x < 1 ||
      n_indicators_x != as.integer(n_indicators_x)) {
    stop("Error: 'n_indicators_x' must be a positive integer.")
  }
  if (!is.numeric(n_indicators_y) || n_indicators_y < 1 ||
      n_indicators_y != as.integer(n_indicators_y)) {
    stop("Error: 'n_indicators_y' must be a positive integer.")
  }
  invariance <- match.arg(invariance)

  m <- ""

  # ==================== MEASUREMENT MODEL ====================
  m <- paste0(m, "# Measurement model (latent factors p_w, q_w)\n")
  for (w in 1:waves) {
    # p_w =~ 1*x{w}_1 + lam_x{j}*x{w}_{j} + ...
    lhs <- paste0("p", w, " =~ 1*x", w, "_1")
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

    lhs <- paste0("q", w, " =~ 1*y", w, "_1")
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
      lab <- if (invariance == "strict") {
        paste0("u_var_x", j)
      } else {
        paste0("u_var_x", j, "_w", w)
      }
      m <- paste0(m, "x", w, "_", j, " ~~ ", lab, "*x", w, "_", j, "\n")
    }
    for (j in 1:n_indicators_y) {
      lab <- if (invariance == "strict") {
        paste0("u_var_y", j)
      } else {
        paste0("u_var_y", j, "_w", w)
      }
      m <- paste0(m, "y", w, "_", j, " ~~ ", lab, "*y", w, "_", j, "\n")
    }
  }

  # ==================== CORRELATED RESIDUALS ACROSS WAVES ====================
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

  # ==================== LATENT MEAN STRUCTURE ====================
  m <- paste0(m, "\n# Latent factor means\n")
  if (estimate_means) {
    m <- paste0(m, "p1 ~ 1\nq1 ~ 1\n")
    if (waves >= 2) {
      for (w in 2:waves) {
        m <- paste0(m, "p", w, " ~ 0*1\nq", w, " ~ 0*1\n")
      }
    }
  } else {
    for (w in 1:waves) {
      m <- paste0(m, "p", w, " ~ 0*1\nq", w, " ~ 0*1\n")
    }
  }

  # ==================== AR + CL ====================
  m <- paste0(m, "\n# Autoregressive and cross-lagged effects\n")
  s_ar <- if (start_values) "start(0.3)*" else ""
  s_cl <- if (start_values) "start(0.1)*" else ""
  for (w in 2:waves) {
    ar_y <- if (constrain_beta) "ar_y" else paste0("ar_y", w)
    ar_x <- if (constrain_beta) "ar_x" else paste0("ar_x", w)
    cl_xy <- if (constrain_omega) "cl_xy" else paste0("cl_xy", w)
    cl_yx <- if (constrain_omega) "cl_yx" else paste0("cl_yx", w)
    m <- paste0(m, "p", w, " ~ ", s_ar, ar_y, "*p", w - 1,
                " + ", s_cl, cl_xy, "*q", w - 1, "\n")
    m <- paste0(m, "q", w, " ~ ", s_ar, ar_x, "*q", w - 1,
                " + ", s_cl, cl_yx, "*p", w - 1, "\n")
  }

  # ==================== DYNAMIC (LATENT) RESIDUAL VARIANCES ====================
  m <- paste0(m, "\n# Dynamic residual variances on p, q\n")
  if (constrain_residual_variances) {
    m <- paste0(m, "p1 ~~ d_var_y1*p1\nq1 ~~ d_var_x1*q1\n")
    for (w in 2:waves) {
      m <- paste0(m, "p", w, " ~~ d_var_y*p", w, "\n")
      m <- paste0(m, "q", w, " ~~ d_var_x*q", w, "\n")
    }
  } else {
    for (w in 1:waves) {
      m <- paste0(m, "p", w, " ~~ p", w, "\nq", w, " ~~ q", w, "\n")
    }
  }

  # ==================== DYNAMIC RESIDUAL COVARIANCES ====================
  m <- paste0(m, "\n# Within-time latent residual covariances\n")
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

  return(m)
}
