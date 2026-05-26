#' @title estimateLCMSR. This is just the GCLPM. Structure the means here.
#' @description Generates lavaan model syntax for the bivariate Latent Curve Model
#'   with Structured Residuals (LCM-SR), also known as the General Cross-Lagged
#'   Model (GCLM).
#'
#' @param waves Integer. Number of waves (time points). Must be >= 3.
#' @param time_scores Numeric vector of time scores for slope factor loadings.
#'   If NULL, uses 0, 1, 2, ..., waves-1. Default is NULL.
#' @param constrain_beta Logical. If TRUE, constrains autoregressive effects to
#'   equality across waves. Default is TRUE.
#' @param constrain_omega Logical. If TRUE, constrains cross-lagged effects to
#'   equality across waves. Default is TRUE.
#' @param constrain_residual_variances Logical. If TRUE, constrains dynamic residual
#'   variances to equality across waves. Default is TRUE.
#' @param constrain_residual_covariances Logical. If TRUE, constrains dynamic residual
#'   covariances to equality across waves. Default is TRUE.
#' @param estimate_quadratic Logical. If TRUE, includes quadratic growth factors.
#'   Requires >= 4 waves. Default is FALSE.
#'
#' @return A character string containing lavaan model syntax.
#'
#' @details
#' The LCM-SR (also called the GCLM or RI-CLPM with structured means) combines
#' latent growth curves with autoregressive and cross-lagged dynamics in the
#' \emph{residuals} from those curves, following the unified framework of
#' Usami, Murayama, & Hamaker (2019).
#'
#' Unlike the ALT model (where AR/CL dynamics operate on the scores directly and
#' growth factors are accumulating), the LCM-SR first removes the growth
#' trajectory, then models the deviations with AR/CL dynamics. This means the
#' growth factors (I, S) have \emph{only direct effects} on scores -- they are
#' true growth factors, not accumulating factors.
#'
#' Parameter labels follow the unified naming convention:
#' \itemize{
#'   \item \code{ar_x}, \code{ar_y}: Autoregressive effects on within-person
#'     deviations from the growth trajectory.
#'   \item \code{cl_xy}, \code{cl_yx}: Cross-lagged effects on within-person
#'     deviations.
#'   \item \code{I_x}, \code{I_y}: Growth intercept factors (direct effects only).
#'   \item \code{S_x}, \code{S_y}: Growth slope factors (direct effects only).
#'   \item \code{d_var_x}, \code{d_var_y}: Dynamic residual variances.
#'   \item \code{d_cov_xy}: Dynamic residual covariance.
#' }
#'
#' Within-person deviations from the growth curve are labeled \code{p} (for X)
#' and \code{q} (for Y), analogous to the RI-CLPM's within-person factors.
#'
#' @references
#' Curran, P. J., Howard, A. L., Bainter, S. A., Lane, S. T., & McGinley, J. S.
#'   (2014). The separation of between-person and within-person components of
#'   individual change over time. \emph{Annual Review of Psychology}, 65, 637-660.
#'
#' Berry, D., & Willoughby, M. T. (2017). On the practical interpretability of
#'   cross-lagged panel models. \emph{Structural Equation Modeling}, 24(3), 455-468.
#'
#' Usami, S., Murayama, K., & Hamaker, E. L. (2019). A unified framework of
#'   longitudinal models to examine reciprocal relations. \emph{Psychological
#'   Methods}, 24(5), 637-657.
#'
#' @examples
#' \dontrun{
#' syntax <- estimateLCMSR(waves = 5)
#' cat(syntax)
#'
#' library(lavaan)
#' fit <- lavaan::sem(syntax, data = my_data, meanstructure = TRUE)
#' summary(fit, fit.measures = TRUE, standardized = TRUE)
#' }
#'
#' @export
estimateLCMSR <- function(waves = 5,
                          time_scores = NULL,
                          constrain_beta = TRUE,
                          constrain_omega = TRUE,
                          constrain_residual_variances = TRUE,
                          constrain_residual_covariances = TRUE,
                          estimate_quadratic = FALSE) {

  if (!is.numeric(waves) || waves <= 0 || waves != as.integer(waves)) {
    stop("Error: Parameter 'waves' must be a positive integer.")
  }
  if (waves < 3) {
    stop("Error: LCM-SR model requires at least 3 waves.")
  }
  if (estimate_quadratic && waves < 4) {
    stop("Error: Quadratic LCM-SR requires at least 4 waves.")
  }

  # Time scores
  if (is.null(time_scores)) {
    time_scores <- 0:(waves - 1)
  } else {
    if (!is.numeric(time_scores) || length(time_scores) != waves) {
      stop("Error: 'time_scores' must be a numeric vector of length equal to 'waves'.")
    }
  }

  if (estimate_quadratic) {
    quad_scores <- time_scores^2
  }

  model <- ""

  # ==================== OBSERVED TO LATENT ====================
  model <- paste0(model, "# Within-person latent deviations\n")
  for (w in 1:waves) {
    model <- paste0(model, "p", w, " =~ 1*x", w, "\n")
    model <- paste0(model, "q", w, " =~ 1*y", w, "\n")
  }

  # Fix observed residuals to 0 (perfect indicators)
  model <- paste0(model, "\n# Fix observed residuals to zero\n")
  for (w in 1:waves) {
    model <- paste0(model, "x", w, " ~~ 0*x", w, "\n")
    model <- paste0(model, "y", w, " ~~ 0*y", w, "\n")
  }

  # Fix observed intercepts to 0
  model <- paste0(model, "\n# Fix observed intercepts to zero\n")
  for (w in 1:waves) {
    model <- paste0(model, "x", w, " ~ 0*1\n")
    model <- paste0(model, "y", w, " ~ 0*1\n")
  }

  # ==================== GROWTH FACTORS (DIRECT EFFECTS ON OBSERVED) ====================
  # Growth factors load on OBSERVED variables, not on p/q. This ensures that
  # p/q capture structured residuals (deviations from growth), and AR/CL
  # operates on those residuals. Per Usami et al. (2019) Eq. 11, growth factors
  # have only direct effects on scores -- they do NOT accumulate through the
  # AR process (unlike the ALT model).
  model <- paste0(model, "\n# Growth factors (intercept and slope, direct effects on observed)\n")

  # I_x: unit loadings on all x
  model <- paste0(model, "I_x =~ 1*x1")
  for (w in 2:waves) {
    model <- paste0(model, " + 1*x", w)
  }
  model <- paste0(model, "\n")

  # S_x: time loadings on all x
  model <- paste0(model, "S_x =~ ", time_scores[1], "*x1")
  for (w in 2:waves) {
    model <- paste0(model, " + ", time_scores[w], "*x", w)
  }
  model <- paste0(model, "\n")

  # Q_x if quadratic
  if (estimate_quadratic) {
    model <- paste0(model, "Q_x =~ ", quad_scores[1], "*x1")
    for (w in 2:waves) {
      model <- paste0(model, " + ", quad_scores[w], "*x", w)
    }
    model <- paste0(model, "\n")
  }

  # I_y: unit loadings on all y
  model <- paste0(model, "I_y =~ 1*y1")
  for (w in 2:waves) {
    model <- paste0(model, " + 1*y", w)
  }
  model <- paste0(model, "\n")

  # S_y: time loadings on all y
  model <- paste0(model, "S_y =~ ", time_scores[1], "*y1")
  for (w in 2:waves) {
    model <- paste0(model, " + ", time_scores[w], "*y", w)
  }
  model <- paste0(model, "\n")

  if (estimate_quadratic) {
    model <- paste0(model, "Q_y =~ ", quad_scores[1], "*y1")
    for (w in 2:waves) {
      model <- paste0(model, " + ", quad_scores[w], "*y", w)
    }
    model <- paste0(model, "\n")
  }

  # ==================== STRUCTURED RESIDUALS: AR AND CL ====================
  model <- paste0(model, "\n# Autoregressive and cross-lagged paths (structured residuals)\n")
  for (w in 2:waves) {
    if (constrain_beta && constrain_omega) {
      model <- paste0(model, "p", w, " ~ ar_x*p", w - 1, " + cl_yx*q", w - 1, "\n")
      model <- paste0(model, "q", w, " ~ ar_y*q", w - 1, " + cl_xy*p", w - 1, "\n")
    } else if (constrain_beta && !constrain_omega) {
      model <- paste0(model, "p", w, " ~ ar_x*p", w - 1, " + cl_yx", w, "*q", w - 1, "\n")
      model <- paste0(model, "q", w, " ~ ar_y*q", w - 1, " + cl_xy", w, "*p", w - 1, "\n")
    } else if (!constrain_beta && constrain_omega) {
      model <- paste0(model, "p", w, " ~ ar_x", w, "*p", w - 1, " + cl_yx*q", w - 1, "\n")
      model <- paste0(model, "q", w, " ~ ar_y", w, "*q", w - 1, " + cl_xy*p", w - 1, "\n")
    } else {
      model <- paste0(model, "p", w, " ~ ar_x", w, "*p", w - 1, " + cl_yx", w, "*q", w - 1, "\n")
      model <- paste0(model, "q", w, " ~ ar_y", w, "*q", w - 1, " + cl_xy", w, "*p", w - 1, "\n")
    }
  }

  # ==================== DYNAMIC RESIDUAL VARIANCES ====================
  model <- paste0(model, "\n# Dynamic residual variances and covariances\n")
  # Wave 1 variances (free, initial conditions after removing growth)
  model <- paste0(model, "p1 ~~ p1\n")
  model <- paste0(model, "q1 ~~ q1\n")
  model <- paste0(model, "p1 ~~ q1\n")

  if (constrain_residual_variances) {
    for (w in 2:waves) {
      model <- paste0(model, "p", w, " ~~ d_var_x*p", w, "\n")
      model <- paste0(model, "q", w, " ~~ d_var_y*q", w, "\n")
    }
  } else {
    for (w in 2:waves) {
      model <- paste0(model, "p", w, " ~~ p", w, "\n")
      model <- paste0(model, "q", w, " ~~ q", w, "\n")
    }
  }

  if (constrain_residual_covariances) {
    for (w in 2:waves) {
      model <- paste0(model, "p", w, " ~~ d_cov_xy*q", w, "\n")
    }
  } else {
    for (w in 2:waves) {
      model <- paste0(model, "p", w, " ~~ q", w, "\n")
    }
  }

  # Fix latent means to 0
  model <- paste0(model, "\n# Latent means fixed to zero\n")
  for (w in 1:waves) {
    model <- paste0(model, "p", w, " ~ 0*1\n")
    model <- paste0(model, "q", w, " ~ 0*1\n")
  }

  # ==================== GROWTH FACTOR MEANS ====================
  model <- paste0(model, "\n# Growth factor means\n")
  model <- paste0(model, "I_x ~ 1\n")
  model <- paste0(model, "S_x ~ 1\n")
  model <- paste0(model, "I_y ~ 1\n")
  model <- paste0(model, "S_y ~ 1\n")
  if (estimate_quadratic) {
    model <- paste0(model, "Q_x ~ 1\n")
    model <- paste0(model, "Q_y ~ 1\n")
  }

  # ==================== GROWTH FACTOR VARIANCES ====================
  model <- paste0(model, "\n# Growth factor variances\n")
  model <- paste0(model, "I_x ~~ I_x\n")
  model <- paste0(model, "S_x ~~ S_x\n")
  model <- paste0(model, "I_y ~~ I_y\n")
  model <- paste0(model, "S_y ~~ S_y\n")
  if (estimate_quadratic) {
    model <- paste0(model, "Q_x ~~ Q_x\n")
    model <- paste0(model, "Q_y ~~ Q_y\n")
  }

  # ==================== GROWTH FACTOR COVARIANCES ====================
  model <- paste0(model, "\n# Growth factor covariances (within variable)\n")
  model <- paste0(model, "I_x ~~ S_x\n")
  model <- paste0(model, "I_y ~~ S_y\n")
  if (estimate_quadratic) {
    model <- paste0(model, "I_x ~~ Q_x\n")
    model <- paste0(model, "S_x ~~ Q_x\n")
    model <- paste0(model, "I_y ~~ Q_y\n")
    model <- paste0(model, "S_y ~~ Q_y\n")
  }

  model <- paste0(model, "\n# Growth factor covariances (between variable)\n")
  model <- paste0(model, "I_x ~~ I_y\n")
  model <- paste0(model, "S_x ~~ S_y\n")
  model <- paste0(model, "I_x ~~ S_y\n")
  model <- paste0(model, "I_y ~~ S_x\n")
  if (estimate_quadratic) {
    model <- paste0(model, "Q_x ~~ Q_y\n")
    model <- paste0(model, "I_x ~~ Q_y\n")
    model <- paste0(model, "I_y ~~ Q_x\n")
    model <- paste0(model, "S_x ~~ Q_y\n")
    model <- paste0(model, "S_y ~~ Q_x\n")
  }

  return(model)
}
