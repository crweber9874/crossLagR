#' @title estimateALT, following Usami, Murayama, & Hamaker (2019) unified frameowk.
#' @description Generates lavaan model syntax for the bivariate Autoregressive Latent
#'   Trajectory (ALT) model.
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
#'
#' @return A character string containing lavaan model syntax.
#'
#' @details
#' The ALT model combines latent growth trajectories (intercept I, slope S) with
#' autoregressive and cross-lagged dynamics, following the unified framework of
#' Usami, Murayama, & Hamaker (2019). Unlike the LCM-SR (which models dynamics
#' in the residuals from growth curves), the ALT model applies AR/CL dynamics
#' directly to the latent scores. The growth factors in the ALT act as
#' \emph{accumulating factors} -- they have both direct effects on scores AND
#' indirect effects through the lagged dynamics.
#'
#' Parameter labels follow the unified naming convention:
#' \itemize{
#'   \item \code{ar_x}, \code{ar_y}: Autoregressive effects on latent scores.
#'     These reflect carry-over that compounds with the growth trajectory.
#'   \item \code{cl_xy}, \code{cl_yx}: Cross-lagged effects on latent scores.
#'   \item \code{I_x}, \code{I_y}: Latent intercept factors (accumulating).
#'   \item \code{S_x}, \code{S_y}: Latent slope factors (accumulating).
#'   \item \code{d_var_x}, \code{d_var_y}: Dynamic residual variances.
#'   \item \code{d_cov_xy}: Dynamic residual covariance.
#' }
#'
#' The ALT conditions on the first observation (wave 1 is exogenous). Growth
#' factor loadings start from wave 2.
#'
#' @references
#' Bollen, K. A., & Curran, P. J. (2006). \emph{Latent curve models: A structural
#'   equation perspective}. Wiley.
#'
#' Curran, P. J., & Bollen, K. A. (2001). The best of both worlds: Combining
#'   autoregressive and latent curve models. In L. M. Collins & A. G. Sayer (Eds.),
#'   \emph{New methods for the analysis of change} (pp. 107-135). APA.
#'
#' Usami, S., Murayama, K., & Hamaker, E. L. (2019). A unified framework of
#'   longitudinal models to examine reciprocal relations. \emph{Psychological
#'   Methods}, 24(5), 637-657.
#'
#' @examples
#' \dontrun{
#' syntax <- estimateALT(waves = 5)
#' cat(syntax)
#'
#' library(lavaan)
#' fit <- lavaan::sem(syntax, data = my_data)
#' summary(fit, fit.measures = TRUE, standardized = TRUE)
#' }
#'
#' @export
estimateALT <- function(waves = 5,
                        time_scores = NULL,
                        constrain_beta = TRUE,
                        constrain_omega = TRUE,
                        constrain_residual_variances = TRUE,
                        constrain_residual_covariances = TRUE) {

  if (!is.numeric(waves) || waves <= 0 || waves != as.integer(waves)) {
    stop("Error: Parameter 'waves' must be a positive integer.")
  }
  if (waves < 3) {
    stop("Error: ALT model requires at least 3 waves.")
  }

  # Time scores for slope loadings (wave 2 onward)
  if (is.null(time_scores)) {
    time_scores <- 0:(waves - 1)
  } else {
    if (!is.numeric(time_scores) || length(time_scores) != waves) {
      stop("Error: 'time_scores' must be a numeric vector of length equal to 'waves'.")
    }
  }

  model <- ""

  # ==================== LATENT VARIABLES ====================
  model <- paste0(model, "# Latent variables (perfect indicators)\n")
  for (w in 1:waves) {
    model <- paste0(model, "p", w, " =~ 1*x", w, "\n")
    model <- paste0(model, "q", w, " =~ 1*y", w, "\n")
  }

  # Fix observed variable residuals to 0
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

  # ==================== GROWTH FACTORS (ACCUMULATING) ====================
  # I and S load on endogenous waves (2:waves) only; wave 1 is exogenous
  model <- paste0(model, "\n# Growth factors (accumulating, load on waves 2:T)\n")

  # I_x: unit loadings on p2, p3, ..., pT
  model <- paste0(model, "I_x =~ 1*p2")
  for (w in 3:waves) {
    model <- paste0(model, " + 1*p", w)
  }
  model <- paste0(model, "\n")

  # S_x: time-varying loadings on p2, p3, ..., pT
  model <- paste0(model, "S_x =~ ", time_scores[2], "*p2")
  for (w in 3:waves) {
    model <- paste0(model, " + ", time_scores[w], "*p", w)
  }
  model <- paste0(model, "\n")

  # I_y: unit loadings on q2, q3, ..., qT
  model <- paste0(model, "I_y =~ 1*q2")
  for (w in 3:waves) {
    model <- paste0(model, " + 1*q", w)
  }
  model <- paste0(model, "\n")

  # S_y: time-varying loadings on q2, q3, ..., qT
  model <- paste0(model, "S_y =~ ", time_scores[2], "*q2")
  for (w in 3:waves) {
    model <- paste0(model, " + ", time_scores[w], "*q", w)
  }
  model <- paste0(model, "\n")

  # ==================== AUTOREGRESSIVE AND CROSS-LAGGED ====================
  model <- paste0(model, "\n# Autoregressive and cross-lagged paths\n")
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

  # ==================== EXOGENOUS WAVE 1 ====================
  model <- paste0(model, "\n# Exogenous wave 1 variances, covariance, and means\n")
  model <- paste0(model, "p1 ~~ p1\n")
  model <- paste0(model, "q1 ~~ q1\n")
  model <- paste0(model, "p1 ~~ q1\n")
  model <- paste0(model, "p1 ~ 1\n")
  model <- paste0(model, "q1 ~ 1\n")

  # ==================== DYNAMIC RESIDUAL VARIANCES ====================
  model <- paste0(model, "\n# Dynamic residual variances and covariances\n")
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

  # Fix latent means to 0 for endogenous waves
  model <- paste0(model, "\n# Fix endogenous latent means to zero\n")
  for (w in 2:waves) {
    model <- paste0(model, "p", w, " ~ 0*1\n")
    model <- paste0(model, "q", w, " ~ 0*1\n")
  }

  # ==================== GROWTH FACTOR VARIANCES AND COVARIANCES ====================
  model <- paste0(model, "\n# Growth factor means\n")
  model <- paste0(model, "I_x ~ 1\n")
  model <- paste0(model, "S_x ~ 1\n")
  model <- paste0(model, "I_y ~ 1\n")
  model <- paste0(model, "S_y ~ 1\n")

  model <- paste0(model, "\n# Growth factor variances\n")
  model <- paste0(model, "I_x ~~ I_x\n")
  model <- paste0(model, "S_x ~~ S_x\n")
  model <- paste0(model, "I_y ~~ I_y\n")
  model <- paste0(model, "S_y ~~ S_y\n")

  model <- paste0(model, "\n# Growth factor covariances (within variable)\n")
  model <- paste0(model, "I_x ~~ S_x\n")
  model <- paste0(model, "I_y ~~ S_y\n")

  model <- paste0(model, "\n# Growth factor covariances (between variable)\n")
  model <- paste0(model, "I_x ~~ I_y\n")
  model <- paste0(model, "S_x ~~ S_y\n")
  model <- paste0(model, "I_x ~~ S_y\n")
  model <- paste0(model, "I_y ~~ S_x\n")

  # ==================== EXOGENOUS-GROWTH COVARIANCES ====================
  model <- paste0(model, "\n# Exogenous wave 1 covariances with growth factors\n")
  model <- paste0(model, "p1 ~~ I_x\n")
  model <- paste0(model, "p1 ~~ S_x\n")
  model <- paste0(model, "p1 ~~ I_y\n")
  model <- paste0(model, "p1 ~~ S_y\n")
  model <- paste0(model, "q1 ~~ I_x\n")
  model <- paste0(model, "q1 ~~ S_x\n")
  model <- paste0(model, "q1 ~~ I_y\n")
  model <- paste0(model, "q1 ~~ S_y\n")

  return(model)
}
