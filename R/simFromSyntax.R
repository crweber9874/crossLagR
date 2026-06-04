#' @title simFromSyntax
#' @description Simulate data from any unified-label crossLagR estimator by
#'   substituting concrete numeric values into the lavaan syntax produced by the
#'   estimator's syntax-generating function, then calling
#'   \code{lavaan::simulateData()}.
#'
#' @details
#' Each crossLagR estimator follows the Usami et al. (2019) unified naming
#' convention: \code{ar_x}, \code{ar_y}, \code{cl_xy}, \code{cl_yx} on the
#' dynamic equation; \code{d_var_x}, \code{d_var_y}, \code{d_cov_xy} on the
#' dynamic residuals. This function takes the estimator name, the chosen number
#' of waves, and concrete population values for those parameters, builds the
#' lavaan syntax, substitutes the labels with fixed numeric values, and
#' simulates data.
#'
#' Other labeled parameters (e.g. growth-factor variances for ALT/LGM/LCMSR,
#' trait variances for TSO, latent FE for BB) are left at the lavaan defaults
#' used by \code{simulateData()}.
#'
#' @param estimator Character. One of "CLPM", "RICLPM", "ALT", "LGM", "LCMSR",
#'   "LCHANGE", "BB", "TSO".
#' @param waves Integer. Number of waves.
#' @param sample_size Integer. Sample size.
#' @param ar_x,ar_y Autoregressive effects (or proportional change for LCS).
#' @param cl_xy,cl_yx Cross-lagged effects (or coupling for LCS).
#' @param d_var_x,d_var_y Dynamic residual variances.
#' @param d_cov_xy Dynamic residual covariance.
#' @param estimator_args Named list. Additional args forwarded to the estimator's
#'   syntax-generating function.
#' @param seed Optional integer seed.
#' @param keep_unified_observed Logical. If TRUE (default), the returned data
#'   frame is renamed so columns x1..xT, y1..yT exist (the conventions used by
#'   the estimateXXX() functions). RICLPM in this package generates data with
#'   columns named x1..xT, y1..yT already; this is a no-op for it.
#'
#' @return data.frame of simulated observed variables (xT, yT columns).
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' dat <- simFromSyntax("ALT", waves = 5, sample_size = 500,
#'                      ar_x = 0.3, ar_y = 0.3,
#'                      cl_xy = 0.2, cl_yx = 0.1)
#' head(dat)
#' }
#'
#' @export
simFromSyntax <- function(estimator,
                          waves = 5,
                          sample_size = 1000,
                          ar_x = 0.3, ar_y = 0.3,
                          cl_xy = 0.1, cl_yx = 0.1,
                          d_var_x = 0.5, d_var_y = 0.5,
                          d_cov_xy = 0,
                          estimator_args = list(),
                          seed = NULL,
                          keep_unified_observed = TRUE) {

  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' is required for simFromSyntax().")
  }

  valid <- c("CLPM", "RICLPM", "ALT", "LGM", "LCMSR", "LCHANGE", "BB", "TSO")
  if (!estimator %in% valid) {
    stop("estimator must be one of: ", paste(valid, collapse = ", "))
  }

  ## Default to bivariate where the estimator supports it, so the DGP
  ## produces both x* and y* columns (required by downstream estimators).
  defaults_per_estimator <- list(
    LGM     = list(variable_type = "bivariate"),
    LCHANGE = list(variable_type = "bivariate")
  )
  args <- c(list(waves = waves),
            defaults_per_estimator[[estimator]],
            as.list(estimator_args))
  args <- args[!duplicated(names(args))]

  fn <- switch(
    estimator,
    "CLPM"    = crossLagR::estimateCLPM,
    "RICLPM"  = crossLagR::estimateRICLPM,
    "ALT"     = crossLagR::estimateALT,
    "LGM"     = crossLagR::estimateLGM,
    "LCMSR"   = crossLagR::estimateLCMSR,
    "LCHANGE" = crossLagR::estimateLChange,
    "BB"      = crossLagR::estimateBollen_and_Brand,
    "TSO"     = crossLagR::estimateTSO
  )

  args <- args[intersect(names(args), names(formals(fn)))]
  syntax <- do.call(fn, args)

  populated <- populate_unified_labels(
    syntax,
    ar_x = ar_x, ar_y = ar_y,
    cl_xy = cl_xy, cl_yx = cl_yx,
    d_var_x = d_var_x, d_var_y = d_var_y,
    d_cov_xy = d_cov_xy
  )

  if (!is.null(seed)) set.seed(seed)

  ## simulateData supports a model with a mix of fixed and free parameters.
  ## sample.nobs sets N; meanstructure = TRUE since some estimators (RICLPM)
  ## free trait means.
  dat <- lavaan::simulateData(
    model        = populated,
    sample.nobs  = sample_size,
    meanstructure = TRUE
  )

  ## Most estimateXXX() functions use x1..xT, y1..yT. Sanity-check + return.
  if (keep_unified_observed) {
    needed <- c(paste0("x", seq_len(waves)), paste0("y", seq_len(waves)))
    missing_obs <- setdiff(needed, names(dat))
    if (length(missing_obs)) {
      warning("simFromSyntax: generated data missing expected columns: ",
              paste(missing_obs, collapse = ", "))
    }
  }

  dat
}

#' @title populate_unified_labels
#' @description Substitute concrete numeric values into a unified-label lavaan
#'   syntax string. Replaces every occurrence of \code{label*} with
#'   \code{<value>*}, which lavaan interprets as a fixed parameter.
#' @param syntax Character. The unmodified lavaan syntax.
#' @param ar_x,ar_y,cl_xy,cl_yx Numeric values for AR/CL labels.
#' @param d_var_x,d_var_y,d_cov_xy Numeric values for dynamic-residual labels.
#' @return Character. The populated syntax.
#' @export
populate_unified_labels <- function(syntax,
                                    ar_x = 0.3, ar_y = 0.3,
                                    cl_xy = 0.1, cl_yx = 0.1,
                                    d_var_x = 0.5, d_var_y = 0.5,
                                    d_cov_xy = 0) {

  fmt <- function(x) sprintf("%.6g", x)

  ## Order matters only so that long names match before short ones; using
  ## word-boundary regex avoids partial matches (e.g. "ar_x" not matching "ar_xx").
  subs <- list(
    "\\bar_x\\*"     = paste0(fmt(ar_x), "*"),
    "\\bar_y\\*"     = paste0(fmt(ar_y), "*"),
    "\\bcl_xy\\*"    = paste0(fmt(cl_xy), "*"),
    "\\bcl_yx\\*"    = paste0(fmt(cl_yx), "*"),
    "\\bd_var_x\\*"  = paste0(fmt(d_var_x), "*"),
    "\\bd_var_y\\*"  = paste0(fmt(d_var_y), "*"),
    "\\bd_cov_xy\\*" = paste0(fmt(d_cov_xy), "*")
  )

  s <- syntax
  for (pat in names(subs)) {
    s <- gsub(pat, subs[[pat]], s, perl = TRUE)
  }
  s
}
