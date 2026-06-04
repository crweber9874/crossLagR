#' @title diagnoseFit
#' @description Inspect a fitted lavaan object for Heywood cases and other
#'   empirical-identification boundary conditions common in panel SEM. Returns
#'   a structured diagnosis with plain-language explanations and remedial
#'   suggestions, and (by default) emits a single condensed warning when any
#'   issues are detected.
#'
#' @param fit A fitted lavaan object (from \code{lavaan::lavaan} or
#'   \code{lavaan::sem}).
#' @param warnings_seen Optional character vector of warning messages captured
#'   during fitting (e.g., via \code{withCallingHandlers}). Some boundary
#'   conditions (NPD covariance, vcov singular) are detectable only from
#'   lavaan's warning text. \code{\link{lavaanDiagnose}} captures these
#'   automatically.
#' @param quiet Logical. If \code{FALSE} (default), emits a single condensed
#'   \code{warning()} summarizing detected issues. Set \code{TRUE} for silent
#'   programmatic use (e.g., inside Monte Carlo loops).
#' @param tol_negative Numeric. Variance estimates below this are flagged as
#'   negative. Default \code{-1e-6} (tiny numerical noise treated as zero).
#' @param tol_std Numeric. Standardized loadings or correlations with
#'   absolute value above this are flagged. Default \code{1.0}.
#'
#' @return An object of class \code{crossLagR_diagnosis} with elements:
#'   \itemize{
#'     \item \code{issues}: list of detected issues; each has \code{code},
#'       \code{severity}, \code{parameters}, \code{message},
#'       \code{explanation}, \code{remedies}.
#'     \item \code{n_issues}: integer count.
#'     \item \code{codes}: character vector of detected codes.
#'     \item \code{summary}: one-line summary string.
#'   }
#'
#' @details
#' Detected issue codes:
#' \itemize{
#'   \item \code{non_convergence}: lavaan optimizer did not converge.
#'   \item \code{negative_variance}: Heywood case proper — a variance estimate
#'     is negative.
#'   \item \code{std_loading_gt_one}: standardized factor loading > 1.
#'   \item \code{std_correlation_gt_one}: latent correlation > 1
#'     (implies a non-PD latent covariance matrix).
#'   \item \code{npd_covariance}: lavaan flagged a non-positive-definite
#'     model-implied or sample covariance matrix.
#'   \item \code{vcov_not_pd}: parameter Hessian/vcov is not PD, so standard
#'     errors are unreliable.
#'   \item \code{iter_limit}: optimizer reached iteration cap without
#'     converging.
#' }
#'
#' Call \code{\link{heywoodHelp}()} to print the full catalog without a fit.
#'
#' @seealso \code{\link{lavaanDiagnose}}, \code{\link{heywoodHelp}}.
#' @export
diagnoseFit <- function(fit,
                        warnings_seen = NULL,
                        quiet = FALSE,
                        tol_negative = -1e-6,
                        tol_std = 1.0) {

  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' is required for diagnoseFit().")
  }
  if (!inherits(fit, "lavaan")) {
    stop("'fit' must be a fitted lavaan object.")
  }

  issues <- list()
  add <- function(iss) issues[[length(issues) + 1L]] <<- iss

  ## ---- convergence ---------------------------------------------------------
  converged <- isTRUE(lavaan::lavInspect(fit, "converged"))
  if (!converged) add(.fill_params(.HEYWOOD_CATALOG$non_convergence))

  ## ---- parameter-estimate based checks -------------------------------------
  pe <- tryCatch(
    lavaan::parameterEstimates(fit, standardized = TRUE),
    error = function(e) NULL
  )

  if (!is.null(pe)) {
    var_rows <- pe[pe$op == "~~" & pe$lhs == pe$rhs, , drop = FALSE]
    neg <- var_rows[is.finite(var_rows$est) & var_rows$est < tol_negative, ,
                    drop = FALSE]
    if (nrow(neg) > 0) {
      add(.fill_params(.HEYWOOD_CATALOG$negative_variance,
                       params = unique(neg$lhs)))
    }

    if ("std.all" %in% names(pe)) {
      load_rows <- pe[pe$op == "=~" & is.finite(pe$std.all), , drop = FALSE]
      big <- load_rows[abs(load_rows$std.all) > tol_std, , drop = FALSE]
      if (nrow(big) > 0) {
        add(.fill_params(.HEYWOOD_CATALOG$std_loading_gt_one,
                         params = paste0(big$lhs, "=~", big$rhs)))
      }

      cor_rows <- pe[pe$op == "~~" & pe$lhs != pe$rhs &
                     is.finite(pe$std.all), , drop = FALSE]
      big_r <- cor_rows[abs(cor_rows$std.all) > tol_std, , drop = FALSE]
      if (nrow(big_r) > 0) {
        add(.fill_params(.HEYWOOD_CATALOG$std_correlation_gt_one,
                         params = paste0(big_r$lhs, "~~", big_r$rhs)))
      }
    }
  }

  ## ---- warning-text based checks -------------------------------------------
  has_code <- function(code) {
    any(vapply(issues, function(x) identical(x$code, code), logical(1)))
  }
  if (length(warnings_seen)) {
    .has <- function(pat) any(grepl(pat, warnings_seen, ignore.case = TRUE))

    if (.has("not positive[- ]definite|positive definite") &&
        !has_code("negative_variance") &&
        !has_code("std_correlation_gt_one")) {
      add(.fill_params(.HEYWOOD_CATALOG$npd_covariance))
    }
    if (.has("variance[- ]covariance matrix of the estimated parameters|vcov")) {
      add(.fill_params(.HEYWOOD_CATALOG$vcov_not_pd))
    }
    if (.has("reached the iteration limit|did not converge") &&
        !has_code("non_convergence")) {
      add(.fill_params(.HEYWOOD_CATALOG$iter_limit))
    }
  }

  codes <- vapply(issues, function(x) x$code, character(1))
  summary_str <- if (length(issues) == 0L) {
    "No Heywood-class issues detected."
  } else {
    paste0(length(issues), " issue(s): ", paste(codes, collapse = ", "), ".")
  }

  out <- structure(
    list(issues = issues, n_issues = length(issues),
         codes = codes, summary = summary_str),
    class = "crossLagR_diagnosis"
  )

  if (!quiet && out$n_issues > 0L) {
    warning(
      "crossLagR: ", out$n_issues, " Heywood-class issue(s) detected (",
      paste(codes, collapse = ", "),
      "). Use print() on the returned diagnosis for explanations and ",
      "remedies, or heywoodHelp() for the full catalog.",
      call. = FALSE
    )
  }

  out
}


#' @title lavaanDiagnose
#' @description Fit a lavaan model and bundle the fit with three crossLagR
#'   diagnostic components: (1) a labelled parameter summary
#'   (\code{\link{crossLagRSummary}}), (2) a Heywood-case diagnosis
#'   (\code{\link{diagnoseFit}}), and (3) an optional between/within ICC
#'   report (\code{\link{iccReport}}). Captures lavaan warnings so the
#'   diagnosis can include warning-only issues like NPD covariances, and
#'   emits a single condensed warning when issues or high ICCs are detected.
#'
#' @param model A lavaan model syntax string.
#' @param data A data frame.
#' @param fit_function Either \code{"lavaan"} (default) or \code{"sem"}.
#' @param icc_vars Optional character vector of column stems (e.g.,
#'   \code{c("x", "y")}) for computing the data-level ICC via
#'   \code{\link{iccReport}}. If \code{NULL} (default), ICC is skipped.
#' @param icc_threshold Numeric. ICC values above this trigger the
#'   crossLagR Shiny-sensitivity warning. Default 0.65.
#' @param ... Additional arguments passed to \code{lavaan::lavaan} or
#'   \code{lavaan::sem}.
#'
#' @return An object of class \code{crossLagR_fit} (a list) with:
#'   \itemize{
#'     \item \code{fit}: the underlying lavaan object.
#'     \item \code{summary}: a \code{crossLagR_summary} (labelled parameters
#'       + fit indices).
#'     \item \code{diagnosis}: a \code{crossLagR_diagnosis} (Heywood cases).
#'     \item \code{icc}: a \code{crossLagR_icc} (NULL if \code{icc_vars} not
#'       supplied).
#'     \item \code{warnings}: character vector of lavaan warnings captured
#'       during fitting.
#'   }
#'   A single condensed \code{warning()} is emitted summarizing Heywood-class
#'   issues and any ICC-threshold flags.
#'
#' @examples
#' \dontrun{
#' syntax <- estimateRICLPM(waves = 4)
#' dat    <- simCLPM(waves = 4, sample_size = 200)$data
#' out    <- lavaanDiagnose(syntax, dat, icc_vars = c("x", "y"))
#' out                         ## prints the full bundled report
#' out$summary                 ## just the labelled estimates + fit
#' out$diagnosis               ## just the Heywood diagnosis
#' out$icc                     ## just the ICC report
#' }
#' @export
lavaanDiagnose <- function(model, data,
                           fit_function = c("lavaan", "sem"),
                           icc_vars = NULL,
                           icc_threshold = 0.65,
                           ...) {
  fit_function <- match.arg(fit_function)
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' is required for lavaanDiagnose().")
  }

  warnings_seen <- character(0)
  fit <- withCallingHandlers({
    if (fit_function == "sem") {
      lavaan::sem(model, data = data, ...)
    } else {
      lavaan::lavaan(model, data = data, ...)
    }
  }, warning = function(w) {
    warnings_seen <<- c(warnings_seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })

  diag <- diagnoseFit(fit, warnings_seen = warnings_seen, quiet = TRUE)
  summ <- tryCatch(crossLagRSummary(fit), error = function(e) NULL)
  icc  <- if (!is.null(icc_vars)) {
    tryCatch(iccReport(data, vars = icc_vars, threshold = icc_threshold,
                       quiet = TRUE),
             error = function(e) NULL)
  } else NULL

  out <- structure(
    list(fit = fit, summary = summ, diagnosis = diag,
         icc = icc, warnings = warnings_seen),
    class = "crossLagR_fit"
  )

  ## Combined condensed warning
  msgs <- character(0)
  if (diag$n_issues > 0L) {
    msgs <- c(msgs, paste0(
      diag$n_issues, " Heywood-class issue(s): ",
      paste(diag$codes, collapse = ", ")
    ))
  }
  if (!is.null(icc) && any(icc$exceeds_threshold)) {
    flagged <- icc[icc$exceeds_threshold, , drop = FALSE]
    msgs <- c(msgs, paste0(
      "ICC > ", icc_threshold, " for ",
      paste(sprintf("%s (%.2f)", flagged$variable, flagged$icc),
            collapse = "; "),
      " -- consider sensitivity check via launch_shiny_sim()"
    ))
  }
  if (length(msgs)) {
    warning("crossLagR: ", paste(msgs, collapse = "; "),
            ". See heywoodHelp() for remedies.", call. = FALSE)
  }
  out
}


#' @export
print.crossLagR_fit <- function(x, ...) {
  if (!is.null(x$summary)) print(x$summary)
  cat("\n")
  if (!is.null(x$diagnosis)) print(x$diagnosis)
  if (!is.null(x$icc)) {
    cat("\n")
    print(x$icc)
  }
  invisible(x)
}

#' @export
summary.crossLagR_fit <- function(object, ...) {
  print(object, ...)
}


#' @title heywoodHelp
#' @description Print the full catalog of Heywood-class issues that
#'   \code{\link{diagnoseFit}} can detect, with plain-language explanations
#'   and remediation suggestions. Useful for browsing the catalog without a
#'   fitted model on hand.
#'
#' @param code Optional character. If supplied, print only the entry for that
#'   code (e.g., \code{"negative_variance"}).
#'
#' @return Invisibly returns the catalog list.
#' @export
heywoodHelp <- function(code = NULL) {
  cat("crossLagR Heywood-case catalog\n",
      strrep("=", 31), "\n\n", sep = "")
  catalog <- if (is.null(code)) {
    .HEYWOOD_CATALOG
  } else {
    if (!code %in% names(.HEYWOOD_CATALOG)) {
      stop("Unknown code '", code, "'. Available: ",
           paste(names(.HEYWOOD_CATALOG), collapse = ", "))
    }
    .HEYWOOD_CATALOG[code]
  }
  for (i in seq_along(catalog)) {
    .print_issue(catalog[[i]], i)
  }
  invisible(catalog)
}


#' @export
print.crossLagR_diagnosis <- function(x, ...) {
  cat(x$summary, "\n", sep = "")
  if (x$n_issues == 0L) return(invisible(x))
  for (i in seq_along(x$issues)) .print_issue(x$issues[[i]], i)
  invisible(x)
}

#' @export
summary.crossLagR_diagnosis <- function(object, ...) {
  cat(object$summary, "\n", sep = "")
  invisible(object)
}


## ----------------------------------------------------------------------------
## Internal: catalog and helpers
## ----------------------------------------------------------------------------

.HEYWOOD_CATALOG <- list(

  non_convergence = list(
    code = "non_convergence",
    severity = "high",
    parameters = character(0),
    message = "Optimizer did not converge.",
    explanation = paste(
      "lavaan's optimizer failed to find a stationary point in the",
      "log-likelihood. This usually reflects an underidentified model,",
      "very poor starting values, or a near-flat likelihood surface."
    ),
    remedies = c(
      "Re-fit with simpler starts: lavaan(..., start = 'simple').",
      "Constrain AR/CL parameters across waves (constrain_beta = TRUE, constrain_omega = TRUE) to reduce empirical underidentification.",
      "Increase the iteration cap: lavaan(..., control = list(iter.max = 5000)).",
      "Estimate in blavaan with weakly-informative priors; MCMC explores the posterior without requiring convergence to a single mode."
    )
  ),

  negative_variance = list(
    code = "negative_variance",
    severity = "high",
    parameters = character(0),
    message = "Negative variance estimate (Heywood case).",
    explanation = paste(
      "lavaan estimated a negative value for a variance parameter.",
      "Variances are non-negative by definition; a negative point estimate",
      "signals the model is over-parameterized relative to the data, or",
      "that the true variance is at the boundary (zero) and sampling noise",
      "has pushed it below."
    ),
    remedies = c(
      "If the negative variance is on a random intercept (I_x, I_y), the data may lack genuine between-person variation in that construct -- consider a plain CLPM (no random intercept) instead of RICLPM.",
      "Inspect the between/within decomposition with withinBetween() -- a near-zero ICC means a within-between model is mis-specified.",
      "Constrain the offending variance to non-negative (e.g., 'I_x ~~ lower(0)*I_x') or fix it to zero.",
      "Estimate in blavaan with a positive-support prior (e.g., gamma(1, 0.5)) on the variance.",
      "Increase sample size and/or number of waves -- boundary estimates are common in small N."
    )
  ),

  std_loading_gt_one = list(
    code = "std_loading_gt_one",
    severity = "medium",
    parameters = character(0),
    message = "Standardized loading exceeds 1.",
    explanation = paste(
      "A standardized factor loading > 1 in absolute value is empirically",
      "impossible under a reflective measurement model -- it implies the",
      "indicator's variance is smaller than the variance inherited from",
      "the latent factor. Usually accompanies a negative residual variance",
      "or a non-PD covariance matrix."
    ),
    remedies = c(
      "Constrain the corresponding residual variance to be non-negative.",
      "Re-check measurement assumptions -- a congeneric or formative spec may be more appropriate.",
      "Estimate in blavaan with priors that bound loadings (e.g., normal(1, 0.3) truncated)."
    )
  ),

  std_correlation_gt_one = list(
    code = "std_correlation_gt_one",
    severity = "high",
    parameters = character(0),
    message = "Latent correlation exceeds 1 (NPD latent covariance).",
    explanation = paste(
      "A correlation > 1 between latent variables signals a non-positive",
      "definite latent covariance matrix. Often arises when two latent",
      "constructs are empirically indistinguishable (collinearity)."
    ),
    remedies = c(
      "Drop one of the redundant latent factors or combine them.",
      "Estimate in blavaan with an LKJ(2) or LKJ(4) prior on the latent correlation matrix to regularize toward identifiability.",
      "Verify that constructs are conceptually distinct -- if not, the model is mis-specified."
    )
  ),

  npd_covariance = list(
    code = "npd_covariance",
    severity = "high",
    parameters = character(0),
    message = "Covariance matrix is not positive definite.",
    explanation = paste(
      "lavaan flagged a non-PD model-implied or sample covariance matrix.",
      "Causes include collinear latent factors, near-singular residual",
      "covariances, or boundary variance estimates."
    ),
    remedies = c(
      "Inspect the model-implied covariance: lavInspect(fit, 'cov.lv').",
      "Constrain redundant residual covariances to zero or equality.",
      "Re-fit with smoothed starting values: lavaan(..., start = 'simple').",
      "Estimate in blavaan with regularizing priors on (co)variances."
    )
  ),

  vcov_not_pd = list(
    code = "vcov_not_pd",
    severity = "medium",
    parameters = character(0),
    message = "Parameter vcov is not PD; SEs are unreliable.",
    explanation = paste(
      "The Hessian at the estimated parameter vector is not positive",
      "definite, so standard errors cannot be trusted. Point estimates",
      "may still be sensible but frequentist inference is compromised."
    ),
    remedies = c(
      "Use bootstrap SEs: lavaan(..., se = 'bootstrap', bootstrap = 1000).",
      "Re-fit from multiple starts; multimodal likelihoods often produce non-PD Hessians.",
      "Switch to blavaan for fully Bayesian uncertainty quantification."
    )
  ),

  iter_limit = list(
    code = "iter_limit",
    severity = "high",
    parameters = character(0),
    message = "Optimizer reached iteration limit without converging.",
    explanation = paste(
      "Either the model is empirically unidentified or the likelihood",
      "surface is very flat near the optimum."
    ),
    remedies = c(
      "Increase iterations: lavaan(..., control = list(iter.max = 5000)).",
      "Constrain parameters across waves (constrain_beta = TRUE, constrain_omega = TRUE).",
      "Reduce model complexity or switch to blavaan."
    )
  )
)

## Inject specific parameter names into a catalog entry; append an "Affected"
## clause to the message so the user sees which parameters are flagged.
.fill_params <- function(iss, params = character(0)) {
  iss$parameters <- params
  if (length(params)) {
    iss$message <- paste0(iss$message, " Affected: ",
                          paste(unique(params), collapse = ", "), ".")
  }
  iss
}

.print_issue <- function(iss, i) {
  cat("[", i, "] ", iss$code, "  (severity: ", iss$severity, ")\n",
      sep = "")
  cat("    ", iss$message, "\n", sep = "")
  cat("    ", iss$explanation, "\n", sep = "")
  if (length(iss$parameters)) {
    cat("    Parameters: ", paste(iss$parameters, collapse = ", "), "\n",
        sep = "")
  }
  cat("    Remedies:\n", sep = "")
  for (r in iss$remedies) cat("      - ", r, "\n", sep = "")
  cat("\n")
}
