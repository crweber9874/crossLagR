#' @title crossLagRSummary
#' @description Produce a tidy, human-readable summary of a lavaan fit using
#'   the crossLagR unified parameter labels (Usami et al., 2019). Each
#'   estimate is annotated with a plain-language description of what the
#'   parameter represents, alongside SE, z-statistic, p-value, and 95\% CI.
#'   Fit indices are returned in a separate compact block.
#'
#' @param fit A fitted lavaan object.
#' @param ci Numeric. Confidence level for intervals. Default 0.95.
#' @param digits Integer. Digits for printed output. Default 3.
#'
#' @return An object of class \code{crossLagR_summary} with elements:
#'   \itemize{
#'     \item \code{parameters}: data.frame of labelled parameters (label,
#'       description, est, se, z, p, ci.lower, ci.upper, std.all).
#'     \item \code{other}: data.frame of remaining (unlabelled) estimated
#'       parameters from the model.
#'     \item \code{fit}: named numeric vector of fit indices (chisq, df,
#'       pvalue, cfi, tli, rmsea, srmr, aic, bic).
#'     \item \code{n_obs}, \code{n_miss}: sample-size bookkeeping (per
#'       CLAUDE.md preference).
#'     \item \code{converged}: logical.
#'   }
#'
#' @details
#' Recognized unified labels (see \code{\link{unifiedLabels}}) include
#' autoregressive/cross-lagged terms (\code{ar_x, ar_y, cl_xy, cl_yx}),
#' random-intercept variances (\code{I_x, I_y, I_z}), latent growth
#' factors (\code{S_x, S_y, Q_x, Q_y}), trait-state decomposition labels
#' (\code{T_var_x, s_var_x}), and dynamic-residual labels
#' (\code{d_var_x, d_var_y, d_cov_xy}). Any label that does not appear in
#' the catalog is shifted to \code{other}.
#'
#' @references
#' Usami, S., Murayama, K., & Hamaker, E. L. (2019). A unified framework of
#' longitudinal models to examine reciprocal relations. \emph{Psychological
#' Methods}, 24(5), 637--657.
#'
#' @examples
#' \dontrun{
#' fit <- lavaan::lavaan(estimateCLPM(waves = 4),
#'                       data = simCLPM(waves = 4, sample_size = 500)$data,
#'                       meanstructure = TRUE)
#' crossLagRSummary(fit)
#' }
#' @export
crossLagRSummary <- function(fit, ci = 0.95, digits = 3) {

  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' is required for crossLagRSummary().")
  }
  if (!inherits(fit, "lavaan")) {
    stop("'fit' must be a fitted lavaan object.")
  }
  if (!is.numeric(ci) || ci <= 0 || ci >= 1) {
    stop("'ci' must be in (0, 1).")
  }

  pe <- lavaan::parameterEstimates(fit, standardized = TRUE,
                                   level = ci, remove.unused = FALSE)

  ## Keep only rows with an estimated label; drop wave-specific suffixes
  ## (e.g., "ar_x2", "ar_x3") down to the base label for description lookup.
  has_label <- nzchar(pe$label %||% "")
  labelled  <- pe[has_label, , drop = FALSE]
  labelled  <- labelled[!duplicated(labelled$label), , drop = FALSE]

  base_lab <- .strip_wave_suffix(labelled$label)
  desc     <- .UNIFIED_LABEL_TABLE$description[match(base_lab,
                                                     .UNIFIED_LABEL_TABLE$label)]
  known    <- !is.na(desc)

  params <- data.frame(
    label       = labelled$label[known],
    description = desc[known],
    est         = labelled$est[known],
    se          = labelled$se[known],
    z           = labelled$z[known],
    p           = labelled$pvalue[known],
    ci.lower    = labelled$ci.lower[known],
    ci.upper    = labelled$ci.upper[known],
    std.all     = labelled$std.all[known],
    stringsAsFactors = FALSE
  )

  other <- data.frame(
    label = labelled$label[!known],
    est   = labelled$est[!known],
    se    = labelled$se[!known],
    z     = labelled$z[!known],
    p     = labelled$pvalue[!known],
    stringsAsFactors = FALSE
  )

  fm <- tryCatch(
    lavaan::fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "tli",
                               "rmsea", "rmsea.ci.lower", "rmsea.ci.upper",
                               "srmr", "aic", "bic")),
    error = function(e) {
      setNames(rep(NA_real_, 11),
               c("chisq","df","pvalue","cfi","tli",
                 "rmsea","rmsea.ci.lower","rmsea.ci.upper",
                 "srmr","aic","bic"))
    }
  )

  n_obs <- tryCatch(as.integer(lavaan::lavInspect(fit, "nobs")),
                    error = function(e) NA_integer_)
  data <- tryCatch(lavaan::lavInspect(fit, "data"),
                   error = function(e) NULL)
  n_miss <- if (is.matrix(data) || is.data.frame(data)) {
    sum(!stats::complete.cases(data))
  } else NA_integer_

  out <- list(
    parameters = params,
    other      = other,
    fit        = as.numeric(fm),
    fit_names  = names(fm),
    n_obs      = n_obs,
    n_miss     = n_miss,
    converged  = isTRUE(lavaan::lavInspect(fit, "converged")),
    ci_level   = ci,
    digits     = digits
  )
  names(out$fit) <- out$fit_names
  out$fit_names <- NULL
  class(out) <- "crossLagR_summary"
  out
}


#' @title unifiedLabels
#' @description Print or return the catalog of unified parameter labels used
#'   across crossLagR estimators, with plain-language descriptions.
#'
#' @param print Logical. If \code{TRUE} (default), pretty-prints the table.
#'   If \code{FALSE}, returns it silently.
#'
#' @return A data.frame with columns \code{label}, \code{description},
#'   \code{group}.
#' @export
unifiedLabels <- function(print = TRUE) {
  if (print) {
    cat("crossLagR unified parameter labels\n",
        strrep("=", 35), "\n", sep = "")
    by_group <- split(.UNIFIED_LABEL_TABLE, .UNIFIED_LABEL_TABLE$group)
    for (g in names(by_group)) {
      cat("\n[", g, "]\n", sep = "")
      sub <- by_group[[g]]
      width <- max(nchar(sub$label)) + 2L
      for (i in seq_len(nrow(sub))) {
        cat("  ", format(sub$label[i], width = width),
            "  ", sub$description[i], "\n", sep = "")
      }
    }
  }
  invisible(.UNIFIED_LABEL_TABLE)
}


#' @export
print.crossLagR_summary <- function(x, ...) {
  digits <- x$digits %||% 3L

  cat("crossLagR model summary",
      "  (CI: ", x$ci_level * 100, "%)\n", sep = "")
  cat(strrep("=", 50), "\n", sep = "")
  cat("Converged: ", isTRUE(x$converged),
      "    N obs: ", x$n_obs %||% NA,
      "    N missing: ", x$n_miss %||% NA, "\n", sep = "")
  cat("\n", sep = "")

  ## --- parameter table ----------------------------------------------------
  if (nrow(x$parameters) == 0L) {
    cat("(no unified-label parameters were estimated)\n")
  } else {
    cat("Parameters (unified labels)\n")
    cat(strrep("-", 50), "\n", sep = "")
    p <- x$parameters
    fp <- function(v) formatC(v, digits = digits, format = "f")
    sig <- vapply(p$p, .sig_stars, character(1))
    tab <- data.frame(
      label       = p$label,
      description = .truncate(p$description, 40),
      est         = fp(p$est),
      se          = fp(p$se),
      z           = fp(p$z),
      p           = fp(p$p),
      `95% CI`    = paste0("[", fp(p$ci.lower), ", ", fp(p$ci.upper), "]"),
      sig         = sig,
      std         = fp(p$std.all),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    print(tab, row.names = FALSE, right = FALSE)
    cat("---\nSignificance: '***' p<.001  '**' p<.01  '*' p<.05  '.' p<.10\n")
  }

  ## --- other (unlabelled) -------------------------------------------------
  if (nrow(x$other) > 0L) {
    cat("\nOther estimated parameters (no unified label):\n")
    cat(strrep("-", 50), "\n", sep = "")
    o <- x$other
    fp <- function(v) formatC(v, digits = digits, format = "f")
    print(data.frame(
      label = o$label, est = fp(o$est), se = fp(o$se),
      z = fp(o$z), p = fp(o$p),
      stringsAsFactors = FALSE
    ), row.names = FALSE, right = FALSE)
  }

  ## --- fit indices --------------------------------------------------------
  cat("\nFit statistics\n")
  cat(strrep("-", 50), "\n", sep = "")
  f <- x$fit
  fp <- function(v) formatC(v, digits = digits, format = "f")
  cat(sprintf("  chi-square = %s   df = %s   p = %s\n",
              fp(f["chisq"]), formatC(f["df"], format = "d"), fp(f["pvalue"])))
  cat(sprintf("  CFI = %s   TLI = %s\n", fp(f["cfi"]), fp(f["tli"])))
  cat(sprintf("  RMSEA = %s   [%s, %s]\n",
              fp(f["rmsea"]), fp(f["rmsea.ci.lower"]),
              fp(f["rmsea.ci.upper"])))
  cat(sprintf("  SRMR = %s\n", fp(f["srmr"])))
  cat(sprintf("  AIC = %s   BIC = %s\n", fp(f["aic"]), fp(f["bic"])))

  invisible(x)
}

#' @export
summary.crossLagR_summary <- function(object, ...) {
  print(object, ...)
}


## ----------------------------------------------------------------------------
## Helpers and catalog
## ----------------------------------------------------------------------------

.strip_wave_suffix <- function(labels) {
  ## ar_x2 -> ar_x; cl_xy3 -> cl_xy; I_x stays I_x.
  sub("([a-zA-Z_]+?)([0-9]+)$", "\\1", labels)
}

.sig_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < .001) return("***")
  if (p < .01)  return("**")
  if (p < .05)  return("*")
  if (p < .10)  return(".")
  ""
}

.truncate <- function(x, n) {
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 1), "..."), x)
}

## Source-of-truth catalog. Each row maps a unified label to a plain-language
## description and a coarse group used for grouped printing in unifiedLabels().
.UNIFIED_LABEL_TABLE <- (function() {
  rows <- list(
    ## Dynamic (autoregressive + cross-lagged) per Usami et al. (2019)
    list("ar_x",     "Autoregressive effect of X (X(t-1) -> X(t))",        "dynamics"),
    list("ar_y",     "Autoregressive effect of Y (Y(t-1) -> Y(t))",        "dynamics"),
    list("ar_z",     "Autoregressive effect of Z (Z(t-1) -> Z(t))",        "dynamics"),
    list("cl_xy",    "Cross-lagged effect of X on Y (X(t-1) -> Y(t))",     "dynamics"),
    list("cl_yx",    "Cross-lagged effect of Y on X (Y(t-1) -> X(t))",     "dynamics"),
    list("cl_xz",    "Cross-lagged effect of X on Z (X(t-1) -> Z(t))",     "dynamics"),
    list("cl_yz",    "Cross-lagged effect of Y on Z (Y(t-1) -> Z(t))",     "dynamics"),
    list("cl_zx",    "Cross-lagged effect of Z on X (Z(t-1) -> X(t))",     "dynamics"),
    list("cl_zy",    "Cross-lagged effect of Z on Y (Z(t-1) -> Y(t))",     "dynamics"),

    ## Dynamic residuals
    list("d_var_x",  "Within-person residual variance of X",               "dynamic_residuals"),
    list("d_var_y",  "Within-person residual variance of Y",               "dynamic_residuals"),
    list("d_cov_xy", "Within-person residual covariance X<->Y",            "dynamic_residuals"),

    ## Random intercepts / between-person factors
    list("I_x",      "Random intercept for X (stable between-person mean)", "random_intercepts"),
    list("I_y",      "Random intercept for Y (stable between-person mean)", "random_intercepts"),
    list("I_z",      "Random intercept for Z (stable between-person mean)", "random_intercepts"),

    ## Latent growth factors
    list("S_x",      "Slope factor for X (rate of linear change)",         "growth"),
    list("S_y",      "Slope factor for Y (rate of linear change)",         "growth"),
    list("Q_x",      "Quadratic factor for X (acceleration)",              "growth"),
    list("Q_y",      "Quadratic factor for Y (acceleration)",              "growth"),
    list("mean_I",   "Mean intercept (univariate growth)",                 "growth"),
    list("mean_S",   "Mean slope (univariate growth)",                     "growth"),
    list("mean_Q",   "Mean quadratic (univariate growth)",                 "growth"),
    list("mean_I_x", "Mean intercept of X",                                "growth"),
    list("mean_I_y", "Mean intercept of Y",                                "growth"),
    list("mean_S_x", "Mean slope of X",                                    "growth"),
    list("mean_S_y", "Mean slope of Y",                                    "growth"),
    list("var_I",    "Intercept variance (between-person)",                "growth"),
    list("var_S",    "Slope variance (between-person rate-of-change)",     "growth"),
    list("var_Q",    "Quadratic variance",                                 "growth"),
    list("var_I_x",  "Intercept variance of X",                            "growth"),
    list("var_I_y",  "Intercept variance of Y",                            "growth"),
    list("var_S_x",  "Slope variance of X",                                "growth"),
    list("var_S_y",  "Slope variance of Y",                                "growth"),
    list("cov_IS",   "Intercept-slope covariance",                         "growth"),
    list("cov_IS_x", "Intercept-slope covariance for X",                   "growth"),
    list("cov_IS_y", "Intercept-slope covariance for Y",                   "growth"),

    ## Trait-state-occasion decomposition
    list("T_var_x",  "Trait variance of X (between-person)",               "trait_state"),
    list("T_var_y",  "Trait variance of Y (between-person)",               "trait_state"),
    list("T_cov_xy", "Trait covariance X<->Y",                             "trait_state"),
    list("s_var_x",  "State variance of X (within-person)",                "trait_state"),
    list("s_var_y",  "State variance of Y (within-person)",                "trait_state"),
    list("s_cov_xy", "Within-time state covariance X<->Y",                 "trait_state"),

    ## Univariate (LChange-style) residual
    list("u_var_x",  "Wave-specific residual variance of X",               "residuals"),
    list("u_var_y",  "Wave-specific residual variance of Y",               "residuals"),

    ## Constant-change / accumulating factor (LChange)
    list("A_x",      "Accumulating-change factor for X",                   "change"),
    list("A_y",      "Accumulating-change factor for Y",                   "change")
  )
  out <- do.call(rbind, lapply(rows, function(r) {
    data.frame(label = r[[1]], description = r[[2]], group = r[[3]],
               stringsAsFactors = FALSE)
  }))
  out
})()
