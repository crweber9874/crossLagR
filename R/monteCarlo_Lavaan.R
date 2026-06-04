#' @title monteCarloLavaan
#' @description Generic Monte Carlo wrapper for any lavaan-syntax estimator in the
#'   crossLagR unified framework (CLPM, RICLPM, ALT, LGM, LCMSR, BollenBrand, TSO,
#'   LChange). Generates data once per trial under a chosen DGP, fits the requested
#'   model, and returns unified-label parameter estimates alongside fit indices.
#'
#' @param estimator Character. One of "CLPM", "RICLPM", "ALT", "LGM", "LCMSR",
#'   "BB", "TSO", "LCHANGE".
#' @param dgp Character. One of "clpm", "riclpm", "clpmu", "lchange".
#' @param trials Integer. Trials per parameter cell.
#' @param waves Integer. Number of waves.
#' @param sample_size Integer. Sample size per simulated dataset.
#' @param param_grid data.frame with one row of DGP parameters (the per-cell row
#'   that \code{run_mc_sims} feeds in). Required columns: stability_p, stability_q,
#'   cross_p, cross_q, variance_p, variance_q, variance_between_x, variance_between_y,
#'   cov_pq. Optional confounder_* columns for dgp = "clpmu".
#' @param estimator_args Named list. Extra arguments passed to the estimator's
#'   syntax-generating function (e.g., \code{constrain_beta = TRUE},
#'   \code{variable_type = "bivariate"}).
#' @param fit_function Character. Either "lavaan" (default; preserves user
#'   constraints exactly) or "sem" (adds defaults like fixed.x = FALSE).
#' @param verbose Logical. Print per-trial progress and errors. Default FALSE.
#'
#' @details
#' Returns one row per trial containing:
#' \itemize{
#'   \item Unified-label point estimates (where present): ar_x, ar_y, cl_xy, cl_yx.
#'   \item Standard errors: se_ar_x, se_ar_y, se_cl_xy, se_cl_yx.
#'   \item Convergence and fit: converged, n_obs, n_miss, chisq, df, pvalue,
#'     cfi, tli, rmsea, srmr, aic, bic.
#'   \item True parameter values from the DGP for direct bias computation:
#'     true_ar_x, true_ar_y, true_cl_xy, true_cl_yx.
#'   \item Bookkeeping: trial, estimator, dgp, error_occurred, error_message.
#' }
#'
#' Per CLAUDE.md: n_obs and n_miss are always reported alongside the estimates.
#'
#' @return data.frame with one row per trial.
#'
#' @export
monteCarloLavaan <- function(estimator,
                             dgp = "clpm",
                             trials = 10,
                             waves = 5,
                             sample_size = 1000,
                             param_grid = NULL,
                             estimator_args = list(),
                             fit_function = c("lavaan", "sem"),
                             verbose = FALSE) {

  fit_function <- match.arg(fit_function)

  valid_estimators <- c("CLPM", "RICLPM", "ALT", "LGM", "LCMSR", "BB", "TSO", "LCHANGE")
  if (!estimator %in% valid_estimators) {
    stop("estimator must be one of: ", paste(valid_estimators, collapse = ", "))
  }

  ## DGPs:
  ##   "clpmu"  → simCLPMu (the only DGP with an unmeasured confounder).
  ##   any unified-label estimator name (CLPM, RICLPM, ALT, LGM, LCMSR,
  ##     LCHANGE, BB, TSO) → simFromSyntax, populating the estimator's own
  ##     lavaan syntax with the cell's true ar/cl/variance values.
  ##
  ## Legacy lowercase codes ("clpm","riclpm","lchange") plus the historic
  ## "_syn" suffixed forms ("CLPM_syn", etc.) are silently normalized to the
  ## canonical estimator-name form below so older callers keep working.
  dgp <- switch(
    dgp,
    "clpm"        = "CLPM",
    "riclpm"      = "RICLPM",
    "lchange"     = "LCHANGE",
    "CLPM_syn"    = "CLPM",
    "RICLPM_syn"  = "RICLPM",
    "LCHANGE_syn" = "LCHANGE",
    dgp
  )
  valid_dgps <- c("CLPM", "RICLPM", "ALT", "LGM", "LCMSR",
                  "LCHANGE", "BB", "TSO", "clpmu")
  if (!dgp %in% valid_dgps) {
    stop("dgp must be one of: ", paste(valid_dgps, collapse = ", "))
  }

  if (is.null(param_grid) || nrow(param_grid) == 0) {
    param_grid <- data.frame(
      stability_p = 0.3, stability_q = 0.3,
      cross_p = 0.1, cross_q = 0.1,
      variance_p = 0.5, variance_q = 0.5,
      variance_between_x = 0.5, variance_between_y = 0.5,
      cov_pq = 0
    )
  }

  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' is required for monteCarloLavaan().")
  }

  ## Build the lavaan model syntax once (depends only on waves + estimator_args)
  syntax <- build_estimator_syntax(estimator, waves, estimator_args)

  results <- vector("list", trials * nrow(param_grid))
  k <- 0L

  for (i in seq_len(nrow(param_grid))) {
    params <- as.list(param_grid[i, , drop = FALSE])

    ## True parameter values (use stability_p/q and cross_p/q as the DGP truths;
    ## these map to unified labels at the within-person level).
    true_vals <- list(
      true_ar_x  = params$stability_p %||% NA_real_,
      true_ar_y  = params$stability_q %||% NA_real_,
      ## Naming convention: cross_p = effect of Y on X (cl_yx);
      ## cross_q = effect of X on Y (cl_xy). Matches simCLPM/simRICLPM.
      true_cl_yx = params$cross_p   %||% NA_real_,
      true_cl_xy = params$cross_q   %||% NA_real_
    )

    for (j in seq_len(trials)) {
      k <- k + 1L

      trial_result <- tryCatch({
        dat <- simulate_dgp(dgp, waves, sample_size, params)
        res <- fit_one(syntax, dat, fit_function)
        ## Attach per-trial ICCs computed directly on the simulated data
        ## (cheap, runs in O(N*T)). Always uses 'x' and 'y' stems since
        ## every simFromSyntax / simCLPM data path produces those columns.
        icc <- tryCatch(
          iccReport(dat, vars = c("x", "y"), threshold = 0.65, quiet = TRUE),
          error = function(e) NULL
        )
        if (!is.null(icc)) {
          res$icc_x         <- icc$icc[icc$variable == "x"][1]
          res$icc_y         <- icc$icc[icc$variable == "y"][1]
          res$icc_exceeds   <- any(icc$exceeds_threshold)
        }
        res
      }, error = function(e) {
        if (verbose) message("[", estimator, "] cell ", i, " trial ", j, ": ", e$message)
        list(error_occurred = TRUE, error_message = as.character(e$message))
      })

      row <- merge_row(
        trial_result,
        list(trial = j, param_combo = i, estimator = estimator, dgp = dgp),
        params, true_vals
      )

      results[[k]] <- row
    }
  }

  out <- dplyr::bind_rows(results)

  ## Single condensed summary if any trials hit Heywood-class issues.
  ## We avoid per-trial warnings because an MC run easily produces hundreds.
  n_hit <- sum(out$heywood_n_issues > 0L, na.rm = TRUE)
  if (n_hit > 0L) {
    codes_seen <- unique(unlist(strsplit(
      out$heywood_codes[!is.na(out$heywood_codes)], ";"
    )))
    warning(
      "crossLagR: ", n_hit, " of ", nrow(out),
      " trials triggered Heywood-class issues (",
      paste(codes_seen, collapse = ", "),
      "). See `heywood_codes` column or call heywoodHelp() for remedies.",
      call. = FALSE
    )
  }
  n_icc <- sum(isTRUE(out$icc_exceeds) | out$icc_exceeds %in% TRUE, na.rm = TRUE)
  if (n_icc > 0L) {
    warning(
      "crossLagR: ", n_icc, " of ", nrow(out),
      " trials had ICC > 0.65 on x or y. Consider exploring sensitivity ",
      "with launch_shiny_sim().",
      call. = FALSE
    )
  }
  out
}

# ---- helpers -----------------------------------------------------------------

`%||%` <- function(a, b) if (!is.null(a)) a else b

## Dispatch table: estimator -> function generating lavaan syntax
build_estimator_syntax <- function(estimator, waves, args = list()) {
  args <- as.list(args)
  args$waves <- waves

  fn <- switch(
    estimator,
    "CLPM"    = crossLagR::estimateCLPM,
    "RICLPM"  = crossLagR::estimateRICLPM,
    "ALT"     = crossLagR::estimateALT,
    "LGM"     = crossLagR::estimateLGM,
    "LCMSR"   = crossLagR::estimateLCMSR,
    "BB"      = crossLagR::estimateBollen_and_Brand,
    "TSO"     = crossLagR::estimateTSO,
    "LCHANGE" = crossLagR::estimateLChange,
    stop("Unknown estimator: ", estimator)
  )

  ## Trim args to those accepted by the target function
  formal_names <- names(formals(fn))
  args <- args[intersect(names(args), formal_names)]

  do.call(fn, args)
}

## Dispatch DGP -> data. After the normalization in monteCarloLavaan() the
## dgp here is one of: "clpmu" (parameterized confounder sim) or any
## unified-label estimator name ("CLPM","RICLPM","ALT","LGM","LCMSR",
## "LCHANGE","BB","TSO"), which route through simFromSyntax.
simulate_dgp <- function(dgp, waves, sample_size, params) {
  if (dgp == "clpmu") {
    sim <- crossLagR::simCLPMu(
      waves                = waves,
      stability_p          = params$stability_p          %||% 0.3,
      stability_q          = params$stability_q          %||% 0.3,
      cross_p              = params$cross_p              %||% 0,
      cross_q              = params$cross_q              %||% 0,
      variance_p           = params$variance_p           %||% 0.5,
      variance_q           = params$variance_q           %||% 0.5,
      cov_pq               = params$cov_pq               %||% 0,
      confounder_p         = params$confounder_p         %||% 0.3,
      confounder_q         = params$confounder_q         %||% 0.3,
      confounder_variance  = params$confounder_variance  %||% 1,
      confounder_stability = params$confounder_stability %||% 0.4,
      sample.nobs          = sample_size
    )
    return(sim$data)
  }

  ## Everything else: simulate from the estimator's populated lavaan syntax.
  if (dgp %in% c("CLPM", "RICLPM", "ALT", "LGM", "LCMSR",
                 "LCHANGE", "BB", "TSO")) {
    return(crossLagR::simFromSyntax(
      estimator   = dgp,
      waves       = waves,
      sample_size = sample_size,
      ar_x        = params$stability_p %||% 0.3,
      ar_y        = params$stability_q %||% 0.3,
      cl_xy       = params$cross_q     %||% 0.1,
      cl_yx       = params$cross_p     %||% 0.1,
      d_var_x     = params$variance_p  %||% 0.5,
      d_var_y     = params$variance_q  %||% 0.5,
      d_cov_xy    = params$cov_pq      %||% 0
    ))
  }

  stop("Unsupported dgp: ", dgp)
}

## Fit and extract unified-label params + fit indices.
##
## Caps the optimizer at 200 iterations: lavaan's nlminb default is 10000,
## which means a pathological cell (e.g. high ICC with mis-specified CLPM) can
## spin for minutes before giving up. 200 is enough for a well-posed model and
## marks the bad ones non-converged quickly. Also captures warnings so the
## caller can see how many fits had Heywood cases / non-positive-definite
## covariance matrices.
fit_one <- function(syntax, data, fit_function = "lavaan") {
  warnings_seen <- character(0)

  fit <- withCallingHandlers({
    if (fit_function == "sem") {
      lavaan::sem(syntax, data = data, meanstructure = TRUE,
                  control = list(iter.max = 500L))
    } else {
      lavaan::lavaan(syntax, data = data, meanstructure = TRUE,
                     control = list(iter.max = 500L))
    }
  }, warning = function(w) {
    warnings_seen <<- c(warnings_seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })

  ## Categorize warnings so the worker can summarize them per cell.
  ## We look for the most common lavaan pathologies.
  has_warn <- function(pattern) any(grepl(pattern, warnings_seen, ignore.case = TRUE))
  warn_npd     <- has_warn("not positive definite|positive[- ]definite")
  warn_heywood <- has_warn("variances are negative|heywood")
  warn_iter    <- has_warn("did not converge|reached the iteration limit")

  converged <- isTRUE(lavaan::lavInspect(fit, "converged"))
  if (!converged) {
    diag <- tryCatch(
      diagnoseFit(fit, warnings_seen = warnings_seen, quiet = TRUE),
      error = function(e) NULL
    )
    return(list(
      error_occurred    = TRUE,
      error_message     = if (warn_iter) "iter_limit" else "non_convergence",
      converged         = FALSE,
      warning_count     = length(warnings_seen),
      warn_npd          = warn_npd,
      warn_heywood      = warn_heywood,
      warn_iter         = warn_iter,
      heywood_n_issues  = if (is.null(diag)) 0L else diag$n_issues,
      heywood_codes     = if (is.null(diag) || diag$n_issues == 0L) NA_character_
                          else paste(diag$codes, collapse = ";")
    ))
  }

  pe <- lavaan::parameterEstimates(fit, standardized = TRUE)

  pick <- function(label, col = "est") {
    rows <- which(pe$label == label)
    if (!length(rows)) return(NA_real_)
    suppressWarnings(as.numeric(pe[[col]][rows[1]]))
  }

  fm <- tryCatch(
    lavaan::fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "tli",
                               "rmsea", "srmr", "aic", "bic")),
    error = function(e) setNames(rep(NA_real_, 9),
                                 c("chisq","df","pvalue","cfi","tli",
                                   "rmsea","srmr","aic","bic"))
  )

  ## n_obs / n_miss per CLAUDE.md preference
  n_obs <- tryCatch(lavaan::lavInspect(fit, "nobs"),
                    error = function(e) nrow(data))
  n_miss <- if (is.data.frame(data)) sum(!stats::complete.cases(data)) else NA_integer_

  diag <- tryCatch(
    diagnoseFit(fit, warnings_seen = warnings_seen, quiet = TRUE),
    error = function(e) NULL
  )
  heywood_n     <- if (is.null(diag)) 0L else diag$n_issues
  heywood_codes <- if (is.null(diag) || heywood_n == 0L) NA_character_
                   else paste(diag$codes, collapse = ";")

  list(
    error_occurred   = FALSE,
    error_message    = NA_character_,
    converged        = TRUE,
    n_obs            = as.integer(n_obs),
    n_miss           = as.integer(n_miss),
    warning_count    = length(warnings_seen),
    warn_npd         = warn_npd,
    warn_heywood     = warn_heywood,
    warn_iter        = warn_iter,
    heywood_n_issues = heywood_n,
    heywood_codes    = heywood_codes,

    ar_x = pick("ar_x"),  se_ar_x = pick("ar_x", "se"),  std_ar_x = pick("ar_x", "std.all"),
    ar_y = pick("ar_y"),  se_ar_y = pick("ar_y", "se"),  std_ar_y = pick("ar_y", "std.all"),
    cl_xy = pick("cl_xy"), se_cl_xy = pick("cl_xy", "se"), std_cl_xy = pick("cl_xy", "std.all"),
    cl_yx = pick("cl_yx"), se_cl_yx = pick("cl_yx", "se"), std_cl_yx = pick("cl_yx", "std.all"),

    chisq  = unname(fm["chisq"]),
    df     = unname(fm["df"]),
    pvalue = unname(fm["pvalue"]),
    cfi    = unname(fm["cfi"]),
    tli    = unname(fm["tli"]),
    rmsea  = unname(fm["rmsea"]),
    srmr   = unname(fm["srmr"]),
    aic    = unname(fm["aic"]),
    bic    = unname(fm["bic"])
  )
}

## Merge a trial result with bookkeeping/parameter columns into a 1-row data.frame
merge_row <- function(trial_result, meta, params, true_vals) {
  base <- list(
    trial         = meta$trial,
    param_combo   = meta$param_combo,
    estimator     = meta$estimator,
    dgp           = meta$dgp,
    error_occurred = isTRUE(trial_result$error_occurred),
    error_message  = trial_result$error_message %||% NA_character_,
    converged     = trial_result$converged %||% FALSE,
    n_obs         = trial_result$n_obs   %||% NA_integer_,
    n_miss        = trial_result$n_miss  %||% NA_integer_,
    warning_count = trial_result$warning_count %||% 0L,
    warn_npd      = trial_result$warn_npd      %||% FALSE,
    warn_heywood  = trial_result$warn_heywood  %||% FALSE,
    warn_iter     = trial_result$warn_iter     %||% FALSE,
    heywood_n_issues = trial_result$heywood_n_issues %||% 0L,
    heywood_codes    = trial_result$heywood_codes    %||% NA_character_,
    icc_x            = trial_result$icc_x            %||% NA_real_,
    icc_y            = trial_result$icc_y            %||% NA_real_,
    icc_exceeds      = trial_result$icc_exceeds      %||% NA,

    ar_x  = trial_result$ar_x  %||% NA_real_,
    ar_y  = trial_result$ar_y  %||% NA_real_,
    cl_xy = trial_result$cl_xy %||% NA_real_,
    cl_yx = trial_result$cl_yx %||% NA_real_,

    se_ar_x  = trial_result$se_ar_x  %||% NA_real_,
    se_ar_y  = trial_result$se_ar_y  %||% NA_real_,
    se_cl_xy = trial_result$se_cl_xy %||% NA_real_,
    se_cl_yx = trial_result$se_cl_yx %||% NA_real_,

    std_ar_x  = trial_result$std_ar_x  %||% NA_real_,
    std_ar_y  = trial_result$std_ar_y  %||% NA_real_,
    std_cl_xy = trial_result$std_cl_xy %||% NA_real_,
    std_cl_yx = trial_result$std_cl_yx %||% NA_real_,

    chisq  = trial_result$chisq  %||% NA_real_,
    df     = trial_result$df     %||% NA_real_,
    pvalue = trial_result$pvalue %||% NA_real_,
    cfi    = trial_result$cfi    %||% NA_real_,
    tli    = trial_result$tli    %||% NA_real_,
    rmsea  = trial_result$rmsea  %||% NA_real_,
    srmr   = trial_result$srmr   %||% NA_real_,
    aic    = trial_result$aic    %||% NA_real_,
    bic    = trial_result$bic    %||% NA_real_
  )

  ## Append DGP parameter columns and true values for bias computation
  base <- c(base, true_vals,
            params[setdiff(names(params), names(base))])

  as.data.frame(base, stringsAsFactors = FALSE)
}
