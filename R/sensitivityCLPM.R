#' @title sensitivityCLPM
#' @description Diagnose how badly a fitted CLPM is biased by between-person
#'   trait variance. Fits an RI-CLPM on the same data to obtain the within-
#'   person truth, computes the data's observed ICC, then simulates CLPM
#'   estimates across a grid of hypothetical ICCs to show the bias curve
#'   relative to that truth. The user's actual CLPM estimate is plotted at
#'   the observed ICC so the bias is immediate and concrete.
#'
#' @param fit A fitted lavaan CLPM object (from
#'   \code{lavaan::lavaan}/\code{sem}) whose syntax was produced by
#'   \code{\link{estimateCLPM}}.
#' @param icc_grid Numeric vector of hypothetical ICC values for the
#'   simulation sweep. Default \code{seq(0, 0.8, by = 0.1)}; the function
#'   adds the observed ICC to the grid automatically.
#' @param n_sims Integer. Replications per ICC level. Default 30 (enough to
#'   see the curve; bump to 100+ for publication-quality SE).
#' @param sample_size Integer. Sample size per simulated dataset. Defaults
#'   to the original fit's \code{lavInspect(fit, "nobs")}.
#' @param waves Integer. Inferred from the fit by default.
#' @param data Optional data.frame to use for the RI-CLPM refit and ICC
#'   computation. If omitted, extracted from the lavaan fit via
#'   \code{lavInspect(fit, "data")}.
#' @param seed Optional integer for reproducibility.
#' @param verbose Logical. Print per-cell progress. Default \code{FALSE}.
#'
#' @return An object of class \code{crossLagR_sensitivity} with elements:
#'   \itemize{
#'     \item \code{estimates}: long data.frame of CLPM AR/CL estimates from
#'       simulation across \code{icc_grid}.
#'     \item \code{clpm_actual}: named numeric vector of the user's actual
#'       CLPM estimates (the biased point of departure).
#'     \item \code{riclpm_truth}: named numeric vector of the RI-CLPM
#'       estimates on the same data --- treated as the within-person
#'       truth (the dashed baseline in the plot).
#'     \item \code{observed_icc}: named numeric of empirical ICCs for x and
#'       y in the data.
#'     \item \code{implied_bias}: named numeric of
#'       \code{clpm_actual - riclpm_truth} --- the estimated bias.
#'   }
#'   Has \code{print()} and \code{plot()} methods.
#'
#' @details
#' \strong{The framing.} The CLPM is known to be biased when the data carry
#' a stable trait component (between-person ICC > 0): lagged predictors
#' correlate with the residual because both inherit the trait
#' (@hamaker2015critique). This function answers the natural follow-up:
#' \emph{by how much is my CLPM biased?}
#'
#' The plot shows three pieces of information:
#' \enumerate{
#'   \item \strong{Dashed horizontal line} = RI-CLPM estimate on the same
#'     data, treated as the within-person truth.
#'   \item \strong{Red curve} = mean CLPM estimate across simulated DGPs
#'     with the RI-CLPM's within-person dynamics but varying ICC --- the
#'     bias trajectory.
#'   \item \strong{Red diamond at the observed ICC} = the user's actual
#'     CLPM estimate. The vertical gap between the diamond and the dashed
#'     line is the bias on \emph{your} data.
#' }
#'
#' \strong{What this is not.} It is not a magical bias correction. The
#' RI-CLPM is treated as truth, which is the standard move but rests on
#' the RI-CLPM's own assumptions (uncorrelated trait/wave-1 within
#' factors, no time-varying confounders). If those are violated --- e.g.
#' the data is Bollen-Brand-shaped --- the RI-CLPM is also biased and the
#' "truth" reference shifts. Triangulate with other estimators
#' (Bollen & Brand, LCM-SR) when the stakes are high.
#'
#' @references
#' Hamaker, E. L., Kuiper, R. M., & Grasman, R. P. P. P. (2015). A critique
#' of the cross-lagged panel model. \emph{Psychological Methods}, 20(1),
#' 102--116.
#'
#' @examples
#' \dontrun{
#' dat <- simRICLPM(waves = 4, sample.nobs = 500,
#'                  var_BX = 1.5, var_BY = 1.5)$data
#' fit <- lavaan::lavaan(estimateCLPM(waves = 4), data = dat,
#'                       meanstructure = TRUE)
#' s <- sensitivityCLPM(fit, icc_grid = seq(0, 0.7, 0.1), n_sims = 50)
#' s          ## prints CLPM-actual, RI-CLPM-truth, observed ICC, bias
#' plot(s)    ## bias curve with observed-ICC marker
#' }
#' @export
sensitivityCLPM <- function(fit,
                            icc_grid = seq(0, 0.8, by = 0.1),
                            n_sims = 30,
                            sample_size = NULL,
                            waves = NULL,
                            data = NULL,
                            seed = NULL,
                            verbose = FALSE) {

  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' is required for sensitivityCLPM().")
  }
  if (!inherits(fit, "lavaan")) {
    stop("'fit' must be a fitted lavaan object.")
  }
  if (!is.numeric(icc_grid) || any(icc_grid < 0 | icc_grid >= 1)) {
    stop("'icc_grid' values must be in [0, 1).")
  }
  if (!isTRUE(lavaan::lavInspect(fit, "converged"))) {
    stop("The supplied lavaan fit did not converge; sensitivity is undefined.")
  }
  if (!is.null(seed)) set.seed(seed)

  ## ---- 1. user's actual (biased) CLPM estimates ---------------------------
  pe <- lavaan::parameterEstimates(fit)
  pick <- function(p, lab) {
    v <- p$est[p$label == lab]
    if (!length(v)) NA_real_ else as.numeric(v[1])
  }
  clpm_actual <- c(
    ar_x  = pick(pe, "ar_x"),  ar_y  = pick(pe, "ar_y"),
    cl_xy = pick(pe, "cl_xy"), cl_yx = pick(pe, "cl_yx")
  )
  if (any(is.na(clpm_actual))) {
    stop("Could not extract unified CLPM labels from fit. ",
         "Was the syntax built with estimateCLPM()?")
  }

  ## ---- 2. waves / N -------------------------------------------------------
  if (is.null(sample_size)) {
    sample_size <- tryCatch(as.integer(lavaan::lavInspect(fit, "nobs")),
                            error = function(e) 500L)
  }
  if (is.null(waves)) {
    obs <- tryCatch(lavaan::lavNames(fit, "ov"), error = function(e) character(0))
    waves <- max(c(0L, suppressWarnings(as.integer(sub("^[xy]", "", obs)))),
                 na.rm = TRUE)
    if (!is.finite(waves) || waves < 3L) {
      stop("Could not infer waves from fit; supply 'waves' explicitly.")
    }
  }

  ## ---- 3. RI-CLPM truth ---------------------------------------------------
  if (is.null(data)) {
    data <- tryCatch(as.data.frame(lavaan::lavInspect(fit, "data")),
                     error = function(e) NULL)
    if (is.null(data)) {
      stop("Could not extract data from fit; supply 'data' explicitly.")
    }
    obs_names <- lavaan::lavNames(fit, "ov")
    if (ncol(data) == length(obs_names)) names(data) <- obs_names
  }

  riclpm_fit <- tryCatch(
    suppressWarnings(lavaan::lavaan(
      crossLagR::estimateRICLPM(waves = waves),
      data = data, meanstructure = TRUE,
      control = list(iter.max = 500L))),
    error = function(e) NULL
  )
  if (is.null(riclpm_fit) ||
      !isTRUE(lavaan::lavInspect(riclpm_fit, "converged"))) {
    stop("RI-CLPM refit did not converge on the supplied data; cannot ",
         "construct a within-person truth reference. Consider ",
         "iccReport(data) and inspecting the ICC structure manually.")
  }
  pe_ri <- lavaan::parameterEstimates(riclpm_fit)
  riclpm_truth <- c(
    ar_x  = pick(pe_ri, "ar_x"),  ar_y  = pick(pe_ri, "ar_y"),
    cl_xy = pick(pe_ri, "cl_xy"), cl_yx = pick(pe_ri, "cl_yx")
  )
  d_var_x  <- pick(pe_ri, "d_var_x")  %||% pick(pe, "d_var_x")
  d_var_y  <- pick(pe_ri, "d_var_y")  %||% pick(pe, "d_var_y")
  d_cov_xy <- pick(pe_ri, "d_cov_xy") %||% pick(pe, "d_cov_xy") %||% 0

  ## ---- 4. observed ICC on the actual data ---------------------------------
  icc_obs_df <- tryCatch(
    crossLagR::iccReport(data, vars = c("x", "y"), quiet = TRUE),
    error = function(e) NULL
  )
  observed_icc <- if (is.null(icc_obs_df)) {
    c(x = NA_real_, y = NA_real_)
  } else {
    c(x = icc_obs_df$icc[icc_obs_df$variable == "x"][1],
      y = icc_obs_df$icc[icc_obs_df$variable == "y"][1])
  }

  ## Add the observed ICC to the grid (clipped to [0, 0.95)) so the plot
  ## marker has somewhere to land.
  icc_grid <- sort(unique(c(icc_grid,
                            pmin(0.95, observed_icc[!is.na(observed_icc)]))))

  ## ---- 5. simulate CLPM bias across icc_grid ------------------------------
  rows <- vector("list", length(icc_grid) * n_sims)
  k <- 0L

  for (icc in icc_grid) {
    var_BX <- if (icc == 0) 0 else icc / (1 - icc) * d_var_x
    var_BY <- if (icc == 0) 0 else icc / (1 - icc) * d_var_y

    for (s in seq_len(n_sims)) {
      k <- k + 1L
      est <- tryCatch({
        dat <- crossLagR::simRICLPM(
          waves        = waves,
          sample.nobs  = sample_size,
          beta_x       = riclpm_truth["ar_x"],
          beta_y       = riclpm_truth["ar_y"],
          omega_xy     = riclpm_truth["cl_xy"],
          omega_yx     = riclpm_truth["cl_yx"],
          var_p        = d_var_x, var_q  = d_var_y,
          cov_pq       = d_cov_xy,
          var_BX       = var_BX,  var_BY = var_BY,
          cov_BXBY     = 0
        )$data
        .sens_fit(dat, "CLPM", waves)
      }, error = function(e) {
        if (verbose) message("[icc=", icc, " sim=", s, "] ", e$message)
        rep(NA_real_, 4L)
      })
      rows[[k]] <- .sens_row(icc, s, "CLPM", est)
    }
    if (verbose) message("Finished ICC = ", icc)
  }

  est_df <- do.call(rbind, rows)
  est_df$riclpm_truth <- riclpm_truth[est_df$parameter]

  implied_bias <- clpm_actual - riclpm_truth

  out <- structure(
    list(
      estimates    = est_df,
      clpm_actual  = clpm_actual,
      riclpm_truth = riclpm_truth,
      observed_icc = observed_icc,
      implied_bias = implied_bias,
      icc_grid     = icc_grid,
      n_sims       = n_sims,
      sample_size  = sample_size,
      waves        = waves
    ),
    class = "crossLagR_sensitivity"
  )
  out
}


#' @export
print.crossLagR_sensitivity <- function(x, digits = 3, ...) {
  cat("crossLagR CLPM bias diagnostic\n",
      strrep("=", 32), "\n", sep = "")

  tab <- data.frame(
    parameter    = names(x$clpm_actual),
    clpm_actual  = round(x$clpm_actual,  digits),
    riclpm_truth = round(x$riclpm_truth, digits),
    implied_bias = round(x$implied_bias, digits),
    pct_inflation = ifelse(abs(x$riclpm_truth) > 1e-6,
                           round(100 * x$implied_bias / x$riclpm_truth, 1),
                           NA_real_),
    stringsAsFactors = FALSE
  )
  print(tab, row.names = FALSE)

  cat("\nObserved ICC in data:  x = ", round(x$observed_icc["x"], digits),
      ",  y = ", round(x$observed_icc["y"], digits), "\n", sep = "")
  cat("Sample size = ", x$sample_size, ",  waves = ", x$waves, ",  ",
      x$n_sims, " sims per ICC cell\n", sep = "")
  cat("\nInterpretation: rows where |implied_bias| is large relative to\n",
      "riclpm_truth indicate the CLPM is materially biased on this data.\n",
      "Call plot() to see the bias curve and where you sit on it.\n",
      sep = "")
  invisible(x)
}


#' @export
plot.crossLagR_sensitivity <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot.crossLagR_sensitivity().")
  }
  est <- x$estimates
  est$parameter <- factor(est$parameter,
                          levels = c("ar_x", "ar_y", "cl_xy", "cl_yx"))

  agg <- stats::aggregate(
    est ~ icc + parameter,
    data = est,
    FUN  = function(v) {
      v <- v[is.finite(v)]
      n <- length(v)
      mu <- mean(v)
      se <- if (n > 1L) stats::sd(v) / sqrt(n) else NA_real_
      c(mean = mu, lo = mu - 1.96 * se, hi = mu + 1.96 * se)
    }
  )
  agg <- do.call(data.frame, agg)
  names(agg)[3:5] <- c("mean", "lo", "hi")

  truth_df <- data.frame(parameter = names(x$riclpm_truth),
                         truth     = unname(x$riclpm_truth))
  truth_df$parameter <- factor(truth_df$parameter,
                               levels = levels(est$parameter))

  ## User's actual CLPM estimate at the observed ICC (use ICC_x for AR_x /
  ## CL paths involving x as the dependent, ICC_y for the others).
  obs_icc_for_param <- function(p) {
    switch(p,
           ar_x  = unname(x$observed_icc["x"]),
           cl_yx = unname(x$observed_icc["x"]),
           ar_y  = unname(x$observed_icc["y"]),
           cl_xy = unname(x$observed_icc["y"]),
           NA_real_)
  }
  marker_df <- data.frame(
    parameter = names(x$clpm_actual),
    icc       = vapply(names(x$clpm_actual), obs_icc_for_param, numeric(1)),
    value     = unname(x$clpm_actual),
    stringsAsFactors = FALSE
  )
  marker_df$parameter <- factor(marker_df$parameter,
                                levels = levels(est$parameter))

  vlines_df <- marker_df[, c("parameter", "icc")]

  ggplot2::ggplot(agg, ggplot2::aes(x = .data$icc, y = .data$mean)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
                         fill = "#b2182b", alpha = 0.20) +
    ggplot2::geom_line(color = "#b2182b", linewidth = 0.9) +
    ggplot2::geom_point(color = "#b2182b", size = 1.6) +
    ggplot2::geom_hline(data = truth_df,
                        ggplot2::aes(yintercept = .data$truth),
                        linetype = "dashed", color = "grey30") +
    ggplot2::geom_vline(data = vlines_df,
                        ggplot2::aes(xintercept = .data$icc),
                        linetype = "dotted", color = "grey50") +
    ggplot2::geom_point(data = marker_df,
                        ggplot2::aes(x = .data$icc, y = .data$value),
                        shape = 23, size = 4, fill = "#b2182b",
                        color = "black", stroke = 0.6,
                        inherit.aes = FALSE) +
    ggplot2::facet_wrap(~ .data$parameter, scales = "free_y") +
    ggplot2::labs(
      x = "ICC (between-person trait variance / total variance)",
      y = "Estimate",
      title    = "CLPM bias diagnostic",
      subtitle = paste0(
        "Dashed = RI-CLPM truth on your data. ",
        "Red curve = mean CLPM at each simulated ICC (",
        x$n_sims, " sims, ribbon = +/- 1.96 SE_MC). ",
        "Diamond = your actual CLPM estimate at the observed ICC."
      )
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"),
                   legend.position = "none")
}


## ----------------------------------------------------------------------------
## Internals
## ----------------------------------------------------------------------------

.sens_fit <- function(dat, estimator, waves) {
  syntax <- switch(
    estimator,
    "CLPM"   = crossLagR::estimateCLPM(waves   = waves),
    "RICLPM" = crossLagR::estimateRICLPM(waves = waves),
    stop("unknown sensitivity estimator: ", estimator)
  )
  fit <- tryCatch(
    suppressWarnings(lavaan::lavaan(syntax, data = dat, meanstructure = TRUE,
                                    control = list(iter.max = 500L))),
    error = function(e) NULL
  )
  if (is.null(fit) || !isTRUE(lavaan::lavInspect(fit, "converged"))) {
    return(rep(NA_real_, 4L))
  }
  pe <- lavaan::parameterEstimates(fit)
  pick <- function(lab) {
    v <- pe$est[pe$label == lab]
    if (!length(v)) NA_real_ else as.numeric(v[1])
  }
  c(pick("ar_x"), pick("ar_y"), pick("cl_xy"), pick("cl_yx"))
}

.sens_row <- function(icc, sim, estimator, est_vec) {
  data.frame(
    icc       = icc,
    sim       = sim,
    estimator = estimator,
    parameter = c("ar_x", "ar_y", "cl_xy", "cl_yx"),
    est       = est_vec,
    stringsAsFactors = FALSE
  )
}
