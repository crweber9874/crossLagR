#' @title sensitivityCLPM
#' @description Sensitivity analysis for a fitted CLPM. Treats the CLPM's
#'   estimated AR/CL parameters as the (hypothetical) within-person truth,
#'   then injects between-person trait variance at a grid of ICC levels and
#'   refits both a CLPM and an RI-CLPM. The output diagnoses how badly a
#'   CLPM would mis-estimate the AR/CL parameters under varying degrees of
#'   stable trait variance --- the bias mechanism formalized by
#'   @hamaker2015critique and derived analytically in the CLPM chapter of
#'   the package book.
#'
#' @param fit A fitted lavaan CLPM object (from \code{lavaan::lavaan} or
#'   \code{lavaan::sem}) whose syntax was produced by
#'   \code{\link{estimateCLPM}}.
#' @param icc_grid Numeric vector. Hypothetical ICC values to test
#'   (proportion of total variance from a stable trait). Default
#'   \code{seq(0, 0.8, by = 0.1)}; values \eqn{\ge 0.9} routinely fail to
#'   converge.
#' @param n_sims Integer. Replications per ICC level. Default 30 (sensible
#'   for a sensitivity check; bump to 100+ for publication).
#' @param sample_size Integer. Sample size per simulated dataset. Default:
#'   match the original fit's \code{lavInspect(fit, "nobs")}.
#' @param waves Integer. Number of waves. Default: inferred from the
#'   original fit's observed variable list.
#' @param seed Optional integer for reproducibility.
#' @param verbose Logical. Print per-cell progress. Default \code{FALSE}.
#'
#' @return An object of class \code{crossLagR_sensitivity} with elements:
#'   \itemize{
#'     \item \code{estimates}: long data.frame of recovered AR/CL estimates
#'       across ICC x sim x estimator (CLPM vs RI-CLPM).
#'     \item \code{original}: named numeric vector of the user's baseline
#'       CLPM estimates (\code{ar_x, ar_y, cl_xy, cl_yx}).
#'     \item \code{icc_grid}, \code{n_sims}, \code{sample_size}, \code{waves}.
#'   }
#'   Has \code{print()} and \code{plot()} methods.
#'
#' @details
#' \strong{What this answers.} "If the true data-generating process for my
#' data has within-person dynamics equal to the CLPM estimates I observed,
#' but ICC of the total variance comes from a stable trait, how biased
#' would my CLPM estimates have been --- and what would an RI-CLPM have
#' recovered instead?"
#'
#' \strong{What this does not answer.} It does not calibrate the user's
#' CLPM estimates against an unknown truth. The CLPM estimates already
#' embed any trait-induced bias present in the original data; treating them
#' as "within-person truth" is a sensitivity device, not a correction. To
#' see the actually-implied within-person dynamics, fit an RI-CLPM directly.
#'
#' \strong{The bias mechanism.} The univariate analytical bias is
#' \deqn{b^{\text{CLPM}}_{xx} \approx \text{ICC}_x + (1 - \text{ICC}_x)\,b^{\text{true}}_{xx}}
#' (see Sec. on bias in the CLPM chapter). For \eqn{b^{\text{true}}_{xx} = 0}
#' the CLPM produces a spurious AR equal to \eqn{\text{ICC}_x}. The
#' multivariate case lacks a clean closed form because the coefficients are
#' partial slopes; this function characterizes it by simulation.
#'
#' @references
#' Hamaker, E. L., Kuiper, R. M., & Grasman, R. P. P. P. (2015). A critique
#' of the cross-lagged panel model. \emph{Psychological Methods}, 20(1),
#' 102--116.
#'
#' @examples
#' \dontrun{
#' dat <- simCLPM(waves = 4, sample_size = 500, seed = 1)$data
#' fit <- lavaan::lavaan(estimateCLPM(waves = 4), data = dat,
#'                       meanstructure = TRUE)
#' s <- sensitivityCLPM(fit, icc_grid = seq(0, 0.6, 0.15), n_sims = 20)
#' s
#' plot(s)
#' }
#' @export
sensitivityCLPM <- function(fit,
                            icc_grid = seq(0, 0.8, by = 0.1),
                            n_sims = 30,
                            sample_size = NULL,
                            waves = NULL,
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

  pe <- lavaan::parameterEstimates(fit)
  pick <- function(lab) {
    v <- pe$est[pe$label == lab]
    if (!length(v)) NA_real_ else as.numeric(v[1])
  }
  orig <- c(
    ar_x  = pick("ar_x"),  ar_y  = pick("ar_y"),
    cl_xy = pick("cl_xy"), cl_yx = pick("cl_yx")
  )
  d_var_x  <- pick("d_var_x")
  d_var_y  <- pick("d_var_y")
  d_cov_xy <- pick("d_cov_xy")

  if (any(is.na(orig)) || is.na(d_var_x) || is.na(d_var_y)) {
    stop("Could not extract CLPM parameters from fit. ",
         "Was the syntax built with estimateCLPM()?")
  }

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

  rows <- vector("list", length(icc_grid) * n_sims)
  k <- 0L

  for (icc in icc_grid) {
    var_BX <- if (icc == 0) 0 else icc / (1 - icc) * d_var_x
    var_BY <- if (icc == 0) 0 else icc / (1 - icc) * d_var_y

    for (s in seq_len(n_sims)) {
      k <- k + 1L
      cell <- tryCatch({
        dat <- crossLagR::simRICLPM(
          waves        = waves,
          sample.nobs  = sample_size,
          beta_x       = orig["ar_x"],  beta_y = orig["ar_y"],
          omega_xy     = orig["cl_xy"], omega_yx = orig["cl_yx"],
          var_p        = d_var_x,       var_q  = d_var_y,
          cov_pq       = d_cov_xy,
          var_BX       = var_BX,        var_BY = var_BY,
          cov_BXBY     = 0
        )$data
        list(
          clpm   = .sens_fit(dat, "CLPM",   waves),
          riclpm = .sens_fit(dat, "RICLPM", waves)
        )
      }, error = function(e) {
        if (verbose) message("[icc=", icc, " sim=", s, "] ", e$message)
        NULL
      })

      if (is.null(cell)) {
        rows[[k]] <- .sens_row(icc, s, "CLPM",   rep(NA_real_, 4))
        next
      }
      rows[[k]] <- rbind(
        .sens_row(icc, s, "CLPM",    cell$clpm),
        .sens_row(icc, s, "RI-CLPM", cell$riclpm)
      )
    }
    if (verbose) message("Finished ICC = ", icc)
  }

  est <- do.call(rbind, rows)
  est$original_est <- orig[est$parameter]

  out <- structure(
    list(
      estimates   = est,
      original    = orig,
      icc_grid    = icc_grid,
      n_sims      = n_sims,
      sample_size = sample_size,
      waves       = waves
    ),
    class = "crossLagR_sensitivity"
  )
  out
}


#' @export
print.crossLagR_sensitivity <- function(x, digits = 3, ...) {
  cat("crossLagR CLPM sensitivity analysis\n",
      strrep("=", 35), "\n", sep = "")
  cat("Baseline CLPM estimates (treated as within-person truth):\n")
  print(round(x$original, digits))
  cat("\nGrid: ICC in {", paste(x$icc_grid, collapse = ", "), "}, ",
      x$n_sims, " sims per cell, N = ", x$sample_size, ", waves = ", x$waves,
      "\n", sep = "")

  agg <- stats::aggregate(
    est ~ icc + estimator + parameter,
    data = x$estimates,
    FUN  = function(v) c(mean = mean(v, na.rm = TRUE),
                         sd   = stats::sd(v, na.rm = TRUE))
  )
  agg <- do.call(data.frame, agg)
  names(agg)[4:5] <- c("mean_est", "sd_est")
  cat("\nMean estimate by ICC x estimator x parameter:\n")
  cat(strrep("-", 50), "\n", sep = "")
  agg$mean_est <- round(agg$mean_est, digits)
  agg$sd_est   <- round(agg$sd_est,   digits)
  print(agg, row.names = FALSE)

  cat("\nCall plot() on this object for the sensitivity diagram.\n")
  invisible(x)
}


#' @export
plot.crossLagR_sensitivity <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot.crossLagR_sensitivity().")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required for plot.crossLagR_sensitivity().")
  }
  est <- x$estimates
  est$parameter <- factor(est$parameter,
                          levels = c("ar_x", "ar_y", "cl_xy", "cl_yx"))

  agg <- stats::aggregate(
    est ~ icc + estimator + parameter,
    data = est,
    FUN  = function(v) c(mean = mean(v, na.rm = TRUE),
                         lo   = stats::quantile(v, 0.10, na.rm = TRUE),
                         hi   = stats::quantile(v, 0.90, na.rm = TRUE))
  )
  agg <- do.call(data.frame, agg)
  names(agg)[4:6] <- c("mean", "lo", "hi")

  baseline <- data.frame(parameter = names(x$original),
                         baseline  = unname(x$original))

  ggplot2::ggplot(agg, ggplot2::aes(x = .data$icc, y = .data$mean,
                                    color = .data$estimator,
                                    fill  = .data$estimator)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
                         alpha = 0.18, color = NA) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::geom_hline(data = baseline,
                        ggplot2::aes(yintercept = .data$baseline),
                        linetype = "dashed", color = "grey30") +
    ggplot2::facet_wrap(~ .data$parameter, scales = "free_y") +
    ggplot2::scale_color_manual(values = c("CLPM"    = "#b2182b",
                                           "RI-CLPM" = "#2166ac")) +
    ggplot2::scale_fill_manual(values  = c("CLPM"    = "#b2182b",
                                           "RI-CLPM" = "#2166ac")) +
    ggplot2::labs(
      x = "Injected ICC (trait variance / total variance)",
      y = "Recovered estimate (mean +/- 10th-90th percentile)",
      color = NULL, fill = NULL,
      title    = "CLPM sensitivity to between-person trait variance",
      subtitle = paste0(
        "Dashed line: original CLPM estimate (preserved as within-person truth). ",
        "N = ", x$sample_size, ", ", x$n_sims, " sims per ICC."
      )
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "top",
                   strip.text      = ggplot2::element_text(face = "bold"))
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
