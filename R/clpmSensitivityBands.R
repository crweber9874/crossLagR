# Simulation-based sensitivity bands for a reported CLPM. Generalises the
# figure in the identification chapter: sweep the between-to-within variance
# ratio, simulate an RI-CLPM at each value with a true cross-lag of zero, fit a
# CLPM to every draw, and show where the reported estimates fall against the
# artifact-only curves.

utils::globalVariables(c(
  "AR", "AR_lo", "AR_hi", "CL", "CL_lo", "CL_hi",
  "parameter", "mean_est", "lo", "hi", "value"
))

#' Styling for the CLPM sensitivity band plot
#'
#' Collects every visual choice in one place so callers can restyle the figure
#' without editing the plotting code. Pass the result to the \code{style}
#' argument of \code{\link{clpmSensitivityBands}}.
#'
#' @param band_fill Fill colour for the simulation ribbon.
#' @param band_alpha Ribbon opacity.
#' @param line_colour Colour of the mean artifact curve.
#' @param line_width Width of the mean artifact curve.
#' @param ref_colour,ref_linetype Colour and line type of the horizontal
#'   reference lines marking the reported estimates.
#' @param point_colour,point_size Colour and size of the markers placed where
#'   an artifact curve crosses its reported value. \code{NULL} suppresses them.
#' @param base_size Base font size passed to \code{ggplot2::theme_minimal()}.
#' @param ar_label,cl_label Facet titles for the two panels.
#' @param x_lab,y_lab Axis titles.
#' @return A named list of styling values.
#' @export
clpmBandStyle <- function(band_fill    = "#C0392B",
                          band_alpha   = 0.18,
                          line_colour  = "#C0392B",
                          line_width   = 0.9,
                          ref_colour   = "grey30",
                          ref_linetype = "dashed",
                          point_colour = "#2C3E50",
                          point_size   = 2,
                          base_size    = 11,
                          ar_label     = "Autoregression",
                          cl_label     = "Cross-lag",
                          x_lab        = "Between-to-within variance ratio (r)",
                          y_lab        = "CLPM estimate under a true cross-lag of zero") {
  list(
    band_fill = band_fill, band_alpha = band_alpha,
    line_colour = line_colour, line_width = line_width,
    ref_colour = ref_colour, ref_linetype = ref_linetype,
    point_colour = point_colour, point_size = point_size,
    base_size = base_size, ar_label = ar_label, cl_label = cl_label,
    x_lab = x_lab, y_lab = y_lab
  )
}

# One simulate-and-fit cycle at a given variance ratio.
.clpm_band_draw <- function(ar_true, bw_ratio, rho_i, n_obs, waves) {
  d <- simRICLPM(
    waves       = waves,
    beta_x      = ar_true, beta_y = ar_true,
    omega_xy    = 0,       omega_yx = 0,
    var_p       = 1 - ar_true^2, var_q = 1 - ar_true^2, cov_pq = 0,
    var_BX      = bw_ratio, var_BY = bw_ratio,
    cov_BXBY    = bw_ratio * rho_i,
    sample.nobs = n_obs
  )
  pe <- lavaan::parameterEstimates(
    lavaan::lavaan(estimateCLPM(waves = waves), data = d$data)
  )
  pe <- pe[pe$label != "", ]
  data.frame(
    ar     = mean(pe$est[pe$label == "ar_x"]),
    cl     = mean(pe$est[pe$label == "cl_xy"]),
    n_obs  = nrow(d$data),
    n_miss = sum(!stats::complete.cases(d$data))
  )
}

#' Simulate CLPM sensitivity bands for a reported estimate
#'
#' Sweeps the between-to-within variance ratio, simulates an RI-CLPM at each
#' value with a true cross-lag of zero, fits a CLPM to every replication, and
#' returns the artifact-only curves with simulation intervals. Horizontal
#' reference lines mark the reported estimates, so the crossings show which
#' amounts of trait variance could have produced the report with no
#' within-person effect at all.
#'
#' The simulation is the expensive part, and it is evaluated lazily: supply a
#' previously returned \code{sim} data frame (or the \code{"sim"} attribute of
#' a returned plot) and no data are generated, so a figure can be restyled or
#' re-faceted without paying for the sweep again.
#'
#' @param reported_ar,reported_cl The published CLPM autoregression and
#'   cross-lag. Drawn as horizontal reference lines.
#' @param ar_true True within-person autoregression used to generate data.
#' @param rho_i Trait correlation used to generate data.
#' @param r_grid Between-to-within variance ratios to sweep.
#' @param reps Replications per grid point.
#' @param n_obs Persons per replication.
#' @param waves Number of waves.
#' @param conf Width of the simulation interval, e.g. \code{0.95}.
#' @param seed Optional integer seed, set once before the sweep.
#' @param sim Optional pre-computed sweep from an earlier call. When supplied
#'   the simulation is skipped entirely.
#' @param plot If \code{FALSE}, return the swept data frame instead of a plot.
#' @param style A list of visual settings from \code{\link{clpmBandStyle}}.
#' @return A \code{ggplot} object with the swept data attached as the
#'   \code{"sim"} attribute, or that data frame when \code{plot = FALSE}. The
#'   data frame carries \code{bw_ratio}, the mean and interval bounds for both
#'   parameters, and \code{n_obs} / \code{n_miss} per grid point.
#' @examples
#' \donttest{
#' p <- clpmSensitivityBands(reported_ar = 0.55, reported_cl = 0.12,
#'                           reps = 5, n_obs = 300, seed = 1)
#' # restyle without re-simulating
#' clpmSensitivityBands(0.55, 0.12, sim = attr(p, "sim"),
#'                      style = clpmBandStyle(band_fill = "steelblue"))
#' }
#' @export
clpmSensitivityBands <- function(reported_ar, reported_cl,
                                 ar_true = 0.30,
                                 rho_i   = 0.50,
                                 r_grid  = c(0.05, 0.25, 0.5, 0.75, 1,
                                             1.5, 2, 3, 4.5, 6, 9),
                                 reps    = 40,
                                 n_obs   = 1000,
                                 waves   = 4,
                                 conf    = 0.95,
                                 seed    = NULL,
                                 sim     = NULL,
                                 plot    = TRUE,
                                 style   = clpmBandStyle()) {
  stopifnot(
    is.numeric(reported_ar), is.numeric(reported_cl),
    length(reported_ar) == 1, length(reported_cl) == 1,
    is.numeric(conf), conf > 0, conf < 1, is.list(style)
  )

  if (is.null(sim)) {
    if (!requireNamespace("lavaan", quietly = TRUE)) {
      stop("Package 'lavaan' is required for clpmSensitivityBands().", call. = FALSE)
    }
    stopifnot(is.numeric(r_grid), length(r_grid) >= 1, all(r_grid > 0),
              reps >= 1, n_obs >= 1, waves >= 2)
    if (!is.null(seed)) set.seed(seed)

    lo_p <- (1 - conf) / 2
    hi_p <- 1 - lo_p
    sim <- do.call(rbind, lapply(r_grid, function(r) {
      draws <- do.call(rbind, lapply(
        seq_len(reps),
        function(i) .clpm_band_draw(ar_true, r, rho_i, n_obs, waves)
      ))
      data.frame(
        bw_ratio = r,
        AR    = mean(draws$ar),
        AR_lo = unname(stats::quantile(draws$ar, lo_p)),
        AR_hi = unname(stats::quantile(draws$ar, hi_p)),
        CL    = mean(draws$cl),
        CL_lo = unname(stats::quantile(draws$cl, lo_p)),
        CL_hi = unname(stats::quantile(draws$cl, hi_p)),
        n_obs  = sum(draws$n_obs),
        n_miss = sum(draws$n_miss)
      )
    }))
  }

  if (!isTRUE(plot)) return(sim)

  bands <- rbind(
    data.frame(bw_ratio = sim$bw_ratio, parameter = style$ar_label,
               mean_est = sim$AR, lo = sim$AR_lo, hi = sim$AR_hi),
    data.frame(bw_ratio = sim$bw_ratio, parameter = style$cl_label,
               mean_est = sim$CL, lo = sim$CL_lo, hi = sim$CL_hi)
  )
  refs <- data.frame(
    parameter = c(style$ar_label, style$cl_label),
    value     = c(reported_ar, reported_cl)
  )

  p <- ggplot2::ggplot(bands, ggplot2::aes(bw_ratio, mean_est)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi),
                         fill = style$band_fill, alpha = style$band_alpha) +
    ggplot2::geom_line(linewidth = style$line_width, colour = style$line_colour) +
    ggplot2::geom_hline(data = refs, ggplot2::aes(yintercept = value),
                        linetype = style$ref_linetype, colour = style$ref_colour) +
    ggplot2::facet_wrap(~parameter, scales = "free_y") +
    ggplot2::labs(x = style$x_lab, y = style$y_lab) +
    ggplot2::theme_minimal(base_size = style$base_size)

  if (!is.null(style$point_colour)) {
    cross <- .clpm_band_crossings(sim, reported_ar, reported_cl, style)
    if (nrow(cross)) {
      p <- p + ggplot2::geom_point(
        data = cross, ggplot2::aes(bw_ratio, mean_est),
        colour = style$point_colour, size = style$point_size
      )
    }
  }

  attr(p, "sim") <- sim
  p
}

# Linear interpolation of the points where each artifact curve meets its
# reported value. Purely for marking the figure.
.clpm_band_crossings <- function(sim, reported_ar, reported_cl, style) {
  one <- function(x, y, target, label) {
    d <- y - target
    idx <- which(d[-length(d)] * d[-1] < 0)
    if (!length(idx)) return(NULL)
    xs <- vapply(idx, function(i) {
      x[i] + (x[i + 1] - x[i]) * (target - y[i]) / (y[i + 1] - y[i])
    }, numeric(1))
    data.frame(bw_ratio = xs, mean_est = target, parameter = label)
  }
  out <- rbind(
    one(sim$bw_ratio, sim$AR, reported_ar, style$ar_label),
    one(sim$bw_ratio, sim$CL, reported_cl, style$cl_label)
  )
  if (is.null(out)) data.frame() else out
}
