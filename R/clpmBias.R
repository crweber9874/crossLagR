# Helpers for the CLPM trait-bias figures in 02_Bias.qmd. cl_bias_peak_ratio,
# plot_cl_bias_peak and plot_cl_bias_surface are called directly from the book
# and are exported; the rest are internal implementation details.
# Functions used only by later chapters (the lavaan bias-landscape machinery
# for 08_liss.qmd) live in R/clpmBiasLandscape.R instead.

utils::globalVariables(c(
  "bw_ratio", "cl_bias", "ar", "ar_bias", "rho_i",
  "reliability", "cl_pop", "ar_bias_pct", "cl_bias_pct"
))

#' Population CLPM slope bias from a random-intercept DGP (general form)
#'
#' CLPM bias is due to the nonzero correlation between the endogenous regressors and the omitted traits.
#' This function calculates the bias matrix for a general k-variable CLPM,
#' given the trait covariance matrix \eqn{\boldsymbol\Omega},
#' the within-person covariance matrix \eqn{\boldsymbol\Psi},
#' and the true within-person dynamics \eqn{\mathbf B^{\text{true}}}.
#'
#' @param omega Trait covariance matrix (k x k).
#' @param psi Within-person covariance matrix (k x k).
#' @param b Matrix of true within-person dynamics (k x k).
#' @return A k x k matrix of CLPM slope biases.
#' @keywords internal
clpm_bias_matrix <- function(omega, psi, b) {
  stopifnot(
    is.matrix(omega), is.matrix(psi), is.matrix(b),
    nrow(omega) == ncol(omega),
    all(dim(omega) == dim(psi)), all(dim(omega) == dim(b))
  )
  k <- nrow(omega)
  (diag(k) - b) %*% omega %*% solve(omega + psi)
}

#' Variance ratio at which the cross-lagged bias peaks
#'
#' @param rho_i Trait correlation \eqn{\rho_I}.
#' @return The ratio \eqn{r^\star = 1/\sqrt{1 - \rho_I^2}}.
#' @export
cl_bias_peak_ratio <- function(rho_i = 0.50) {
  1 / sqrt(1 - rho_i^2)
}

#' Plot the cross-lagged bias against the between-to-within variance ratio
#'
#' Traces the CL bias as the variance ratio grows, marking the peak at
#' \eqn{r^\star = 1/\sqrt{1 - \rho_I^2}}. Parameters come from the simulation
#' setup.
#'
#' @param ar True autoregressive parameter.
#' @param rho_i Trait correlation.
#' @param r_max Largest variance ratio to plot.
#' @param n Number of grid points.
#' @return A \code{ggplot} object.
#' @export
plot_cl_bias_peak <- function(ar = 0.30, rho_i = 0.50, r_max = 12, n = 500) {
  r_star <- cl_bias_peak_ratio(rho_i)
  df <- clpm_artifact_grid(seq(0, r_max, length.out = n), ar = ar, rho_i = rho_i)
  peak <- clpm_artifact_grid(r_star, ar = ar, rho_i = rho_i)$cl_bias

  ggplot2::ggplot(df, ggplot2::aes(bw_ratio, cl_bias)) +
    ggplot2::geom_line(color = "#C0392B", linewidth = 1) +
    ggplot2::geom_vline(xintercept = r_star, linetype = "dashed", color = "grey55") +
    ggplot2::geom_point(
      data = data.frame(bw_ratio = r_star, cl_bias = peak),
      color = "#C0392B", size = 2.6
    ) +
    ggplot2::annotate(
      "text",
      x = r_star + 0.35, y = peak, hjust = 0, size = 3.4, color = "grey25",
      label = sprintf("peak: r* = %.2f\nCL bias = %.3f", r_star, peak)
    ) +
    ggplot2::scale_x_continuous(breaks = seq(0, r_max, by = 2)) +
    ggplot2::labs(
      x = "Between-to-within variance ratio  (r)",
      y = "Cross-lagged bias"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}

#' Plot the cross-lagged bias surface over variance ratio and true AR
#'
#'
#' @param rho_i Trait correlation.
#' @param ar_range Length-2 numeric range of true AR values for the y axis.
#' @param r_max Largest variance ratio to plot.
#' @param n_r,n_ar Grid resolution along the ratio and AR axes.
#' @param palette Either a viridis option (e.g. \code{"mako"}, \code{"rocket"},
#'   \code{"magma"}, \code{"viridis"}) or \code{"grey"} / \code{"greys"} for
#'   a greyscale fill.
#' @param style One of \code{"raster"} (smooth gradient) or \code{"topo"}
#'   (filled iso-bias bands, like a topographic map).
#' @param contours Logical; overlay iso-bias contour lines (always on in
#'   \code{"topo"} style).
#' @param bins Number of contour bands or lines.
#' @return A \code{ggplot} object.
#' @export
plot_cl_bias_surface <- function(rho_i = 0.50, ar_range = c(0, 0.95),
                                 r_max = 12, n_r = 250, n_ar = 200,
                                 palette = "mako",
                                 style = c("raster", "topo"),
                                 contours = TRUE, bins = 12) {
  style <- match.arg(style)
  is_grey <- tolower(palette) %in% c("grey", "greys", "gray", "grays")
  ridge_col <- if (is_grey || style == "topo") "#1A1A1A" else "white"

  r_star <- cl_bias_peak_ratio(rho_i)
  grid <- expand.grid(
    bw_ratio = seq(0, r_max, length.out = n_r),
    ar       = seq(ar_range[1], ar_range[2], length.out = n_ar)
  )
  df <- clpm_artifact_grid(grid$bw_ratio, ar = grid$ar, rho_i = rho_i)

  p <- ggplot2::ggplot(df, ggplot2::aes(bw_ratio, ar, z = cl_bias))

  if (style == "topo") {
    p <- p +
      ggplot2::geom_contour_filled(bins = bins) +
      ggplot2::geom_contour(
        color = ridge_col, alpha = 0.6, linewidth = 0.25, bins = bins
      )
    if (is_grey) {
      p <- p +
        ggplot2::scale_fill_grey(
          name = "CL bias",
          start = 0.95, end = 0.20
        )
    } else {
      p <- p +
        ggplot2::scale_fill_viridis_d(name = "CL bias", option = palette)
    }
  } else {
    p <- p + ggplot2::geom_raster(ggplot2::aes(fill = cl_bias),
      interpolate = TRUE
    )
    if (isTRUE(contours)) {
      p <- p +
        ggplot2::geom_contour(
          color = ridge_col, alpha = 0.4, linewidth = 0.25, bins = bins
        )
    }
    if (is_grey) {
      p <- p +
        ggplot2::scale_fill_gradient(
          name = "CL bias",
          low = "white", high = "grey15"
        )
    } else {
      p <- p +
        ggplot2::scale_fill_viridis_c(name = "CL bias", option = palette)
    }
  }

  p +
    ggplot2::geom_vline(
      xintercept = r_star, linetype = "dashed",
      color = ridge_col, linewidth = 0.6
    ) +
    ggplot2::annotate(
      "text",
      x = r_star + 0.25, y = ar_range[2] - 0.05, hjust = 0,
      color = ridge_col, size = 3.2,
      label = sprintf("ridge at r* = %.2f", r_star)
    ) +
    ggplot2::scale_x_continuous(breaks = seq(0, r_max, by = 2)) +
    ggplot2::labs(
      x = "Between-to-within variance ratio  (r)",
      y = "True autoregressive parameter  (b)"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}

#' Cross-lagged and autoregressive bias of the CLPM
#'
#' .
#'
#' The closed form is
#' \deqn{\mathrm{CL\ bias} = (1-b)\,\rho_I\, r / [(1-\rho_I^2) r^2 + 2r + 1],}
#' \deqn{\mathrm{AR\ bias} = (1-b)\,[(1-\rho_I^2) r^2 + r] / [(1-\rho_I^2) r^2 + 2r + 1],}
#' with \eqn{r} the between-to-within variance ratio. It assumes equal trait
#' variances, equal autoregressions, no true cross-lag, and diagonal within
#' errors terms.
#'
#' @param bw_ratio Between-to-within trait variance ratio
#'   \eqn{r = \sigma^2_I / \sigma^2_{\mathrm{within}}}. Non-negative; Vector -- See Bias Chapter
#' @param ar True autoregressive parameter \eqn{b}. A vector -- see Bias Chapter
#' @param rho_i Trait correlation \eqn{\rho_I}, in \eqn{[-1, 1]}. Vectorised.
#' @param expand If \code{TRUE}, evaluate on the full crossing (grid) of the
#'   supplied values; if \code{FALSE} (default), recycle them to a common length.
#' @return A data frame, one row per evaluated combination, with columns
#'   \code{bw_ratio}, \code{ar}, \code{rho_i}, \code{icc}, \code{ar_bias},
#'   \code{ar_bias_pct}, \code{cl_bias}, \code{cl_bias_pct}. Percentage
#'   columns express bias as a share of the assumed truth (\code{100 * bias /
#'   truth}); since the true CL is 0 here, \code{cl_bias_pct} is always
#'   \code{NA} (division by zero).
#' @examples
#' # Cross-lagged bias across a range of variance ratios
#' clpmBias(bw_ratio = c(0.5, 1, 2, 4, 9), ar = 0.3, rho_i = 0.5)
#'
#' # Vary the true autoregression and the trait correlation on a grid
#' clpmBias(bw_ratio = c(1, 2), ar = c(0, 0.3, 0.6), rho_i = 0.5, expand = TRUE)
#' @export
clpmBias <- function(bw_ratio = 1, ar = 0.30, rho_i = 0.50, expand = FALSE) {
  stopifnot(
    is.numeric(bw_ratio), is.numeric(ar), is.numeric(rho_i),
    is.logical(expand)
  )
  if (any(bw_ratio < 0)) stop("`bw_ratio` must be non-negative.", call. = FALSE)
  if (any(abs(rho_i) > 1)) stop("`rho_i` must lie in [-1, 1].", call. = FALSE)

  if (isTRUE(expand)) {
    g <- expand.grid(bw_ratio = bw_ratio, ar = ar, rho_i = rho_i)
    bw_ratio <- g$bw_ratio
    ar <- g$ar
    rho_i <- g$rho_i
  }

  out <- clpm_artifact_grid(bw_ratio, ar = ar, rho_i = rho_i)[c(
    "bw_ratio", "ar", "rho_i", "ar_bias", "cl_bias", "ar_bias_pct", "cl_bias_pct"
  )]
  out$icc <- out$bw_ratio / (1 + out$bw_ratio)
  out[c(
    "bw_ratio", "ar", "rho_i", "icc",
    "ar_bias", "ar_bias_pct",
    "cl_bias", "cl_bias_pct"
  )]
}

# ---- As-if artifact: pop CLPM estimate under an RI-CLPM truth ---------------
#
# Generalises the bias grid to a non-zero TRUE within-person cross-lag c.
# Under a stationary RI-CLPM truth with symmetric B = [[b, c], [c, b]],
# innovation Psi_w = sigma_w^2 * I, trait Omega = omega_trait * [[1, rho_i],
# [rho_i, 1]], and optional classical measurement error of reliability w,
# returns the population CLPM regression slope (the estimate a CLPM analyst
# would obtain in the limit). The user-facing question is "as if" sensitivity:
# at different combinations of trait-variance share and assumed true (ar, cl),
# how much of an observed CLPM cross-lag is artifact?
#
# Eigen-decomp of symmetric B with U = (1/sqrt(2))[[1,1],[1,-1]] gives
# eigenvalues lambda1 = b+c, lambda2 = b-c, so the stationary within
# covariance has equal diagonals and an off-diagonal driven by the difference
# between (1 - lambda1^2)^-1 and (1 - lambda2^2)^-1. We parameterise so
# stationary within DIAGONAL variance = 1, matching the no-cross-lag case.

#' Stationary within covariance for the symmetric two-series RI-CLPM
#'
#' @param ar True autoregression (scalar or vector recyclable).
#' @param cl True within-person cross-lag.
#' @return List with \code{psi_diag} (= 1 by convention), \code{rho_w}
#'   (within-person stationary correlation induced by \code{cl}), and
#'   \code{sigma2_w} (innovation variance compatible with \code{psi_diag = 1}).
#' @keywords internal
.psi_stationary_symmetric <- function(ar, cl) {
  lp <- ar + cl
  lm <- ar - cl
  if (any(abs(lp) >= 1 | abs(lm) >= 1)) {
    stop("Non-stationary B: require |ar + cl| < 1 and |ar - cl| < 1.",
      call. = FALSE
    )
  }
  ipp <- 1 / (1 - lp^2)
  imm <- 1 / (1 - lm^2)
  sigma2_w <- 2 / (ipp + imm)
  rho_w <- (ipp - imm) / (ipp + imm)
  list(psi_diag = rep(1, length(ar)), rho_w = rho_w, sigma2_w = sigma2_w)
}

#' Population CLPM estimate under an RI-CLPM truth (symmetric two-series)
#'
#' Closed form for the population CLPM slope a researcher would obtain when
#' the true DGP is a stationary RI-CLPM with symmetric within-person dynamics
#' \eqn{B = [[ar, cl], [cl, ar]]}. The CLPM analyst ignores the trait
#' decomposition and recovers
#' \deqn{\Phi^{\text{pop}} = \Sigma_{\text{lag1}}\,\Sigma_{\text{var}}^{-1},}
#' where \eqn{\Sigma_{\text{var}} = \Omega + \Psi + \Sigma_u} (with optional
#' measurement error) and \eqn{\Sigma_{\text{lag1}} = \Omega + B\Psi}.
#' The stationary within \eqn{\Psi} is induced by \code{ar} and \code{cl}; see
#' \code{.psi_stationary_symmetric()}.
#'
#'
#' @param bw_ratio Between-to-within variance ratio (trait variance over
#'   stationary within diagonal variance). Vectorised.
#' @param ar True autoregressive parameter. Vectorised.
#' @param cl True within-person cross-lag. Vectorised.
#' @param rho_i Trait correlation. Vectorised.
#' @param reliability Indicator reliability \eqn{\omega \in (0, 1]}. Vectorised.
#' @return Data frame with the population CLPM estimates \code{ar_pop},
#'   \code{cl_pop} and the corresponding biases \code{ar_bias}, \code{cl_bias}.
#' @keywords internal
clpm_artifact_grid <- function(bw_ratio, ar = 0.30, cl = 0,
                               rho_i = 0.50, reliability = 1) {
  args <- data.frame(
    bw_ratio = bw_ratio, ar = ar, cl = cl,
    rho_i = rho_i, reliability = reliability
  )
  psi <- .psi_stationary_symmetric(args$ar, args$cl)
  rho_w <- psi$rho_w

  r <- args$bw_ratio
  b <- args$ar
  cc <- args$cl
  p <- args$rho_i
  w <- args$reliability
  u <- (r + 1) * (1 - w) / w

  a <- r + b + cc * rho_w
  g <- r * p + b * rho_w + cc
  v <- r + 1 + u
  h <- r * p + rho_w
  det_s <- v^2 - h^2

  ar_pop <- (a * v - g * h) / det_s
  cl_pop <- (g * v - a * h) / det_s
  ar_bias <- ar_pop - b
  cl_bias <- cl_pop - cc

  # Percentage bias relative to the assumed truth: `bias / truth * 100`.
  # When the truth is zero the ratio is undefined (division by zero), so it
  # is NA -- the same convention used for ar_bias_pct. We recycle b and cc to
  # ar_bias's length first so the masks behave correctly when scalars are
  # passed.
  b_full <- rep_len(b, length(ar_bias))
  cc_full <- rep_len(cc, length(cl_bias))
  ar_bias_pct <- 100 * ar_bias / b_full
  ar_bias_pct[b_full == 0] <- NA_real_
  cl_bias_pct <- 100 * cl_bias / cc_full
  cl_bias_pct[cc_full == 0] <- NA_real_

  data.frame(
    bw_ratio    = r,
    ar          = b,
    cl          = cc,
    rho_i       = p,
    reliability = w,
    rho_w       = rho_w,
    ar_pop      = ar_pop,
    cl_pop      = cl_pop,
    ar_bias     = ar_bias,
    cl_bias     = cl_bias,
    ar_bias_pct = ar_bias_pct,
    cl_bias_pct = cl_bias_pct
  )
}

#' Tipping-point trait-variance ratio that produces a target CLPM cross-lag
#'
#' For a given true within-person cross-lag \code{cl}, finds the
#' between-to-within variance ratio at which the population CLPM cross-lag
#' equals \code{target}. Returns \code{NA} when no non-negative solution
#' exists in \code{interval} (e.g., \code{target} is unreachable).
#'
#' \code{cl_pop(r)} in \code{clpm_artifact_grid()} is a ratio of two
#' quadratics in \code{r} (numerator \code{g(r)v(r) - a(r)h(r)}, denominator
#' \code{v(r)^2 - h(r)^2}, with \code{a, g, v, h} all linear in \code{r}).
#' Setting \code{cl_pop(r) = target} and clearing the denominator therefore
#' collapses to a single quadratic \eqn{A r^2 + Br + C = 0} in \code{r}, so
#' the tipping point has a closed-form solution via the quadratic formula --
#' no root-finding needed. When two non-negative roots exist (the curve can
#' cross a target twice, e.g. once on the way up and once on the way down),
#' we report the smaller one: the least trait variance able to produce the
#' target as pure artifact.
#'
#' @param target Numeric scalar, the observed (or hypothesised) CLPM cross-lag.
#' @param ar True autoregressive parameter.
#' @param cl True within-person cross-lag.
#' @param rho_i Trait correlation.
#' @param reliability Indicator reliability.
#' @param interval Length-2 numeric bound on acceptable \code{bw_ratio}
#'   solutions.
#' @return Scalar \code{bw_ratio} (the tipping point), or \code{NA_real_}.
#' @keywords internal
clpm_artifact_tipping <- function(target, ar = 0.30, cl = 0, rho_i = 0.50,
                                  reliability = 1, interval = c(0, 50)) {
  psi <- .psi_stationary_symmetric(ar, cl)
  rho_w <- psi$rho_w
  b <- ar
  cc <- cl
  p <- rho_i
  w <- reliability

  # a(r) = r + a0, g(r) = p*r + g0, v(r) = r/w + v0, h(r) = p*r + rho_w
  a0 <- b + cc * rho_w
  g0 <- b * rho_w + cc
  v0 <- 1 / w

  # cl_pop(r) - target == 0  <=>  [g(r)v(r) - a(r)h(r)] - target*[v(r)^2 - h(r)^2] == 0,
  # a quadratic A r^2 + B r + C in r.
  a_coef <- (p / w - p) - target * (1 / w^2 - p^2)
  b_coef <- (p * v0 + g0 / w - rho_w - p * a0) -
    target * (2 * v0 / w - 2 * p * rho_w)
  c_coef <- (g0 * v0 - a0 * rho_w) - target * (v0^2 - rho_w^2)

  roots <- if (abs(a_coef) < 1e-12) {
    if (abs(b_coef) < 1e-12) numeric(0) else -c_coef / b_coef
  } else {
    disc <- b_coef^2 - 4 * a_coef * c_coef
    if (disc < 0) {
      return(NA_real_)
    }
    (-b_coef + c(1, -1) * sqrt(disc)) / (2 * a_coef)
  }

  roots <- roots[
    is.finite(roots) & roots >= interval[1] - 1e-9 & roots <= interval[2] + 1e-9
  ]
  if (length(roots) == 0) {
    return(NA_real_)
  }
  max(min(roots), 0)
}

#' "As-if" CLPM artifact: population CLPM estimate under an RI-CLPM truth
#'
#' Holds true within-person dynamics (\code{ar}, \code{cl}) fixed and sweeps
#' the trait-variance share (\code{bw_ratio}) and optional indicator
#' reliability, returning the population CLPM cross-lag and autoregression a
#' CLPM analyst would obtain. The framing is "as-if": how much of an observed
#' CLPM cross-lag could be artifact produced by trait variance alone, given a
#' hypothesised true within-person cross-lag?
#'
#' The closed form generalises \code{clpmBias()} to a non-zero true cross-lag
#' and optional indicator measurement error; the stationary within covariance
#' is derived from the symmetric eigen-decomposition of
#' \eqn{B = [[ar, cl], [cl, ar]]}, requiring \eqn{|ar + cl| < 1} and
#' \eqn{|ar - cl| < 1} for stationarity.
#'
#' @param bw_ratio Between-to-within trait variance ratio. Vectorised.
#' @param ar True autoregressive parameter. Vectorised.
#' @param rho_i Trait correlation. Vectorised.
#' @param reliability Indicator reliability \eqn{\omega \in (0, 1]}. Vectorised.
#' @param expand If \code{TRUE}, evaluate on the full crossing of the
#'   supplied inputs; if \code{FALSE} (default), recycle them.
#' @param cl True within-person cross-lag. Vectorised.
#' @param verbose If \code{TRUE}, emits a one-line \code{message()} echoing the
#'   user-supplied inputs (the "generating syntax") so a reader of a knitted
#'   chunk can see which scenario produced the table. Defaults to \code{FALSE}.
#' @return Data frame with columns \code{bw_ratio}, \code{ar}, \code{cl},
#'   \code{rho_i}, \code{reliability}, \code{icc}, \code{rho_w}, \code{ar_pop},
#'   \code{cl_pop}, \code{ar_bias}, \code{ar_bias_pct}, \code{cl_bias},
#'   \code{cl_bias_pct}. Percentage columns are bias relative to the assumed
#'   truth (\eqn{100 \times \text{bias}/\text{truth}}); when a truth is
#'   zero, its percentage column is \code{NA} (division by zero), same
#'   convention for both \code{ar_bias_pct} and \code{cl_bias_pct}.
#' @examples
#' # Population CLPM cross-lag at increasing trait shares, true cl = 0
#' clpmCLBias(bw_ratio = c(0, 1, 2, 4, 8), ar = 0.3, cl = 0, rho_i = 0.5)
#'
#' # Compare across hypothesised true cross-lags
#' clpmCLBias(
#'   bw_ratio = c(0, 2, 4), ar = 0.3,
#'   cl = c(0, 0.05, 0.10), rho_i = 0.5, expand = TRUE
#' )
#' @export
clpmCLBias <- function(bw_ratio = 1, ar = 0.30, cl = 0, rho_i = 0.50,
                       reliability = 1, expand = FALSE, verbose = FALSE) {
  stopifnot(
    is.numeric(bw_ratio), is.numeric(ar), is.numeric(cl),
    is.numeric(rho_i), is.numeric(reliability),
    is.logical(expand), is.logical(verbose)
  )
  if (any(bw_ratio < 0)) stop("`bw_ratio` must be non-negative.", call. = FALSE)
  if (any(abs(rho_i) > 1)) stop("`rho_i` must lie in [-1, 1].", call. = FALSE)
  if (any(reliability <= 0 | reliability > 1)) {
    stop("`reliability` must lie in (0, 1].", call. = FALSE)
  }

  if (isTRUE(verbose)) {
    fmt <- function(x) {
      if (length(x) == 1) {
        format(x, trim = TRUE)
      } else {
        paste0("(", paste(format(x, trim = TRUE), collapse = ", "), ")")
      }
    }
    message(sprintf(
      paste0(
        "You specified a between-to-within ratio of %s, ",
        "an autoregressive parameter of %s, ",
        "a cross-lag parameter of %s, ",
        "a trait correlation of %s, ",
        "and a reliability of %s."
      ),
      fmt(bw_ratio), fmt(ar), fmt(cl), fmt(rho_i), fmt(reliability)
    ))
  }

  if (isTRUE(expand)) {
    g <- expand.grid(
      bw_ratio = bw_ratio, ar = ar, cl = cl,
      rho_i = rho_i, reliability = reliability
    )
    bw_ratio <- g$bw_ratio
    ar <- g$ar
    cl <- g$cl
    rho_i <- g$rho_i
    reliability <- g$reliability
  }

  out <- clpm_artifact_grid(bw_ratio,
    ar = ar, cl = cl,
    rho_i = rho_i, reliability = reliability
  )
  out$icc <- out$bw_ratio / (1 + out$bw_ratio)
  out[c(
    "bw_ratio", "ar", "cl", "rho_i", "reliability", "icc",
    "rho_w", "ar_pop", "cl_pop",
    "ar_bias", "ar_bias_pct",
    "cl_bias", "cl_bias_pct"
  )]
}

#' Tipping-point trait share for an observed CLPM cross-lag
#'
#' Convenience wrapper around \code{clpm_artifact_tipping()}. For each
#' hypothesised true within-person cross-lag in \code{cl}, returns the
#' between-to-within variance ratio at which the population CLPM cross-lag
#' would equal \code{observed}. Useful for "as-if" sensitivity statements
#' (e.g., "to attribute the observed cross-lag to trait variance alone we
#' would need r = ..., implying ICC = ...").
#'
#' @param observed Numeric scalar, the observed CLPM cross-lag.
#' @param ar True autoregressive parameter (scalar).
#' @param cl Numeric vector of hypothesised true within-person cross-lags.
#' @param rho_i Trait correlation (scalar).
#' @param reliability Indicator reliability (scalar).
#' @param interval Search interval for \code{bw_ratio}.
#' @return Data frame with one row per \code{cl}: the tipping-point
#'   \code{bw_ratio} and implied trait \code{icc}.
#' @examples
#' clpmArtifactTipping(
#'   observed = 0.08, ar = 0.3, cl = c(0, 0.02, 0.05),
#'   rho_i = 0.5
#' )
#' @export
clpmArtifactTipping <- function(observed, ar = 0.30, cl = 0, rho_i = 0.50,
                                reliability = 1, interval = c(0, 50)) {
  stopifnot(
    length(observed) == 1L, length(ar) == 1L, length(rho_i) == 1L,
    length(reliability) == 1L
  )
  r_star <- vapply(cl, function(cc) {
    clpm_artifact_tipping(
      target = observed, ar = ar, cl = cc,
      rho_i = rho_i, reliability = reliability,
      interval = interval
    )
  }, numeric(1))
  data.frame(
    cl       = cl,
    bw_ratio = r_star,
    icc      = r_star / (1 + r_star)
  )
}

