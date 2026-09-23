# Inverse of the CLPM bias map: given the AR and CL a CLPM analyst reports,
# recover the RI-CLPM parameters consistent with them under a true-CL-of-zero
# null. Companion to clpmBias() / clpmCLBias(), which run the map forwards.

#' Recover RI-CLPM parameters implied by reported CLPM estimates
#'
#' A published CLPM reports an autoregression and a cross-lag but rarely an
#' ICC, so the between-to-within variance ratio \eqn{r} is unknown. This
#' function asks the inverse question: if the true within-person cross-lag were
#' zero, what \eqn{(b, r)} would reproduce the reported estimates?
#'
#' The forward map has three unknowns (\eqn{b}, \eqn{r}, \eqn{\rho_I}) and the
#' analyst supplies only two numbers, so the system is underdetermined.
#' Fixing \eqn{\rho_I} identifies \eqn{(b, r)} exactly, so the function sweeps
#' a grid of assumed trait correlations and returns one solution per value --
#' a sensitivity curve, not a point estimate.
#'
#' @details
#' The solution is closed form. Under the symmetric parameterisation the trait
#' covariance \eqn{\boldsymbol\Omega = r[[1, \rho_I], [\rho_I, 1]]} has fixed
#' eigenvectors \eqn{(1,1)} and \eqn{(1,-1)} whatever \eqn{r} and \eqn{\rho_I}
#' are, and \eqn{\mathbf B^{\text{true}} = b\mathbf I} is diagonal in every
#' basis. The whole system therefore diagonalises in one known basis and the
#' two equations decouple. Writing \eqn{\lambda_1 = A + C} and
#' \eqn{\lambda_2 = A - C} for the eigenvalues of the reported matrix, each
#' satisfies
#' \deqn{\lambda_i = b + (1-b)\,\omega_i/(\omega_i + 1)
#'       \quad\Longrightarrow\quad \omega_i = (\lambda_i - b)/(1 - \lambda_i),}
#' with \eqn{\omega_1 = r(1+\rho_I)} and \eqn{\omega_2 = r(1-\rho_I)}.
#' Requiring both to imply the same \eqn{r} is linear in \eqn{b}, so no root
#' finding is needed. The algebraic solution is unconstrained, so solutions
#' outside the admissible region (\eqn{|b| \ge 1} or \eqn{r \le 0}) are
#' returned as \code{NA}.
#'
#' That shortcut requires a symmetric reported matrix. When the reported
#' autoregressions or cross-lags differ across series, \eqn{(\mathbf I -
#' \mathbf B)} is no longer a scalar multiple of the identity, the
#' eigenvectors rotate away from \eqn{(1, \pm 1)}, and the system does not
#' decouple. The asymmetric case has four unknowns (\eqn{b_x}, \eqn{b_y},
#' \eqn{r_x}, \eqn{r_y}) and four reported entries, so it remains exactly
#' determined given \eqn{\rho_I}, but it is solved numerically: for a
#' candidate \eqn{(r_x, r_y)} each series yields two expressions for its true
#' autoregression, and the solver drives their difference to zero.
#'
#' @param ar_obs Reported CLPM autoregression. Length 1 for the symmetric
#'   case, or length 2 as \code{c(ar_x, ar_y)} for the asymmetric case.
#' @param cl_obs Reported CLPM cross-lag. Length 1 for the symmetric case, or
#'   length 2 as \code{c(cl_xy, cl_yx)} for the asymmetric case. Must be
#'   non-zero.
#' @param rho_i Assumed trait correlation(s). Vectorised; values at exactly
#'   \eqn{\pm 1} are invalid.
#' @param b_max Largest admissible absolute value for the true autoregression.
#'   Solutions with \eqn{|b| \ge} \code{b_max} are reported as \code{NA}.
#' @return A data frame with one row per \code{rho_i}. For symmetric input the
#'   columns are \code{rho_i}, \code{ar_implied}, \code{bw_ratio_implied}
#'   and \code{icc_implied} (\eqn{r/(1+r)}). For asymmetric input the
#'   per-series columns \code{ar_implied_x}, \code{ar_implied_y},
#'   \code{bw_ratio_implied_x}, \code{bw_ratio_implied_y},
#'   \code{icc_implied_x} and \code{icc_implied_y} are returned instead.
#'   Rows with no admissible solution return \code{NA}. The \emph{implied}
#'   naming is deliberate: these are not estimates of the truth, but what the
#'   truth would have to be for the reported values to arise under the assumed
#'   \code{rho_i} and a true cross-lag of zero.
#' @examples
#' # An observed AR of 0.69 and CL of 0.09, swept over assumed rho_I
#' clpmInvert(ar_obs = 0.692, cl_obs = 0.092, rho_i = c(0.3, 0.5, 0.7))
#' @export
clpmInvert <- function(ar_obs, cl_obs, rho_i = seq(0.1, 0.9, by = 0.1),
                       b_max = 1, tol = 1e-8) {
  stopifnot(
    is.numeric(ar_obs), is.numeric(cl_obs), is.numeric(rho_i),
    length(ar_obs) %in% 1:2, length(cl_obs) %in% 1:2, is.numeric(b_max)
  )
  if (any(cl_obs == 0)) stop("`cl_obs` must be non-zero.", call. = FALSE)
  if (any(abs(rho_i) >= 1)) stop("`rho_i` must lie strictly within (-1, 1).", call. = FALSE)

  ar_obs <- rep_len(ar_obs, 2)
  cl_obs <- rep_len(cl_obs, 2)
  symmetric <- isTRUE(all.equal(ar_obs[1], ar_obs[2], tolerance = tol)) &&
    isTRUE(all.equal(cl_obs[1], cl_obs[2], tolerance = tol))

  if (!symmetric) {
    return(.clpm_invert_asym(ar_obs, cl_obs, rho_i, b_max))
  }
  ar_obs <- ar_obs[1]
  cl_obs <- cl_obs[1]

  # Eigenvalues of the reported (symmetric) coefficient matrix.
  lambda1 <- ar_obs + cl_obs
  lambda2 <- ar_obs - cl_obs

  # Requiring both eigen-equations to imply the same r is linear in b.
  a1 <- (1 - lambda2) * (1 - rho_i)
  a2 <- (1 - lambda1) * (1 + rho_i)
  b  <- (lambda2 * a2 - lambda1 * a1) / (a2 - a1)

  omega1 <- (lambda1 - b) / (1 - lambda1)
  omega2 <- (lambda2 - b) / (1 - lambda2)
  r      <- omega1 / (1 + rho_i)

  # Reject the inadmissible algebraic solutions.
  consistent <- is.finite(omega2) & is.finite(r) &
    abs(omega2 / (1 - rho_i) - r) < 1e-8
  ok <- is.finite(b) & is.finite(r) & consistent & r > 0 & abs(b) < b_max
  b[!ok] <- NA_real_
  r[!ok] <- NA_real_

  data.frame(
    rho_i            = rho_i,
    ar_implied       = b,
    bw_ratio_implied = r,
    icc_implied      = r / (1 + r)
  )
}

# Asymmetric reported matrix: B_hat = [[ar_x, cl_xy], [cl_yx, ar_y]] with
# Omega = [[r_x, rho*sqrt(r_x r_y)], [rho*sqrt(r_x r_y), r_y]] and Psi = I.
# Four unknowns, four reported entries, solved numerically per rho_i.
.clpm_invert_asym <- function(ar_obs, cl_obs, rho_i, b_max) {
  V_of <- function(rx, ry, p) {
    Om <- matrix(c(rx, p * sqrt(rx * ry), p * sqrt(rx * ry), ry), 2)
    Om %*% solve(Om + diag(2))
  }
  # For a candidate (r_x, r_y) each series gives two expressions for its true
  # autoregression -- one from the AR entry, one from the CL entry. At the
  # solution they agree.
  resid <- function(par, p) {
    rx <- exp(par[1])
    ry <- exp(par[2])
    V <- tryCatch(V_of(rx, ry, p), error = function(e) NULL)
    if (is.null(V) || !all(is.finite(V))) return(c(1e6, 1e6))
    c((1 - cl_obs[1] / V[1, 2]) - (ar_obs[1] - V[1, 1]) / (1 - V[1, 1]),
      (1 - cl_obs[2] / V[2, 1]) - (ar_obs[2] - V[2, 2]) / (1 - V[2, 2]))
  }

  solve_one <- function(p) {
    best <- NULL
    for (st in list(c(0, 0), c(1, 1), c(-1, -1), c(1, -1), c(-1, 1))) {
      o <- tryCatch(
        stats::optim(st, function(z) sum(resid(z, p)^2), method = "BFGS",
                     control = list(reltol = 1e-14, maxit = 500)),
        error = function(e) NULL
      )
      if (!is.null(o) && (is.null(best) || o$value < best$value)) best <- o
    }
    if (is.null(best) || best$value > 1e-10) return(rep(NA_real_, 4))
    rx <- exp(best$par[1])
    ry <- exp(best$par[2])
    V <- V_of(rx, ry, p)
    bx <- 1 - cl_obs[1] / V[1, 2]
    by <- 1 - cl_obs[2] / V[2, 1]
    ok <- all(is.finite(c(rx, ry, bx, by))) && rx > 0 && ry > 0 &&
      abs(bx) < b_max && abs(by) < b_max
    if (!ok) return(rep(NA_real_, 4))
    c(bx, by, rx, ry)
  }

  sol <- vapply(rho_i, solve_one, numeric(4))
  data.frame(
    rho_i              = rho_i,
    ar_implied_x       = sol[1, ],
    ar_implied_y       = sol[2, ],
    bw_ratio_implied_x = sol[3, ],
    bw_ratio_implied_y = sol[4, ],
    icc_implied_x      = sol[3, ] / (1 + sol[3, ]),
    icc_implied_y      = sol[4, ] / (1 + sol[4, ])
  )
}
