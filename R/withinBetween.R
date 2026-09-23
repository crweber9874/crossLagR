#' @title withinBetween
#' @description Decompose panel variables into between-person and within-person
#'   variance components.
#'
#' @param data A long-format panel data frame.
#' @param id Character. Name of the unit (person) identifier column.
#' @param time Character. Name of the time/wave column.
#' @param vars Character vector of variable names to decompose.
#' @param method Either \code{"dplyr"} (default) for an empirical decomposition
#'   via group-mean centering, or \code{"tso"} to fit a bivariate Trait-State-Occasion
#'   model with lavaan and read the structural variance estimates. The
#'   \code{"tso"} option requires exactly two variables.
#'
#' @return A data frame with one row per variable and columns
#'   \code{variable}, \code{between_var}, \code{within_var}, \code{total_var},
#'   \code{icc}, \code{method}.
#'
#' @details
#' The \code{"dplyr"} method estimates between-person variance as the variance
#' of person-specific means and within-person variance as the variance of
#' deviations from those means. The \code{"tso"} method fits the bivariate
#' Trait-State-Occasion model from \code{\link{estimateTSO}} and reads
#' \code{T_var_x}, \code{T_var_y} (between) and \code{s_var_x}, \code{s_var_y}
#' (within) from the lavaan parameter table.
#'
#' @examples
#' \dontrun{
#' withinBetween(upenn_long, id = "identifier", time = "wave",
#'               vars = c("party_identification", "political_legitimacy"))
#'
#' withinBetween(upenn_long, id = "identifier", time = "wave",
#'               vars = c("party_identification", "political_legitimacy"),
#'               method = "tso")
#' }
#' @export
withinBetween <- function(data, id, time, vars, method = c("dplyr", "tso")) {
  method <- match.arg(method)
  stopifnot(is.character(vars), length(vars) >= 1)
  stopifnot(is.character(id), is.character(time))

  if (method == "dplyr") {
    rows <- lapply(vars, function(v) {
      d <- data[, c(id, time, v)]
      d <- d[stats::complete.cases(d), , drop = FALSE]
      person_means <- tapply(d[[v]], d[[id]], mean, na.rm = TRUE)
      d$mean_i <- person_means[as.character(d[[id]])]
      d$within <- d[[v]] - d$mean_i
      between_var <- stats::var(person_means, na.rm = TRUE)
      within_var  <- stats::var(d$within, na.rm = TRUE)
      data.frame(
        variable    = v,
        between_var = between_var,
        within_var  = within_var,
        total_var   = between_var + within_var,
        icc         = between_var / (between_var + within_var),
        method      = "dplyr",
        stringsAsFactors = FALSE
      )
    })
    return(do.call(rbind, rows))
  }

  # method == "tso"
  if (length(vars) != 2L) {
    stop("method = 'tso' requires exactly two variables (bivariate TSO).")
  }
  waves <- sort(unique(data[[time]]))
  n_waves <- length(waves)
  if (n_waves < 3) stop("TSO requires at least 3 waves.")

  wide <- unique(data[, id, drop = FALSE])
  for (i in seq_along(vars)) {
    prefix <- c("x", "y")[i]
    sub <- data[, c(id, time, vars[i])]
    sub$.wave_idx <- as.integer(factor(sub[[time]], levels = waves))
    w <- stats::reshape(
      sub[, c(id, ".wave_idx", vars[i])],
      idvar     = id,
      timevar   = ".wave_idx",
      direction = "wide",
      sep       = ""
    )
    names(w) <- sub(paste0("^", vars[i]), prefix, names(w))
    wide <- merge(wide, w, by = id, all.x = TRUE)
  }

  syntax <- estimateTSO(waves = n_waves)
  fit <- lavaan::lavaan(syntax, data = wide, missing = "ML",
                        meanstructure = TRUE)
  pe <- lavaan::parameterEstimates(fit)
  grab <- function(lab) pe$est[pe$label == lab][1]

  between <- c(grab("T_var_x"), grab("T_var_y"))
  within  <- c(grab("s_var_x"), grab("s_var_y"))
  data.frame(
    variable    = vars,
    between_var = between,
    within_var  = within,
    total_var   = between + within,
    icc         = between / (between + within),
    method      = "tso",
    stringsAsFactors = FALSE
  )
}


#' @title simWithinBetween
#' @description Parametric bootstrap of between- and within-person variance
#'   components implied by a \code{\link{withinBetween}} fit.
#'
#' @param wb A data frame returned by \code{withinBetween()}.
#' @param n_units Integer. Number of units (persons) per simulated panel.
#'   Defaults to a sensible value matching typical panels.
#' @param n_waves Integer. Number of waves per simulated panel.
#' @param n_sims Integer. Number of bootstrap replicates. Default 500.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A data frame with columns \code{variable}, \code{sim},
#'   \code{between_var}, \code{within_var}.
#'
#' @details
#' For each variable, simulates \code{n_sims} panels from
#' \eqn{y_{it} = \alpha_i + \varepsilon_{it}} with
#' \eqn{\alpha_i \sim N(0, \sigma^2_B)} and
#' \eqn{\varepsilon_{it} \sim N(0, \sigma^2_W)}, then recomputes the empirical
#' between- and within-person variance components in each replicate. The
#' resulting distribution describes the sampling variability of the
#' decomposition under the implied model.
#'
#' @export
simWithinBetween <- function(wb, n_units, n_waves, n_sims = 500, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (missing(n_units) || missing(n_waves)) {
    stop("Both 'n_units' and 'n_waves' must be supplied.")
  }

  rows <- lapply(seq_len(nrow(wb)), function(i) {
    bv <- wb$between_var[i]
    wv <- wb$within_var[i]
    v  <- wb$variable[i]
    sims <- lapply(seq_len(n_sims), function(s) {
      alpha <- stats::rnorm(n_units, 0, sqrt(max(bv, 0)))
      eps   <- matrix(stats::rnorm(n_units * n_waves, 0, sqrt(max(wv, 0))),
                      n_units, n_waves)
      y <- alpha + eps
      pm <- rowMeans(y)
      data.frame(
        variable    = v,
        sim         = s,
        between_var = stats::var(pm),
        within_var  = stats::var(as.vector(y - pm)),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, sims)
  })
  do.call(rbind, rows)
}
