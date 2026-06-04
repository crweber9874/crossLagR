#' @title iccReport
#' @description Compute the intraclass correlation (ICC) and between-/within-
#'   person variance decomposition for panel data, in either wide or long
#'   format. Designed to be run alongside a lavaan fit so the user sees the
#'   data-level ICC structure before interpreting model output.
#'
#' @param data A data.frame.
#' @param vars Character vector. In \code{format = "wide"}, each entry is a
#'   column-name stem (e.g., \code{"x"}) matched against columns \code{x1},
#'   \code{x2}, ..., \code{xT}. In \code{format = "long"}, each entry is a
#'   variable column name in the long data.
#' @param format Either \code{"wide"} (default) or \code{"long"}.
#' @param id Character. Long format only: name of the unit-ID column.
#' @param time Character. Long format only: name of the time/wave column.
#' @param threshold Numeric. ICC values above this are flagged in
#'   \code{exceeds_threshold} and trigger a warning. Default \code{0.65}.
#' @param quiet Logical. If \code{FALSE} (default), emits a single condensed
#'   \code{warning()} when any ICC exceeds the threshold, suggesting
#'   sensitivity analysis with \code{\link{launch_shiny_sim}}. Set \code{TRUE}
#'   for silent programmatic use.
#'
#' @return A data.frame of class \code{crossLagR_icc} with one row per
#'   variable and columns:
#'   \itemize{
#'     \item \code{variable}, \code{between_var}, \code{within_var},
#'       \code{total_var}, \code{icc}
#'     \item \code{n_obs} (complete cases used), \code{n_waves}
#'     \item \code{exceeds_threshold} (logical)
#'   }
#'   The threshold is stored as \code{attr(x, "threshold")}.
#'
#' @details
#' Computed empirically as
#' \deqn{ICC = \frac{Var(\bar{x}_i)}{Var(\bar{x}_i) + Var(x_{it} - \bar{x}_i)}}
#' where \eqn{\bar{x}_i} is the person mean. This matches the
#' \code{method = "dplyr"} decomposition in \code{\link{withinBetween}}.
#'
#' Why the threshold matters: the CLPM does not separate stable between-person
#' differences from within-person dynamics, so when between-person variance
#' dominates (high ICC), CLPM and RICLPM estimates can diverge substantially
#' (Hamaker, Kuiper & Grasman, 2015). The 0.65 default is a rough heuristic
#' marking the region where the modeling choice begins to matter materially;
#' there is no formal cutoff in the literature, so adjust \code{threshold}
#' for your application.
#'
#' @references
#' Hamaker, E. L., Kuiper, R. M., & Grasman, R. P. P. P. (2015). A critique
#' of the cross-lagged panel model. \emph{Psychological Methods}, 20(1),
#' 102--116.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' dat <- simRICLPM(waves = 4, sample.nobs = 300,
#'                  var_BX = 2.0, var_BY = 2.0)$data
#' iccReport(dat, vars = c("x", "y"))
#' }
#'
#' @seealso \code{\link{withinBetween}}, \code{\link{launch_shiny_sim}},
#'   \code{\link{lavaanDiagnose}}.
#' @export
iccReport <- function(data,
                      vars,
                      format = c("wide", "long"),
                      id = NULL,
                      time = NULL,
                      threshold = 0.65,
                      quiet = FALSE) {

  format <- match.arg(format)
  if (!is.data.frame(data)) stop("'data' must be a data.frame.")
  if (!is.character(vars) || !length(vars)) {
    stop("'vars' must be a non-empty character vector.")
  }
  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold)) {
    stop("'threshold' must be a single numeric value.")
  }

  out <- if (format == "wide") {
    .icc_from_wide(data, vars)
  } else {
    if (is.null(id) || is.null(time)) {
      stop("format = 'long' requires both 'id' and 'time' arguments.")
    }
    .icc_from_long(data, id = id, time = time, vars = vars)
  }

  out$exceeds_threshold <- !is.na(out$icc) & out$icc > threshold

  if (!quiet && any(out$exceeds_threshold)) {
    flagged <- out[out$exceeds_threshold, , drop = FALSE]
    detail <- paste(
      sprintf("%s (ICC = %.2f)", flagged$variable, flagged$icc),
      collapse = "; "
    )
    warning(
      "crossLagR: ICC exceeds ", threshold, " for: ", detail, ". ",
      "When between-person variance dominates, CLPM and RICLPM estimates ",
      "can diverge materially (Hamaker, Kuiper & Grasman, 2015). ",
      "Consider exploring sensitivity with the bundled simulator: ",
      "launch_shiny_sim().",
      call. = FALSE
    )
  }

  class(out) <- c("crossLagR_icc", "data.frame")
  attr(out, "threshold") <- threshold
  out
}


#' @export
print.crossLagR_icc <- function(x, digits = 3, ...) {
  thr <- attr(x, "threshold")
  cat("crossLagR ICC report  (threshold = ", thr, ")\n", sep = "")
  cat(strrep("-", 41), "\n", sep = "")
  body <- as.data.frame(x)
  num <- vapply(body, is.numeric, logical(1))
  body[num] <- lapply(body[num], function(v) {
    if (is.integer(v)) v else round(v, digits)
  })
  print(body, row.names = FALSE)
  if (any(x$exceeds_threshold)) {
    flagged <- x$variable[x$exceeds_threshold]
    cat("\n[!] ICC > ", thr, " for: ",
        paste(flagged, collapse = ", "), "\n", sep = "")
    cat("    High ICC means between-person variance dominates; CLPM vs RICLPM\n",
        "    estimates can diverge materially (Hamaker et al., 2015).\n",
        "    Explore sensitivity with launch_shiny_sim().\n", sep = "")
  }
  invisible(x)
}


## ----------------------------------------------------------------------------
## Internal computation
## ----------------------------------------------------------------------------

.icc_from_wide <- function(data, vars) {
  rows <- lapply(vars, function(stem) {
    pat <- paste0("^", stem, "[0-9]+$")
    cols <- grep(pat, names(data), value = TRUE)
    if (!length(cols)) {
      warning("iccReport: no columns matching '", stem,
              "[0-9]+' in data.", call. = FALSE)
      return(.icc_empty_row(stem, 0L))
    }
    wave_num <- as.integer(sub(stem, "", cols))
    cols <- cols[order(wave_num)]
    mat  <- as.matrix(data[, cols, drop = FALSE])
    keep <- stats::complete.cases(mat)
    mat  <- mat[keep, , drop = FALSE]
    if (nrow(mat) < 2L) return(.icc_empty_row(stem, length(cols), nrow(mat)))

    person_mean <- rowMeans(mat)
    between_var <- stats::var(person_mean)
    within_var  <- stats::var(as.vector(mat - person_mean))
    total       <- between_var + within_var
    data.frame(
      variable    = stem,
      between_var = between_var,
      within_var  = within_var,
      total_var   = total,
      icc         = if (total > 0) between_var / total else NA_real_,
      n_obs       = nrow(mat),
      n_waves     = length(cols),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.icc_from_long <- function(data, id, time, vars) {
  wb <- withinBetween(data, id = id, time = time, vars = vars,
                      method = "dplyr")
  data.frame(
    variable    = wb$variable,
    between_var = wb$between_var,
    within_var  = wb$within_var,
    total_var   = wb$total_var,
    icc         = wb$icc,
    n_obs       = vapply(wb$variable, function(v) {
      sum(stats::complete.cases(data[, c(id, time, v)]))
    }, integer(1)),
    n_waves     = length(unique(data[[time]])),
    stringsAsFactors = FALSE
  )
}

.icc_empty_row <- function(stem, n_waves, n_obs = 0L) {
  data.frame(
    variable    = stem,
    between_var = NA_real_,
    within_var  = NA_real_,
    total_var   = NA_real_,
    icc         = NA_real_,
    n_obs       = as.integer(n_obs),
    n_waves     = as.integer(n_waves),
    stringsAsFactors = FALSE
  )
}
