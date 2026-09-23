# Small internal helper shared across the monteCarlo* runners. Previously
# copy-pasted verbatim into 4 separate files (monteCarlo_CLPM.R,
# monteCarlo_RICLPM.R, monteCarlo_LChange.R, monteCarloAllisonChamberlainFI.R);
# only monteCarlo_LChange.R actually called it, the other three copies were
# dead code. Kept the un-prefixed name (matching the original copies) so
# monteCarlo_LChange.R's existing call sites needed no changes.

#' Safely extract a named coefficient, returning NA if it isn't present
#' @param coef_name Character. Name of the coefficient to extract.
#' @param coeffs Named numeric vector of coefficients (e.g. from `coef()`).
#' @return Numeric scalar, or `NA_real_` if `coef_name` is not in `coeffs`.
#' @keywords internal
get_coef_safe <- function(coef_name, coeffs) {
  if (coef_name %in% names(coeffs)) {
    as.numeric(coeffs[coef_name])
  } else {
    NA_real_
  }
}
