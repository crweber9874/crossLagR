#' @title generate_riclpm_indicator_lists_for_estimation
#' @description Generates lists of indicator variable names in the format required by the estimateRICLPMl function.
#'
#' @param waves The number of waves (time points) in the model.
#' @param num_indicators The number of indicators per construct (X, Y, and optionally Z) at each wave.
#' @param has_z Logical, indicating whether the model includes a Z variable. Defaults to FALSE.
#' @param x_prefix The prefix for the X variable indicators. Defaults to "x".
#' @param y_prefix The prefix for the Y variable indicators. Defaults to "y".
#' @param z_prefix The prefix for the Z variable indicators. Defaults to "z".
#' @param indicator_separator The separator between the prefix, indicator number, and wave number. Defaults to "_".
#'
#' @return A list containing the `time_varying_x`, `time_varying_y`, and optionally `time_varying_z` lists, formatted for `estimateRICLPMl`.
#'
#' @examples
#' indicator_lists_2waves_3indicators_for_est <- generate_riclpm_indicator_lists_for_estimation(waves = 2, num_indicators = 3)
#' print(indicator_lists_2waves_3indicators_for_est)
#'
#' indicator_lists_3waves_2indicators_withZ_for_est <- generate_riclpm_indicator_lists_for_estimation(waves = 3, num_indicators = 2, has_z = TRUE)
#' print(indicator_lists_3waves_2indicators_withZ_for_est)
#'
#' indicator_lists_2waves_2indicators_custom_prefix_for_est <- generate_riclpm_indicator_lists_for_estimation(waves = 2, num_indicators = 2, x_prefix = "varX", y_prefix = "varY")
#' print(indicator_lists_2waves_2indicators_custom_prefix_for_est)
#'
#' @export
generate_riclpm_indicator_lists_for_estimation <- function(waves,
                                                           num_indicators,
                                                           has_z = FALSE,
                                                           x_prefix = "x",
                                                           y_prefix = "y",
                                                           z_prefix = "z",
                                                           indicator_separator = "_") {
  if (!is.numeric(waves) || waves <= 0 || waves != as.integer(waves)) {
    stop("❌ Error: 'waves' must be a positive integer.")
  }
  if (!is.numeric(num_indicators) || num_indicators <= 0 || num_indicators != as.integer(num_indicators)) {
    stop("❌ Error: 'num_indicators' must be a positive integer.")
  }
  if (!is.logical(has_z)) {
    stop("❌ Error: 'has_z' must be a logical value (TRUE or FALSE).")
  }
  if (!is.character(x_prefix) || length(x_prefix) != 1) {
    stop("❌ Error: 'x_prefix' must be a single character string.")
  }
  if (!is.character(y_prefix) || length(y_prefix) != 1) {
    stop("❌ Error: 'y_prefix' must be a single character string.")
  }
  if (!is.character(z_prefix) || length(z_prefix) != 1) {
    stop("❌ Error: 'z_prefix' must be a single character string.")
  }
  if (!is.character(indicator_separator) || length(indicator_separator) != 1) {
    stop("❌ Error: 'indicator_separator' must be a single character string.")
  }

  time_varying_x <- lapply(1:num_indicators, function(i) {
    paste0(x_prefix, i, indicator_separator, 1:waves)
  })

  time_varying_y <- lapply(1:num_indicators, function(i) {
    paste0(y_prefix, i, indicator_separator, 1:waves)
  })

  result_list <- list(
    time_varying_x = time_varying_x,
    time_varying_y = time_varying_y
  )

  if (has_z) {
    time_varying_z <- lapply(1:num_indicators, function(i) {
      paste0(z_prefix, i, indicator_separator, 1:waves)
    })
    result_list$time_varying_z <- time_varying_z
  }

  return(result_list)
}
