#' @title estimateFI
#' @description Estimate Fixed Individual Effects model using within-person deviations
#'
#' This function fits two OLS regression models using within-person deviations to control
#' for time-invariant unmeasured confounders. This approach is equivalent to including
#' individual fixed effects in a panel data model.
#'
#' @param data Long-format data frame with columns for within-person deviations and lagged variables.
#'   Must include: within.x, within.y, xlagw, ylagw (created by reshape_long_sim_cr function).
#' @param y_outcome Character string specifying the Y outcome variable name. Default is "within.y".
#' @param x_outcome Character string specifying the X outcome variable name. Default is "within.x".
#' @param y_lag Character string specifying the lagged Y predictor name. Default is "ylagw".
#' @param x_lag Character string specifying the lagged X predictor name. Default is "xlagw".
#' @param return_models Logical. If TRUE, returns the full lm objects. If FALSE, returns
#'   just the coefficients. Default is FALSE.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{xlag_x}: Autoregressive effect of X (X(t-1) -> X(t))
#'     \item \code{ylag_x}: Cross-lagged effect of Y on X (Y(t-1) -> X(t))
#'     \item \code{xlag_y}: Cross-lagged effect of X on Y (X(t-1) -> Y(t))
#'     \item \code{ylag_y}: Autoregressive effect of Y (Y(t-1) -> Y(t))
#'     \item \code{converged}: Always TRUE for OLS (included for consistency)
#'     \item \code{models}: If return_models=TRUE, includes the fitted lm objects
#'   }
#'
#' @details
#' Fixed Individual Effects estimation works by:
#' \enumerate{
#'   \item Computing within-person deviations (person-mean centering)
#'   \item Fitting OLS regressions on these deviations with lagged within-person deviations
#'   \item This removes all time-invariant between-person confounders
#' }
#'
#' The approach is particularly effective for controlling time-invariant confounders,
#' but may still be biased by time-varying confounders.
#'
#' @examples
#' \dontrun{
#' # Generate some data
#' sim_data <- simRICLPM(waves = 4, sample.nobs = 500)$data
#' long_data <- reshape_long_sim_cr(sim_data)
#'
#' # Estimate Fixed Individual Effects model
#' fi_results <- estimateFI(long_data)
#' print(fi_results)
#'
#' # Get full model objects
#' fi_full <- estimateFI(long_data, return_models = TRUE)
#' summary(fi_full$models$x_model)
#' summary(fi_full$models$y_model)
#' }
#'
#' @import stats
#' @export
estimateFI <- function(data,
                       y_outcome = "within.y",
                       x_outcome = "within.x",
                       y_lag = "ylagw",
                       x_lag = "xlagw",
                       return_models = FALSE) {

  # Input validation
  required_cols <- c(y_outcome, x_outcome, y_lag, x_lag)
  missing_cols <- required_cols[!required_cols %in% names(data)]

  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", "),
               "\nMake sure data is from reshape_long_sim_cr() function."))
  }

  # Remove rows with missing values for the analysis
  analysis_data <- data[complete.cases(data[required_cols]), ]

  if (nrow(analysis_data) == 0) {
    stop("No complete cases available for analysis after removing missing values.")
  }

  # Check for sufficient variation
  if (var(analysis_data[[y_lag]], na.rm = TRUE) == 0 ||
      var(analysis_data[[x_lag]], na.rm = TRUE) == 0) {
    warning("Zero variance detected in lagged predictors. Results may be unreliable.")
  }

  # Fit Fixed Individual Effects models
  tryCatch({
    # Y equation: within.y ~ xlagw + ylagw
    formula_y <- as.formula(paste(y_outcome, "~", x_lag, "+", y_lag))
    model_fit_y <- lm(formula_y, data = analysis_data)

    # X equation: within.x ~ xlagw + ylagw
    formula_x <- as.formula(paste(x_outcome, "~", x_lag, "+", y_lag))
    model_fit_x <- lm(formula_x, data = analysis_data)

    # Extract coefficients safely
    get_coef_safe <- function(model, coef_name) {
      coefs <- coefficients(model)
      if (coef_name %in% names(coefs) && !is.na(coefs[coef_name])) {
        return(as.numeric(coefs[coef_name]))
      } else {
        return(NA_real_)
      }
    }

    # Extract coefficients
    xlag_x <- get_coef_safe(model_fit_x, x_lag)   # Autoregressive X effect
    ylag_x <- get_coef_safe(model_fit_x, y_lag)   # Cross-lagged Y -> X effect
    xlag_y <- get_coef_safe(model_fit_y, x_lag)   # Cross-lagged X -> Y effect
    ylag_y <- get_coef_safe(model_fit_y, y_lag)   # Autoregressive Y effect

    # Prepare results
    results <- list(
      xlag_x = xlag_x,
      ylag_x = ylag_x,
      xlag_y = xlag_y,
      ylag_y = ylag_y,
      converged = TRUE,  # OLS always "converges"
      n_obs = nrow(analysis_data),
      method = "Fixed Individual Effects (Within-Person Deviations)"
    )

    # Add full models if requested
    if (return_models) {
      results$models <- list(
        x_model = model_fit_x,
        y_model = model_fit_y,
        x_formula = formula_x,
        y_formula = formula_y
      )
    }

    return(results)

  }, error = function(e) {
    # Return NA results if fitting fails
    results <- list(
      xlag_x = NA_real_,
      ylag_x = NA_real_,
      xlag_y = NA_real_,
      ylag_y = NA_real_,
      converged = FALSE,
      n_obs = nrow(analysis_data),
      method = "Fixed Individual Effects (Within-Person Deviations)",
      error_message = as.character(e$message)
    )

    if (return_models) {
      results$models <- NULL
    }

    warning(paste("Fixed Individual Effects estimation failed:", e$message))
    return(results)
  })
}
