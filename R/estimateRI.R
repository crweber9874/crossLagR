#' @title estimateRI
#' @description Estimate Random Intercepts model using lmer() from lme4 package
#'
#' This function fits a random intercepts model using lme4::lmer(), where individual-specific
#' intercepts are treated as random effects. This approach is intermediate between OLS
#' (which pools all individuals) and Fixed Individual Effects (which estimates separate
#' intercepts for each individual).
#'
#' @param data Long-format data frame with columns for variables and lagged variables.
#'   Must include: id, x, y, xlag, ylag (created by reshape_long_sim_cr function).
#' @param y_outcome Character string specifying the Y outcome variable name. Default is "y".
#' @param x_outcome Character string specifying the X outcome variable name. Default is "x".
#' @param y_lag Character string specifying the lagged Y predictor name. Default is "ylag".
#' @param x_lag Character string specifying the lagged X predictor name. Default is "xlag".
#' @param id_var Character string specifying the ID variable name. Default is "id".
#' @param return_models Logical. If TRUE, returns the full lmer objects. If FALSE, returns
#'   just the coefficients. Default is FALSE.
#' @param control_params List of control parameters to pass to lmer(). Default uses
#'   lmerControl with optimizer = "bobyqa" and optCtrl = list(maxfun = 20000).
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{xlag_x}: Autoregressive effect of X (X(t-1) -> X(t))
#'     \item \code{ylag_x}: Cross-lagged effect of Y on X (Y(t-1) -> X(t))
#'     \item \code{xlag_y}: Cross-lagged effect of X on Y (X(t-1) -> Y(t))
#'     \item \code{ylag_y}: Autoregressive effect of Y (Y(t-1) -> Y(t))
#'     \item \code{converged}: Whether both models converged successfully
#'     \item \code{n_obs}: Number of observations used in analysis
#'     \item \code{n_groups}: Number of groups (individuals) in the analysis
#'     \item \code{models}: If return_models=TRUE, includes the fitted lmer objects
#'     \item \code{random_effects}: Random intercept variances and residual variances
#'   }
#'
#' @details
#' Random Intercepts estimation works by:
#' \enumerate{
#'   \item Fitting mixed-effects models with random intercepts for each individual
#'   \item Allowing for correlation within individuals while estimating population-level effects
#'   \item Assuming individual intercepts are drawn from a normal distribution
#' }
#'
#' The approach handles unbalanced data well and provides more efficient estimates than
#' Fixed Individual Effects when the random intercept assumption is appropriate.
#'
#' Model specifications:
#' \itemize{
#'   \item Y equation: \code{y ~ xlag + ylag + (1|id)}
#'   \item X equation: \code{x ~ xlag + ylag + (1|id)}
#' }
#'
#' @examples
#' \dontrun{
#' # Generate some data
#' sim_data <- simRICLPM(waves = 4, sample.nobs = 500)$data
#' long_data <- reshape_long_sim_cr(sim_data)
#'
#' # Estimate Random Intercepts model
#' ri_results <- estimateRI(long_data)
#' print(ri_results)
#'
#' # Get full model objects
#' ri_full <- estimateRI(long_data, return_models = TRUE)
#' summary(ri_full$models$x_model)
#' summary(ri_full$models$y_model)
#'
#' # Check random effects
#' ri_full$random_effects
#' }
#'
#' @import lme4 stats
#' @export
estimateRI <- function(data,
                       y_outcome = "y",
                       x_outcome = "x",
                       y_lag = "ylag",
                       x_lag = "xlag",
                       id_var = "id",
                       return_models = FALSE,
                       control_params = lme4::lmerControl(optimizer = "bobyqa",
                                                          optCtrl = list(maxfun = 20000))) {

  # Load required libraries
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required but not installed. Please install it with: install.packages('lme4')")
  }

  # Input validation
  required_cols <- c(y_outcome, x_outcome, y_lag, x_lag, id_var)
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

  # Check for sufficient groups
  n_groups <- length(unique(analysis_data[[id_var]]))
  if (n_groups < 3) {
    warning("Very few groups (individuals) detected. Random intercepts model may not be appropriate.")
  }

  # Fit Random Intercepts models
  tryCatch({
    # Y equation: y ~ xlag + ylag + (1|id)
    formula_y <- as.formula(paste(y_outcome, "~", x_lag, "+", y_lag, "+ (1|", id_var, ")"))
    model_fit_y <- lme4::lmer(formula_y, data = analysis_data, control = control_params)

    # X equation: x ~ xlag + ylag + (1|id)
    formula_x <- as.formula(paste(x_outcome, "~", x_lag, "+", y_lag, "+ (1|", id_var, ")"))
    model_fit_x <- lme4::lmer(formula_x, data = analysis_data, control = control_params)

    # Check convergence
    converged_y <- !is.null(model_fit_y@optinfo$conv$lme4) &&
      model_fit_y@optinfo$conv$lme4$code == 0
    converged_x <- !is.null(model_fit_x@optinfo$conv$lme4) &&
      model_fit_x@optinfo$conv$lme4$code == 0
    converged <- converged_y && converged_x

    # Extract coefficients safely
    get_coef_safe <- function(model, coef_name) {
      fixed_effects <- lme4::fixef(model)
      if (coef_name %in% names(fixed_effects) && !is.na(fixed_effects[coef_name])) {
        return(as.numeric(fixed_effects[coef_name]))
      } else {
        return(NA_real_)
      }
    }

    # Extract coefficients
    xlag_x <- get_coef_safe(model_fit_x, x_lag)   # Autoregressive X effect
    ylag_x <- get_coef_safe(model_fit_x, y_lag)   # Cross-lagged Y -> X effect
    xlag_y <- get_coef_safe(model_fit_y, x_lag)   # Cross-lagged X -> Y effect
    ylag_y <- get_coef_safe(model_fit_y, y_lag)   # Autoregressive Y effect

    # Extract random effects information
    random_effects <- list(
      x_model = list(
        random_intercept_var = as.numeric(lme4::VarCorr(model_fit_x)[[id_var]][1, 1]),
        residual_var = attr(lme4::VarCorr(model_fit_x), "sc")^2
      ),
      y_model = list(
        random_intercept_var = as.numeric(lme4::VarCorr(model_fit_y)[[id_var]][1, 1]),
        residual_var = attr(lme4::VarCorr(model_fit_y), "sc")^2
      )
    )

    # Prepare results
    results <- list(
      xlag_x = xlag_x,
      ylag_x = ylag_x,
      xlag_y = xlag_y,
      ylag_y = ylag_y,
      converged = converged,
      n_obs = nrow(analysis_data),
      n_groups = n_groups,
      method = "Random Intercepts (lmer)",
      random_effects = random_effects
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
      n_groups = n_groups,
      method = "Random Intercepts (lmer)",
      error_message = as.character(e$message),
      random_effects = list(
        x_model = list(random_intercept_var = NA_real_, residual_var = NA_real_),
        y_model = list(random_intercept_var = NA_real_, residual_var = NA_real_)
      )
    )

    if (return_models) {
      results$models <- NULL
    }

    warning(paste("Random Intercepts estimation failed:", e$message))
    return(results)
  })
}
