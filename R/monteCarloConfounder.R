#' @title monteCarloConfounder
#' @description Monte Carlo simulation for Cross-Lagged Panel Model with unmeasured confounder.
#' Tests how different estimators perform when data is generated with an unmeasured confounder
#' that affects both X and Y variables across time.
#'
#' @param trials The number of trials for the Monte Carlo simulation.
#' @param waves The number of waves (time points) in the model.
#' @param estimator The estimation method to use. Options: "OLS", "RICLPM", "CLPM".
#' @param stability_p The stability parameter for the X variable (autoregressive effect).
#' @param stability_q The stability parameter for the Y variable (autoregressive effect).
#' @param cross_p The cross-lagged effect of Y on X at the next time point.
#' @param cross_q The cross-lagged effect of X on Y at the next time point.
#' @param variance_p The variance of the p latent variable.
#' @param variance_q The variance of the q latent variable.
#' @param cov_pq The covariance between X and Y within the same time point.
#' @param confounder_p The effect of the confounder on X variables.
#' @param confounder_q The effect of the confounder on Y variables.
#' @param confounder_variance The variance of the confounder.
#' @param confounder_stability The stability parameter for the confounder (autoregressive effect).
#' @param sample_size The sample size for each simulation.
#' @param include_confounder Logical. Whether to include the confounder in data generation.
#' @param ... Additional arguments.
#'
#' @return A data frame containing the results of the simulation with columns for:
#'   estimated coefficients, true parameters, trial number, estimator used, and bias measures.
#'
#' @details
#' This function generates data using simCLPMu (which includes an unmeasured confounder)
#' and then fits models that may or may not account for this confounder. This allows
#' testing of bias introduced by unmeasured confounders.
#'
#' The function tests the following scenario:
#' 1. Data is generated with an unmeasured confounder affecting both X and Y
#' 2. Models are fitted that ignore this confounder
#' 3. Bias in parameter estimates is measured
#'
#' @examples
#' \dontrun{
#' # Test bias from unmeasured confounder
#' results <- monteCarloConfounder(
#'   trials = 100,
#'   waves = 4,
#'   estimator = "RICLPM",
#'   confounder_p = 0.3,
#'   confounder_q = 0.3,
#'   sample_size = 1000
#' )
#'
#' # Compare multiple estimators
#' estimators <- c("OLS", "RICLPM", "CLPM")
#' results_list <- lapply(estimators, function(est) {
#'   monteCarloConfounder(trials = 50, estimator = est)
#' })
#' }
#'
#' @import tidyr dplyr lavaan
#' @export
monteCarloConfounder <- function(
    trials = 10,
    waves = 5,
    estimator = "RICLPM",
    stability_p = 0.2,
    stability_q = 0.2,
    cross_p = 0.1,
    cross_q = 0.1,
    variance_p = 1,
    variance_q = 1,
    cov_pq = 0.1,
    confounder_p = 0.3,
    confounder_q = 0.3,
    confounder_variance = 1,
    confounder_stability = 0.4,
    sample_size = 2500,
    include_confounder = TRUE,
    ...
) {

  # Validate inputs
  valid_estimators <- c("OLS", "RICLPM", "CLPM")
  if (!estimator %in% valid_estimators) {
    stop("estimator must be one of: ", paste(valid_estimators, collapse = ", "))
  }

  # Create parameter grid - for now just single combination
  # Could be expanded to test multiple parameter combinations
  simulation_parameters <- data.frame(
    stability_p = stability_p,
    stability_q = stability_q,
    cross_p = cross_p,
    cross_q = cross_q,
    variance_p = variance_p,
    variance_q = variance_q,
    cov_pq = cov_pq,
    confounder_p = confounder_p,
    confounder_q = confounder_q,
    confounder_variance = confounder_variance,
    confounder_stability = confounder_stability
  )

  # Initialize results list
  results <- list()

  # Loop through parameter combinations and trials
  for (i in 1:nrow(simulation_parameters)) {
    for (j in 1:trials) {
      params <- simulation_parameters[i, ]

      # Generate data with confounder using simCLPMu
      dat <- simCLPMu(
        waves = waves,
        stability_p = params$stability_p,
        stability_q = params$stability_q,
        cross_p = params$cross_p,
        cross_q = params$cross_q,
        variance_p = params$variance_p,
        variance_q = params$variance_q,
        cov_pq = params$cov_pq,
        include_confounder = include_confounder,
        confounder_p = params$confounder_p,
        confounder_q = params$confounder_q,
        confounder_variance = params$confounder_variance,
        confounder_stability = params$confounder_stability,
        sample.nobs = sample_size,
        ...
      )$data

      # Fit model using specified estimator (ignoring the confounder)
      tryCatch({
        if (estimator == "OLS") {
          # OLS estimation using reshaped data
          dat_long <- reshape_long_sim_cr(dat) %>%
            na.omit()

          model_fit_y <- lm(y ~ xlag + ylag, dat_long)
          model_fit_x <- lm(x ~ xlag + ylag, dat_long)

          cross_lag_y <- model_fit_y$coefficients[["xlag"]]
          auto_regressive_y <- model_fit_y$coefficients[["ylag"]]
          cross_lag_x <- model_fit_x$coefficients[["ylag"]]
          auto_regressive_x <- model_fit_x$coefficients[["xlag"]]

        } else if (estimator == "RICLPM") {
          # RICLPM estimation
          model <- lavaan::lavaan(
            estimateRICLPM(
              time_varying_x = paste0("x", 1:waves),
              time_varying_y = paste0("y", 1:waves),
              waves = waves
            ),
            data = dat
          )

          cross_lag_y <- coef(model)["cl_yeqn"]
          auto_regressive_y <- coef(model)["ar_yeqn"]
          cross_lag_x <- coef(model)["cl_xeqn"]
          auto_regressive_x <- coef(model)["ar_xeqn"]

        } else if (estimator == "CLPM") {
          # CLPM estimation
          model <- lavaan::lavaan(
            estimateCLPM(waves = waves),
            data = dat
          )

          cross_lag_y <- coef(model)["cl_yeqn"]
          auto_regressive_y <- coef(model)["ar_yeqn"]
          cross_lag_x <- coef(model)["cl_xeqn"]
          auto_regressive_x <- coef(model)["ar_xeqn"]
        }

        # Calculate bias (difference between estimated and true values)
        bias_cross_p <- as.numeric(cross_lag_x) - params$cross_p
        bias_cross_q <- as.numeric(cross_lag_y) - params$cross_q
        bias_stability_p <- as.numeric(auto_regressive_x) - params$stability_p
        bias_stability_q <- as.numeric(auto_regressive_y) - params$stability_q

        # Store results
        results[[length(results) + 1]] <- data.frame(
          # True parameter values
          true_stability_p = params$stability_p,
          true_stability_q = params$stability_q,
          true_cross_p = params$cross_p,
          true_cross_q = params$cross_q,
          true_confounder_p = params$confounder_p,
          true_confounder_q = params$confounder_q,
          true_confounder_variance = params$confounder_variance,
          true_confounder_stability = params$confounder_stability,

          # Estimated values
          est_stability_p = as.numeric(auto_regressive_x),
          est_stability_q = as.numeric(auto_regressive_y),
          est_cross_p = as.numeric(cross_lag_x),
          est_cross_q = as.numeric(cross_lag_y),

          # Bias measures
          bias_stability_p = bias_stability_p,
          bias_stability_q = bias_stability_q,
          bias_cross_p = bias_cross_p,
          bias_cross_q = bias_cross_q,

          # Meta information
          trial = j,
          estimator = estimator,
          dgp = "clpmu_confounder",
          waves = waves,
          sample_size = sample_size,
          include_confounder = include_confounder,

          stringsAsFactors = FALSE
        )

      }, error = function(e) {
        # Store error information
        results[[length(results) + 1]] <- data.frame(
          trial = j,
          estimator = estimator,
          dgp = "clpmu_confounder",
          error_occurred = TRUE,
          error_message = as.character(e$message),
          stringsAsFactors = FALSE
        )
      })
    }
  }

  # Combine results
  results_df <- do.call(rbind, results)

  return(results_df)
}
