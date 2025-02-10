#' @title monteCarlo_OLS
#' @description A somewhat strange function. It's common to estimate the cross-lagged panel model with two OLS regression models.
#' allows the user to specify different Data Generating conditions, and apply the OLS model to the data.
#' The output is a data frame that includes the estimated coefficients across these trials
#' @param trials The number of trials for the Monte Carlo simulation.
#' @param waves The number of waves (time points) in the model.
#' @param model Specify whether an "ols" or "clpm" model.
#' @param data_generation Specify the data generation method: "clpm", "ri-clpm"
#' @param proportion.change_x Proportion of change in x.
#' @param proportion.change_y Proportion of change in y.
#'
#' @import tidyr dplyr

#'
#' @return A data frame containing the results of the simulation.
#'
#' @export
monteCarlo_OLS <- function(
    trials = 10,
    waves = 5,
    data_generation = c("ri-clpm", "ri-clpms"),
    variance.between.x = 0.5,
    variance.between.y = 0.5,
    within_person_stability_x = 0,
    within_person_stability_y = 0,
    cross_lag_x = 0,
    cross_lag_y = 0,
    variance_x = 1,
    variance_y = 1,
    sample_size = 2500,
    cov.between = 0.5,
    confound = 0.5,
    ...
) {
  data_generation <- match.arg(data_generation)

    if (data_generation == "ri-clpm") {
      # Create a data frame to store simulation parameters
      simulation_parameters <- expand.grid(
            variance.between.x = variance.between.x,
            variance.between.y = variance.between.y,
            cov.between = cov.between) %>% as.data.frame()

      # Initialize an empty list to store results
      results <- list()

      # Loop through each combination of parameters and trials
      for (i in 1:nrow(simulation_parameters)) {
        for (j in 1:trials) {
          params <- simulation_parameters[i, ]

          dat <- simulate_riclpm(
            waves = waves,
            stability.p = 0,
            stability.q = 0,
            cross.p = 0,
            cross.q = 0,
            sample.nobs = sample_size,
            variance.p = 1,
            variance.q = 1,
            variance.between.x = params$variance.between.x,
            variance.between.y = params$variance.between.y,
            cov.between = params$cov.between
          )$data %>%
            reshape_long_sim_cr() %>% as.data.frame() %>% na.omit()

          # Fit the OLS model and extract coefficients
          model_fit_y <- lm(y ~ xlag + ylag, dat)
          model_fit_x <- lm(x ~ xlag + ylag, dat)

          xvalue_x <-  model_fit_x$coefficients[["xlag"]]
          yvalue_x <-  model_fit_x$coefficients[["ylag"]]

          xvalue_y <-  model_fit_y$coefficients[["xlag"]]
          yvalue_y <-  model_fit_y$coefficients[["ylag"]]

          # model_fit <- lm(y ~ -1 + xlag + ylag, dat)
          # xvalue <- coef(model_fit)[["xlag"]]
          # yvalue <- coef(model_fit)[["ylag"]]

          # Store results in the list
          results[[length(results) + 1]] <- data.frame(
            variance.between.x = params$variance.between.x,
            variance.between.y = params$variance.between.y,
            cov.between = params$cov.between,
            xlag_x = xvalue_x,
            ylag_x = yvalue_x,
            xlag_y = xvalue_y,
            ylag_y = yvalue_y,
            trial = j,
            model = "ri-clpm",
            data_generation = data_generation
          )
        }
      }
      results_df <- do.call(rbind, results)
    }
    if (data_generation == "ri-clpms") {
    # A confounder affecting the intercepts
    simulation_parameters <- expand.grid(
      variance.between.x = variance.between.x,
      variance.between.y = variance.between.y,
      cov.between = cov.between,
      confound = confound) %>% as.data.frame()

    # Initialize an empty list to store results
    results <- list()

    # Loop through each combination of parameters and trials
    for (i in 1:nrow(simulation_parameters)) {
      for (j in 1:trials) {
        params <- simulation_parameters[i, ]

        dat <- simulate_riclpm_confound(
          waves = waves,
          confound = params$confound,
          stability.p = 0,
          stability.q = 0,
          cross.p = 0,
          cross.q = 0,
          sample.nobs = sample_size,
          variance.p = 1,
          variance.q = 1,
          variance.between.x = params$variance.between.x,
          variance.between.y = params$variance.between.y,
          cov.between = params$cov.between
        )$data %>%
          reshape_long_sim_cr() %>% as.data.frame() %>% na.omit()

        # Fit the OLS model and extract coefficients
        model_fit_y <- lm(y ~ xlag + ylag, dat)
        model_fit_x <- lm(x ~ xlag + ylag, dat)

        xvalue_x <-  model_fit_x$coefficients[["xlag"]]
        yvalue_x <-  model_fit_x$coefficients[["ylag"]]

        xvalue_y <-  model_fit_y$coefficients[["xlag"]]
        yvalue_y <-  model_fit_y$coefficients[["ylag"]]

        results[[length(results) + 1]] <- data.frame(
          variance.between.x = params$variance.between.x,
          variance.between.y = params$variance.between.y,
          cov.between = params$cov.between,
          confound = params$confound,
          xlag_x = xvalue_x,
          ylag_x = yvalue_x,
          xlag_y = xvalue_y,
          ylag_y = yvalue_y,
          trial = j,
          model = "ri-clpms",
          data_generation = data_generation
        )
      }
    }
    results_df <- do.call(rbind, results)
  }
    return(results_df)

    }



