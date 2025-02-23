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
    variance.between.x = 0.5,
    variance.between.y = 0.5,
    stability.q = 0.25,
    stability.p = 0.25,
    cross.p = 0,
    cross.q = 0,
    variance.p = 1,
    variance.q = 1,
    sample_size = 2500,
    ...
) {

    simulation_parameters <- expand.grid(
            variance.between.x = variance.between.x,
            variance.between.y = variance.between.y,
            stability.p = stability.p,
            stability.q = stability.q,
            cross.p = cross.p,
            cross.q = cross.q,
            variance.p = variance.p,
            variance.q = variance.q) %>% as.data.frame()

      # Initialize an empty list to store results
      results <- list()

      # Loop through each combination of parameters and trials
      for (i in 1:nrow(simulation_parameters)) {
        for (j in 1:trials) {
          params <- simulation_parameters[i, ]

          dat <- simRICLPM(
            waves = waves,
            sample.nobs = 2500,
            stability.p =  params$stability.p,
            stability.q =  params$stability.q,
            cross.p = params$cross.p,
            cross.q = params$cross.q,
            variance.p = params$variance.p,
            variance.q = params$variance.q,
            variance.between.x = params$variance.between.x,
            variance.between.y = params$variance.between.y,
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
            stability.p =  params$stability.p,
            stability.q =  params$stability.q,
            cross.p = params$cross.p,
            cross.q = params$cross.q,
            variance.p = params$variance.p,
            variance.q = params$variance.q,
            xlag_x = xvalue_x,
            ylag_x = yvalue_x,
            xlag_y = xvalue_y,
            ylag_y = yvalue_y,
            trial = j
          )
        }
      }
      results_df <- do.call(rbind, results)

    return(results_df)

    }



