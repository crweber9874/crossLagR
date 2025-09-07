#' @title monteCarloRICLPM
#' @description A somewhat strange function. It's common to estimate the cross-lagged panel model with two OLS regression models.
#' allows the user to specify different Data Generating conditions, and apply the OLS model to the data.
#' The output is a data frame that includes the estimated coefficients across these trials
#' @param trials The number of trials for the Monte Carlo simulation.
#' @param waves The number of waves (time points) in the model.
#' @param model Specify whether an "ols" or "clpm" model.
#' @param data_generation Specify the data generation method: "clpm", "ri-clpm"
#' @param proportion_change_x Proportion of change in x.
#' @param proportion_change_y Proportion of change in y.
#'
#' @import tidyr dplyr lavaan
#' @importFrom dplyr %>% select mutate row_number group_by ungroup

#'
#' @return A data frame containing the results of the simulation.
#'
#' @export
monteCarloLChange<- function(
    trials = 10,
    waves = 5,
    variance_between_x = 0.5,
    variance_between_y = 0.5,
    stability_q = 0.25,
    stability_p = 0.25,
    cross_p = 0,
    cross_q = 0,
    variance_p = 1,
    variance_q = 1,
    sample_size = 2500,
    dgp = "riclpm",
    ...
) {
  if(dgp == "riclpm") {
    simulation_parameters <- expand.grid(
      variance_between_x = variance_between_x,
      variance_between_y = variance_between_y,
      stability_p = stability_p,
      stability_q = stability_q,
      cross_p = cross_p,
      cross_q = cross_q,
      variance_p = variance_p,
      variance_q = variance_q) %>% as.data.frame()

    # Initialize an empty list to store results
    results <- list()

    # Loop through each combination of parameters and trials
    for (i in 1:nrow(simulation_parameters)) {
      for (j in 1:trials) {
        params <- simulation_parameters[i, ]

        dat <- simRICLPM(
          waves = waves,
          stability_p =  params$stability_p,
          stability_q =  params$stability_q,
          cross_p = params$cross_p,
          cross_q = params$cross_q,
          variance_p = params$variance_p,
          variance_q = params$variance_q,
          variance_between_x = params$variance_between_x,
          variance_between_y = params$variance_between_y,
          sample.nobs = sample_size
        )$data -> df

        suppressWarnings({

          lavaan::lavaan(
            latentChange(
              waves = waves,
              model_type = "dual_change"),
            meanstructure = FALSE,
            estimator = "ML",
            missing = "fiml",
            fixed.x = FALSE,
            mimic="mplus",
            control=list(iter.max=1000),
            verbose=FALSE,
            data = df) -> model
        })

        cross_lag_y = coef(model)["change.x"]
        auto_regressive_y = coef(model)["ar_y"]

        cross_lag_x = coef(model)["change.y"]
        auto_regressive_x = coef(model)["ar_x"]

        # Store results in the list
        results[[length(results) + 1]] <- data.frame(
          variance_between_x = params$variance_between_x,
          variance_between_y = params$variance_between_y,
          stability_p =  params$stability_p,
          stability_q =  params$stability_q,
          cross_p = params$cross_p,
          cross_q = params$cross_q,
          variance_p = params$variance_p,
          variance_q = params$variance_q,
          xlag_x = auto_regressive_x,
          ylag_x = cross_lag_x,
          xlag_y = cross_lag_y,
          ylag_y = auto_regressive_y,
          trial = j,
          estimator = "lchange",
          dgp = "riclpm"

        )
      }
    }
    results_df <- do.call(rbind, results)

  }
  else if(dgp == "clpm") {
    simulation_parameters <- expand.grid(
      stability_p = stability_p,
      stability_q = stability_q,
      cross_p = cross_p,
      cross_q = cross_q,
      variance_p = variance_p,
      variance_q = variance_q) %>% as.data.frame()

    # Initialize an empty list to store results
    results <- list()

    # Loop through each combination of parameters and trials
    for (i in 1:nrow(simulation_parameters)) {
      for (j in 1:trials) {
        params <- simulation_parameters[i, ]

        dat <- simCLPM(
          waves = waves,
          stability_p =  params$stability_p,
          stability_q =  params$stability_q,
          cross_p = params$cross_p,
          cross_q = params$cross_q,
          variance_p = params$variance_p,
          variance_q = params$variance_q,
          sample.nobs = sample_size
        )$data -> df

        lavaan::lavaan(
          estimateRICLPM(
            time_varying_x = paste0("x", c(1:waves)),
            time_varying_y = paste0("y", c(1:waves)),
            waves = waves),
          data = df) -> model

        cross_lag_y = lavaan::coef(model)["cl_xeqn"]
        auto_regressive_y = lavaan::coef(model)["ar_yeqn"]

        cross_lag_x = lavaan::coef(model)["cl_yeqn"]
        auto_regressive_x = lavaan::coef(model)["ar_xeqn"]

        # Store results in the list
        results[[length(results) + 1]] <- data.frame(
          stability_p =  params$stability_p,
          stability_q =  params$stability_q,
          cross_p = params$cross_p,
          cross_q = params$cross_q,
          variance_p = params$variance_p,
          variance_q = params$variance_q,
          xlag_x = auto_regressive_x,
          ylag_x = cross_lag_x,
          xlag_y = cross_lag_y,
          ylag_y = auto_regressive_y,
          trial = j,
          estimator = "riclpm",
          dgp = "clpm"

        )
      }
    }
    results_df <- do.call(rbind, results)

  }

  return(results_df)

}

