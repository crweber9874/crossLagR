#' @title monteCarloCLPM
#' @description A somewhat strange function. It's common to estimate the cross-lagged panel model with two OLS regression models.
#' allows the user to specify different Data Generating conditions, and apply the OLS model to the data.
#' The output is a data frame that includes the estimated coefficients across these trials
#' @param trials The number of trials for the Monte Carlo simulation.
#' @param waves The number of waves (time points) in the model.
#' @param model Specify whether an "ols" or "clpm" model.
#' @param data_generation Specify the data generation method: "clpm", "ri-clpm"
#' @param proportion_change_x Proportion of change in x.
#' @param proportion_change_y Proportion of change in y.
#' @param confounder_type Character string specifying confounder type for clpmu dgp.
#'   Must be one of "time_variant" (default) or "time_invariant".
#'
#' @import tidyr
#' @import dplyr
#' @import lavaan
#' @importFrom dplyr %>% select mutate row_number group_by ungroup

#'
#' @return A data frame containing the results of the simulation.
#'
#' @export
monteCarloCLPM <- function(
    trials = 10,
    waves = 5,
    stability_q = 0.25,
    stability_p = 0.25,
    variance_between_y = 0.5,
    variance_between_x = 0.5,
    cross_p = 0,
    cross_q = 0,
    variance_p = 1,
    variance_q = 1,
    sample_size = 2500,
    dgp = "riclpm",
    # Add confounder parameters
    confounder_p = 0.3,
    confounder_q = 0.3,
    confounder_variance = 1,
    confounder_stability = 0.4,
    include_confounder = TRUE,
    confounder_type = "time_variant",  # NEW PARAMETER
    ...
) {
  # Validate confounder_type
  if (!confounder_type %in% c("time_variant", "time_invariant")) {
    stop("confounder_type must be either 'time_variant' or 'time_invariant'")
  }

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

        tryCatch({
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
          )$data

          model <- lavaan::lavaan(
            estimateCLPM(waves = waves),
            data = dat
          )

          cross_lag_y <- coef(model)["cl_yeqn"]
          auto_regressive_y <- coef(model)["ar_yeqn"]
          cross_lag_x <- coef(model)["cl_xeqn"]
          auto_regressive_x <- coef(model)["ar_xeqn"]

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
            xlag_x = as.numeric(auto_regressive_x),
            ylag_x = as.numeric(cross_lag_x),
            xlag_y = as.numeric(cross_lag_y),
            ylag_y = as.numeric(auto_regressive_y),
            trial = j,
            estimator = "clpm",
            dgp = "riclpm"
          )
        }, error = function(e) {
          # Store NA results for failed trials
          results[[length(results) + 1]] <<- data.frame(
            variance_between_x = params$variance_between_x,
            variance_between_y = params$variance_between_y,
            stability_p =  params$stability_p,
            stability_q =  params$stability_q,
            cross_p = params$cross_p,
            cross_q = params$cross_q,
            variance_p = params$variance_p,
            variance_q = params$variance_q,
            xlag_x = NA,
            ylag_x = NA,
            xlag_y = NA,
            ylag_y = NA,
            trial = j,
            estimator = "clpm",
            dgp = "riclpm"
          )
        })
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

        tryCatch({
          dat <- simCLPM(
            waves = waves,
            stability_p =  params$stability_p,
            stability_q =  params$stability_q,
            cross_p = params$cross_p,
            cross_q = params$cross_q,
            variance_p = params$variance_p,
            variance_q = params$variance_q,
            sample.nobs = sample_size
          )$data

          model <- lavaan::lavaan(
            estimateCLPM(waves = waves),
            data = dat
          )

          cross_lag_y <- coef(model)["cl_yeqn"]
          auto_regressive_y <- coef(model)["ar_yeqn"]
          cross_lag_x <- coef(model)["cl_xeqn"]
          auto_regressive_x <- coef(model)["ar_xeqn"]

          # Store results in the list
          results[[length(results) + 1]] <- data.frame(
            stability_p =  params$stability_p,
            stability_q =  params$stability_q,
            cross_p = params$cross_p,
            cross_q = params$cross_q,
            variance_p = params$variance_p,
            variance_q = params$variance_q,
            xlag_x = as.numeric(auto_regressive_x),
            ylag_x = as.numeric(cross_lag_x),
            xlag_y = as.numeric(cross_lag_y),
            ylag_y = as.numeric(auto_regressive_y),
            trial = j,
            estimator = "clpm",
            dgp = "clpm"
          )
        }, error = function(e) {
          # Store NA results for failed trials
          results[[length(results) + 1]] <<- data.frame(
            stability_p =  params$stability_p,
            stability_q =  params$stability_q,
            cross_p = params$cross_p,
            cross_q = params$cross_q,
            variance_p = params$variance_p,
            variance_q = params$variance_q,
            xlag_x = NA,
            ylag_x = NA,
            xlag_y = NA,
            ylag_y = NA,
            trial = j,
            estimator = "clpm",
            dgp = "clpm"
          )
        })
      }
    }
    results_df <- do.call(rbind, results)

  }
  else if(dgp == "clpmu") {  # UPDATED SECTION FOR CONFOUNDER MODEL
    # Build simulation parameters based on confounder type
    if (confounder_type == "time_variant") {
      simulation_parameters <- expand.grid(
        stability_p = stability_p,
        stability_q = stability_q,
        cross_p = cross_p,
        cross_q = cross_q,
        variance_p = variance_p,
        variance_q = variance_q,
        confounder_p = confounder_p,
        confounder_q = confounder_q,
        confounder_variance = confounder_variance,
        confounder_stability = confounder_stability,
        include_confounder = include_confounder
      ) %>% as.data.frame()
    } else {  # time_invariant
      simulation_parameters <- expand.grid(
        stability_p = stability_p,
        stability_q = stability_q,
        cross_p = cross_p,
        cross_q = cross_q,
        variance_p = variance_p,
        variance_q = variance_q,
        confounder_p = confounder_p,
        confounder_q = confounder_q,
        confounder_variance = confounder_variance,
        include_confounder = include_confounder
        # Note: No confounder_stability for time-invariant
      ) %>% as.data.frame()
    }

    # Initialize an empty list to store results
    results <- list()

    # Loop through each combination of parameters and trials
    for (i in 1:nrow(simulation_parameters)) {
      for (j in 1:trials) {
        params <- simulation_parameters[i, ]

        tryCatch({
          # Choose data generation function based on confounder type
          if (confounder_type == "time_variant") {
            dat <- simCLPMu(
              waves = waves,
              stability_p = params$stability_p,
              stability_q = params$stability_q,
              cross_p = params$cross_p,
              cross_q = params$cross_q,
              variance_p = params$variance_p,
              variance_q = params$variance_q,
              cov_pq = 0.1,
              include_confounder = params$include_confounder,
              confounder_p = params$confounder_p,
              confounder_q = params$confounder_q,
              confounder_variance = params$confounder_variance,
              confounder_stability = params$confounder_stability,
              sample.nobs = sample_size
            )$data
          } else {  # time_invariant
            dat <- simCLPM_timeInvariantU(
              waves = waves,
              stability_p = params$stability_p,
              stability_q = params$stability_q,
              cross_p = params$cross_p,
              cross_q = params$cross_q,
              variance_p = params$variance_p,
              variance_q = params$variance_q,
              cov_pq = 0.1,
              include_confounder = params$include_confounder,
              confounder_p = params$confounder_p,
              confounder_q = params$confounder_q,
              confounder_variance = params$confounder_variance,
              sample.nobs = sample_size
            )$data
          }

          model <- lavaan::lavaan(
            estimateCLPM(waves = waves),
            data = dat
          )

          cross_lag_y <- coef(model)["cl_yeqn"]
          auto_regressive_y <- coef(model)["ar_yeqn"]
          cross_lag_x <- coef(model)["cl_xeqn"]
          auto_regressive_x <- coef(model)["ar_xeqn"]

          # Store results in the list - include confounder type info
          result_row <- data.frame(
            stability_p = params$stability_p,
            stability_q = params$stability_q,
            cross_p = params$cross_p,
            cross_q = params$cross_q,
            variance_p = params$variance_p,
            variance_q = params$variance_q,
            confounder_p = params$confounder_p,
            confounder_q = params$confounder_q,
            confounder_variance = params$confounder_variance,
            include_confounder = params$include_confounder,
            confounder_type = confounder_type,
            xlag_x = as.numeric(auto_regressive_x),
            ylag_x = as.numeric(cross_lag_x),
            xlag_y = as.numeric(cross_lag_y),
            ylag_y = as.numeric(auto_regressive_y),
            trial = j,
            estimator = "clpm",
            dgp = "clpmu"
          )

          # Add confounder_stability only for time_variant confounders
          if (confounder_type == "time_variant") {
            result_row$confounder_stability <- params$confounder_stability
          } else {
            result_row$confounder_stability <- NA
          }

          results[[length(results) + 1]] <- result_row

        }, error = function(e) {
          # Store NA results for failed trials
          result_row <- data.frame(
            stability_p = params$stability_p,
            stability_q = params$stability_q,
            cross_p = params$cross_p,
            cross_q = params$cross_q,
            variance_p = params$variance_p,
            variance_q = params$variance_q,
            confounder_p = params$confounder_p,
            confounder_q = params$confounder_q,
            confounder_variance = params$confounder_variance,
            include_confounder = params$include_confounder,
            confounder_type = confounder_type,
            xlag_x = NA,
            ylag_x = NA,
            xlag_y = NA,
            ylag_y = NA,
            trial = j,
            estimator = "clpm",
            dgp = "clpmu"
          )

          # Add confounder_stability
          if (confounder_type == "time_variant") {
            result_row$confounder_stability <- params$confounder_stability
          } else {
            result_row$confounder_stability <- NA
          }

          results[[length(results) + 1]] <<- result_row
        })
      }
    }
    results_df <- do.call(rbind, results)
  }

  return(results_df)
}
