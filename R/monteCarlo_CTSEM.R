#' @title monteCarloCTSEM
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
#' @import tidyr dplyr lavaan
#' @importFrom dplyr %>% select mutate row_number group_by ungroup

#'
#' @return A data frame containing the results of the simulation.
#'
#' @export
monteCarloCTSEM <- function(
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
    dgp = "riclpm",
    cores = 10,
    ...
) {
  if(dgp == "riclpm") {
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
          sample.nobs = sample_size,
          stability.p =  params$stability.p,
          stability.q =  params$stability.q,
          cross.p = params$cross.p,
          cross.q = params$cross.q,
          variance.p = params$variance.p,
          variance.q = params$variance.q,
          variance.between.x = params$variance.between.x,
          variance.between.y = params$variance.between.y,
        )$data %>%
          reshape_long_sim_cr() %>%
          select(id, wave, x, y) -> ct_sem_data

        names(ct_sem_data) = c("id", "time", "x", "y")
        model = estimateCTSEM(data = ct_sem_data,
)


        cross_lag_y =   model[[2]]$yx
        auto_regressive_y =  model[[2]]$yy

        cross_lag_x =   model[[2]]$xy
        auto_regressive_x = model[[2]]$xx

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
          xlag_x = auto_regressive_x,
          ylag_x = cross_lag_x,
          xlag_y = cross_lag_y,
          ylag_y = auto_regressive_y,
          trial = j,
          estimator = "ctsem",
          dgp = "riclpm"
        )
      }
    }
    results_df <- do.call(rbind, results)

  }
  else if(dgp == "clpm") {
    simulation_parameters <- expand.grid(
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

        dat <- simCLPM(
          waves = waves,
          sample.nobs = sample_size,
          stability.p =  params$stability.p,
          stability.q =  params$stability.q,
          cross.p = params$cross.p,
          cross.q = params$cross.q,
          variance.p = params$variance.p,
          variance.q = params$variance.q,
        )$data  %>%
          reshape_long_sim_cr() %>%
          select(id, wave, x, y) -> ct_sem_data

        names(ct_sem_data) = c("id", "time", "x", "y")
        model = estimateCTSEM(data = ct_sem_data)


        cross_lag_y =   model[[2]]$yx
        auto_regressive_y =  model[[2]]$yy

        cross_lag_x =   model[[2]]$xy
        auto_regressive_x = model[[2]]$xx

        # Store results in the list
        results[[length(results) + 1]] <- data.frame(
          stability.p =  params$stability.p,
          stability.q =  params$stability.q,
          cross.p = params$cross.p,
          cross.q = params$cross.q,
          variance.p = params$variance.p,
          variance.q = params$variance.q,
          xlag_x = auto_regressive_x,
          ylag_x = cross_lag_x,
          xlag_y = cross_lag_y,
          ylag_y = auto_regressive_y,
          trial = j,
          estimator = "ctsem",
          dgp = "clpm"

        )
      }
    }
    results_df <- do.call(rbind, results)
  }

  return(results_df)

}



