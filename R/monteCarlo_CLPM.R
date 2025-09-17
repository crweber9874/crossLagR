#' @title monteCarloCLPM
#' @description Monte Carlo simulation for Cross-Lagged Panel Model
#'
#' @param trials The number of trials for the Monte Carlo simulation.
#' @param waves The number of waves (time points) in the model.
#' @param stability_q Autoregressive effect for q (y) variable
#' @param stability_p Autoregressive effect for p (x) variable
#' @param variance_between_y Between-person variance for y variable
#' @param variance_between_x Between-person variance for x variable
#' @param cross_p Cross-lagged effect of y on x
#' @param cross_q Cross-lagged effect of x on y
#' @param variance_p Within-person variance for p variable
#' @param variance_q Within-person variance for q variable
#' @param sample_size Sample size for simulation
#' @param dgp Data generation process: "riclpm", "clpm", or "clpmu"
#' @param verbose Whether to print progress messages
#' @param confounder_type Character string specifying confounder type for clpmu dgp.
#'   Must be one of "time_variant" (default) or "time_invariant".
#' @param confounder_p The effect of the confounder on the X variable.
#' @param confounder_q The effect of the confounder on the Y variable.
#' @param confounder_variance The variance of the confounder.
#' @param confounder_stability The stability of the time-variant confounder.
#'
#' @import tidyr dplyr lavaan
#' @return A data frame containing the results of the simulation.
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
    verbose = FALSE,
    confounder_type = "time_variant",
    confounder_p = 0.3,
    confounder_q = 0.3,
    confounder_variance = 1,
    confounder_stability = 0.4,
    cov_pq = 0.1,
    ...
) {

  # Validate confounder_type
  if (!confounder_type %in% c("time_variant", "time_invariant")) {
    stop("confounder_type must be either 'time_variant' or 'time_invariant'")
  }

  # Initialize an empty list to store results
  results <- list()

  if (dgp == "riclpm") {
    simulation_parameters <- expand.grid(
      variance_between_x = variance_between_x,
      variance_between_y = variance_between_y,
      stability_p = stability_p,
      stability_q = stability_q,
      cross_p = cross_p,
      cross_q = cross_q,
      variance_p = variance_p,
      variance_q = variance_q
    ) %>% as.data.frame()

    # Loop through each combination of parameters and trials
    for (i in 1:nrow(simulation_parameters)) {
      for (j in 1:trials) {
        params <- simulation_parameters[i, ]
        if (verbose) {
          message(paste0("DGP: ", dgp, " | Simulation ", i, " of ", nrow(simulation_parameters), " | Trial ", j, " of ", trials))
        }
        tryCatch({
          # Generate data
          dat <- simRICLPM(
            waves = waves,
            stability_p = params$stability_p,
            stability_q = params$stability_q,
            cross_p = params$cross_p,
            cross_q = params$cross_q,
            variance_p = params$variance_p,
            variance_q = params$variance_q,
            variance_between_x = params$variance_between_x,
            variance_between_y = params$variance_between_y,
            sample.nobs = sample_size
          )$data

          # Estimate model
          model <- lavaan::lavaan(
            estimateCLPM(waves = waves),
            data = dat
          )

          # Extract coefficients from the parameter table - CORRECTED MAPPING
          param_table <- lavaan::parameterEstimates(model)

          # IMPORTANT: The estimateCLPM function uses confusing parameter names:
          # - "ar_yeqn" is actually the X autoregressive parameter
          # - "ar_xeqn" is actually the Y autoregressive parameter
          # - "cl_yeqn" is actually the Y→X cross-lagged parameter
          # - "cl_xeqn" is actually the X→Y cross-lagged parameter

          ar_x <- param_table[param_table$label == "ar_yeqn", "est"]      # X autoregressive
          ar_y <- param_table[param_table$label == "ar_xeqn", "est"]      # Y autoregressive
          cl_x_to_y <- param_table[param_table$label == "cl_xeqn", "est"] # X→Y cross-lagged
          cl_y_to_x <- param_table[param_table$label == "cl_yeqn", "est"] # Y→X cross-lagged

          # Store results in the list with CORRECT mapping
          results[[length(results) + 1]] <- data.frame(
            variance_between_x = params$variance_between_x,
            variance_between_y = params$variance_between_y,
            stability_p = params$stability_p,
            stability_q = params$stability_q,
            cross_p = params$cross_p,
            cross_q = params$cross_q,
            variance_p = params$variance_p,
            variance_q = params$variance_q,
            xlag_x = as.numeric(ar_x),        # X autoregressive (corrected)
            ylag_x = as.numeric(cl_y_to_x),   # Y→X cross-lagged (corrected)
            xlag_y = as.numeric(cl_x_to_y),   # X→Y cross-lagged (corrected)
            ylag_y = as.numeric(ar_y),        # Y autoregressive (corrected)
            trial = j,
            estimator = "clpm",
            dgp = "riclpm"
          )

        }, error = function(e) {
          # Store NA results for failed trials
          results[[length(results) + 1]] <<- data.frame(
            variance_between_x = params$variance_between_x,
            variance_between_y = params$variance_between_y,
            stability_p = params$stability_p,
            stability_q = params$stability_q,
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

  } else if (dgp == "clpm") {
    simulation_parameters <- expand.grid(
      stability_p = stability_p,
      stability_q = stability_q,
      cross_p = cross_p,
      cross_q = cross_q,
      variance_p = variance_p,
      variance_q = variance_q
    ) %>% as.data.frame()

    for (i in 1:nrow(simulation_parameters)) {
      for (j in 1:trials) {
        params <- simulation_parameters[i, ]

        tryCatch({
          dat <- simCLPM(
            waves = waves,
            stability_p = params$stability_p,
            stability_q = params$stability_q,
            cross_p = params$cross_p,
            cross_q = params$cross_q,
            variance_p = params$variance_p,
            variance_q = params$variance_q,
            cov_pq = cov_pq,
            sample.nobs = sample_size
          )$data

          model <- lavaan::lavaan(
            estimateCLPM(waves = waves),
            data = dat
          )

          # Extract coefficients with CORRECTED mapping
          param_table <- lavaan::parameterEstimates(model)

          ar_x <- param_table[param_table$label == "ar_yeqn", "est"]      # X autoregressive
          ar_y <- param_table[param_table$label == "ar_xeqn", "est"]      # Y autoregressive
          cl_x_to_y <- param_table[param_table$label == "cl_xeqn", "est"] # X→Y cross-lagged
          cl_y_to_x <- param_table[param_table$label == "cl_yeqn", "est"] # Y→X cross-lagged

          results[[length(results) + 1]] <- data.frame(
            stability_p = params$stability_p,
            stability_q = params$stability_q,
            cross_p = params$cross_p,
            cross_q = params$cross_q,
            variance_p = params$variance_p,
            variance_q = params$variance_q,
            xlag_x = as.numeric(ar_x),        # X autoregressive (corrected)
            ylag_x = as.numeric(cl_y_to_x),   # Y→X cross-lagged (corrected)
            xlag_y = as.numeric(cl_x_to_y),   # X→Y cross-lagged (corrected)
            ylag_y = as.numeric(ar_y),        # Y autoregressive (corrected)
            trial = j,
            estimator = "clpm",
            dgp = "clpm"
          )
        }, error = function(e) {
          results[[length(results) + 1]] <<- data.frame(
            stability_p = params$stability_p,
            stability_q = params$stability_q,
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

  } else if (dgp == "clpmu") {
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
        confounder_stability = confounder_stability
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
        confounder_variance = confounder_variance
      ) %>% as.data.frame()
    }

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
              cov_pq = cov_pq,
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
              cov_pq = cov_pq,
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

          # Extract coefficients with CORRECTED mapping
          param_table <- lavaan::parameterEstimates(model)

          ar_x <- param_table[param_table$label == "ar_yeqn", "est"]      # X autoregressive
          ar_y <- param_table[param_table$label == "ar_xeqn", "est"]      # Y autoregressive
          cl_x_to_y <- param_table[param_table$label == "cl_xeqn", "est"] # X→Y cross-lagged
          cl_y_to_x <- param_table[param_table$label == "cl_yeqn", "est"] # Y→X cross-lagged

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
            include_confounder = TRUE, # Confounding is always included here
            confounder_type = confounder_type,
            xlag_x = as.numeric(ar_x),        # X autoregressive (corrected)
            ylag_x = as.numeric(cl_y_to_x),   # Y→X cross-lagged (corrected)
            xlag_y = as.numeric(cl_x_to_y),   # X→Y cross-lagged (corrected)
            ylag_y = as.numeric(ar_y),        # Y autoregressive (corrected)
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
            include_confounder = TRUE, # Confounding is always included here
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
