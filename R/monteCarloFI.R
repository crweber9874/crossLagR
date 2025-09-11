#' @title monteCarloFixed
#' @description Monte Carlo simulation for Fixed Individual Effects estimation
#' Allows the user to specify different Data Generating conditions, and apply the FI model to the data.
#' The output is a data frame that includes the estimated coefficients across these trials
#' @param trials The number of trials for the Monte Carlo simulation.
#' @param waves The number of waves (time points) in the model.
#' @param dgp Specify the data generation method: "riclpm", "clpm", "clpmu"
#' @param confounder_type Character string specifying confounder type for clpmu dgp.
#'   Must be one of "time_variant" (default) or "time_invariant".
#' @param stability_p Autoregressive effect for p (x) variable
#' @param stability_q Autoregressive effect for q (y) variable
#' @param cross_p Cross-lagged effect of y on x
#' @param cross_q Cross-lagged effect of x on y
#' @param variance_p Within-person variance for p variable
#' @param variance_q Within-person variance for q variable
#' @param variance_between_x Between-person variance for x variable (for RICLPM data)
#' @param variance_between_y Between-person variance for y variable (for RICLPM data)
#' @param sample_size Sample size for simulation
#' @param confounder_p Effect of confounder on x variables (for clpmu)
#' @param confounder_q Effect of confounder on y variables (for clpmu)
#' @param confounder_variance Variance of the confounder (for clpmu)
#' @param confounder_stability Autoregressive effect of confounder (for clpmu)
#' @param include_confounder Whether to include confounder in clpmu model
#' @param ... Additional arguments
#'
#' @return A data frame containing the results of the simulation.
#'
#' @import dplyr tidyr stats
#' @export
monteCarloFixed <- function(
    trials = 10,
    waves = 5,
    variance_between_x = 0.5,
    variance_between_y = 0.5,
    stability_q = 0.75,
    stability_p = 0.75,
    cross_p = 0.5,
    cross_q = 0.5,
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
    confounder_type = "time_variant",
    ...
) {
  # Load required libraries explicitly
  library(dplyr)

  # Validate inputs
  valid_dgps <- c("riclpm", "clpm", "clpmu")
  if (!dgp %in% valid_dgps) {
    stop("dgp must be one of: ", paste(valid_dgps, collapse = ", "))
  }

  if (!confounder_type %in% c("time_variant", "time_invariant")) {
    stop("confounder_type must be either 'time_variant' or 'time_invariant'")
  }

  # Build simulation parameters based on dgp and confounder type
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
    )
    simulation_parameters <- as.data.frame(simulation_parameters)
  } else if (dgp == "clpm") {
    simulation_parameters <- expand.grid(
      stability_p = stability_p,
      stability_q = stability_q,
      cross_p = cross_p,
      cross_q = cross_q,
      variance_p = variance_p,
      variance_q = variance_q
    )
    simulation_parameters <- as.data.frame(simulation_parameters)
  } else if (dgp == "clpmu") {
    # Build parameters based on confounder type
    if (confounder_type == "time_variant") {
      simulation_parameters <- expand.grid(
        variance_between_x = variance_between_x,
        variance_between_y = variance_between_y,
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
      )
      simulation_parameters <- as.data.frame(simulation_parameters)
    } else { # time_invariant
      simulation_parameters <- expand.grid(
        variance_between_x = variance_between_x,
        variance_between_y = variance_between_y,
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
      )
      simulation_parameters <- as.data.frame(simulation_parameters)
    }
  }

  # Initialize results list
  results <- list()

  # Loop through parameter combinations and trials
  for (i in 1:nrow(simulation_parameters)) {
    for (j in 1:trials) {
      params <- simulation_parameters[i, ]

      # Step 1: Generate data based on dgp
      if (dgp == "riclpm") {
        raw_data <- simRICLPM(
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
      } else if (dgp == "clpm") {
        raw_data <- simCLPM(
          waves = waves,
          stability_p = params$stability_p,
          stability_q = params$stability_q,
          cross_p = params$cross_p,
          cross_q = params$cross_q,
          variance_p = params$variance_p,
          variance_q = params$variance_q,
          sample.nobs = sample_size
        )$data
      } else if (dgp == "clpmu") {
        # Choose data generation function based on confounder type
        if (confounder_type == "time_variant") {
          raw_data <- simCLPMu(
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
        } else { # time_invariant
          raw_data <- simCLPM_timeInvariantU(
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
      }

      # Step 2: Reshape to long format (creates within-person deviations)
      long_data <- reshape_long_sim_cr(raw_data)
      long_data <- as.data.frame(long_data)
      long_data <- na.omit(long_data)

      # Step 3: Fit Fixed Individual Effects model
      model_fit_y <- lm(within.y ~ xlagw + ylagw, data = long_data)
      model_fit_x <- lm(within.x ~ xlagw + ylagw, data = long_data)

      # Extract coefficients safely
      get_coef_safe <- function(model, coef_name) {
        coefs <- coefficients(model)
        if (coef_name %in% names(coefs) && !is.na(coefs[coef_name])) {
          return(as.numeric(coefs[coef_name]))
        } else {
          return(NA_real_)
        }
      }

      xlag_x <- get_coef_safe(model_fit_x, "xlagw")   # Autoregressive X effect
      ylag_x <- get_coef_safe(model_fit_x, "ylagw")   # Cross-lagged Y -> X effect
      xlag_y <- get_coef_safe(model_fit_y, "xlagw")   # Cross-lagged X -> Y effect
      ylag_y <- get_coef_safe(model_fit_y, "ylagw")   # Autoregressive Y effect

      # Step 4: Build result row with parameter info
      result_row <- data.frame(
        # Core coefficients (consistent across all dgps)
        xlag_x = xlag_x,
        ylag_x = ylag_x,
        xlag_y = xlag_y,
        ylag_y = ylag_y,
        converged = TRUE,
        n_obs = nrow(long_data),
        trial = j,
        param_combo = i,
        estimator = "olsfi",
        dgp = dgp,
        stringsAsFactors = FALSE
      )

      # Add dgp-specific parameters
      if (dgp == "riclpm") {
        result_row$variance_between_x <- params$variance_between_x
        result_row$variance_between_y <- params$variance_between_y
        result_row$stability_p <- params$stability_p
        result_row$stability_q <- params$stability_q
        result_row$cross_p <- params$cross_p
        result_row$cross_q <- params$cross_q
        result_row$variance_p <- params$variance_p
        result_row$variance_q <- params$variance_q
      } else if (dgp == "clpm") {
        result_row$stability_p <- params$stability_p
        result_row$stability_q <- params$stability_q
        result_row$cross_p <- params$cross_p
        result_row$cross_q <- params$cross_q
        result_row$variance_p <- params$variance_p
        result_row$variance_q <- params$variance_q
      } else if (dgp == "clpmu") {
        result_row$variance_between_x <- params$variance_between_x
        result_row$variance_between_y <- params$variance_between_y
        result_row$stability_p <- params$stability_p
        result_row$stability_q <- params$stability_q
        result_row$cross_p <- params$cross_p
        result_row$cross_q <- params$cross_q
        result_row$variance_p <- params$variance_p
        result_row$variance_q <- params$variance_q
        result_row$confounder_p <- params$confounder_p
        result_row$confounder_q <- params$confounder_q
        result_row$confounder_variance <- params$confounder_variance
        result_row$include_confounder <- params$include_confounder
        result_row$confounder_type <- confounder_type

        # Add confounder_stability only for time_variant
        if (confounder_type == "time_variant") {
          result_row$confounder_stability <- params$confounder_stability
        } else {
          result_row$confounder_stability <- NA
        }
      }

      results[[length(results) + 1]] <- result_row
    }
  }

  # Combine all results
  results_df <- do.call(rbind, results)
  row.names(results_df) <- NULL

  return(results_df)
}
