#' @title run_mc_sims
#' @description Run Monte Carlo simulations across different estimators and parameter combinations
#'
#' This function provides a unified interface for running Monte Carlo simulations across multiple
#' estimators (OLS, RICLPM, CLPM, CTSEM, FI, LCHANGE, RI) with varying parameter combinations.
#' It handles parameter grid expansion, error catching, and result compilation automatically.
#'
#' @param estimator Character string specifying which estimator to use. Must be one of:
#'   \itemize{
#'     \item "OLS" - Ordinary Least Squares estimation
#'     \item "RICLPM" - Random Intercept Cross-Lagged Panel Model
#'     \item "CLPM" - Cross-Lagged Panel Model
#'     \item "CTSEM" - Continuous Time Structural Equation Model
#'     \item "FI" - Fixed Individual effects (First Differences)
#'     \item "LCHANGE" - Latent Change Score Model
#'     \item "RI" - Random Intercepts using lmer
#'   }
#' @param riclpm_type Character string specifying which RICLPM variant to use (only relevant when estimator = "RICLPM").
#'   Must be one of:
#'   \itemize{
#'     \item "riclpm" - Regular RICLPM with autoregressive and cross-lagged effects (default)
#'     \item "riclpm_nolag" - RICLPM without autoregressive effects (cross-lagged only)
#'   }
#' @param param_grid Data frame or NULL. If NULL, uses default parameter grid. If provided,
#'   should contain columns for parameters to vary across simulations. Missing parameters
#'   will be filled with defaults.
#' @param trials Integer. Number of simulation trials to run for each parameter combination.
#'   Default is 10.
#' @param waves Integer. Number of time waves in the simulated data. Default is 3.
#' @param sample_size Integer. Sample size for each simulated dataset. Default is 1000.
#' @param verbose Logical. If TRUE, prints progress information during simulation.
#'   Default is TRUE.
#' @param data_generation Character string specifying the data generation process. Must be one of:
#'   \itemize{
#'     \item "riclpm" - Random Intercept Cross-Lagged Panel Model data generation
#'     \item "clpm" - Cross-Lagged Panel Model data generation
#'     \item "clpmu" - CLPM with unmeasured confounder (uses simCLPMu)
#'   }
#' @param lchange_type Character string specifying latent change model type (only relevant when estimator = "LCHANGE").
#'   Must be one of:
#'   \itemize{
#'     \item "dual_change" - Bivariate latent change model (default)
#'     \item "latent_change" - Single variable latent change model
#'   }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates input parameters
#'   \item Creates or validates the parameter grid
#'   \item Adds default values for missing parameters
#'   \item Loops through each parameter combination
#'   \item Calls the appropriate Monte Carlo function
#'   \item Handles errors and compiles results
#'   \item Returns a combined data frame of all results
#' }
#'
#' Default parameter values:
#' \itemize{
#'   \item stability_p = 0.2 (autoregressive effect for X)
#'   \item stability_q = 0.5 (autoregressive effect for Y)
#'   \item cross_p = 0.0 (cross-lagged effect Y -> X)
#'   \item cross_q = 0.0 (cross-lagged effect X -> Y)
#'   \item variance_q = 0.5 (within-person variance for Y)
#'   \item variance_p = 0.5 (within-person variance for X)
#'   \item variance_between_x = 0.5 (between-person variance for X)
#'   \item variance_between_y = 0.5 (between-person variance for Y)
#'   \item cov_pq = 0 (within-time covariance between X and Y)
#' }
#'
#' When data_generation = "clpmu", additional confounder parameters:
#' \itemize{
#'   \item confounder_p = 0.3 (effect of confounder on X)
#'   \item confounder_q = 0.3 (effect of confounder on Y)
#'   \item confounder_variance = 1 (variance of confounder)
#'   \item confounder_stability = 0.4 (autoregressive effect of confounder)
#'   \item include_confounder = TRUE (whether to include confounder)
#' }
#'
#' @return A data frame containing simulation results with the following columns:
#' \itemize{
#'   \item Parameter values used in simulation
#'   \item Estimated coefficients (varies by estimator)
#'   \item trial - Trial number
#'   \item estimator - Which estimator was used
#'   \item dgp - Data generation process used
#'   \item error_occurred - TRUE if an error occurred in that trial
#'   \item error_message - Error message if error_occurred is TRUE
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' results <- run_mc_sims(
#'   estimator = "RICLPM",
#'   trials = 5,
#'   waves = 3,
#'   sample_size = 500
#' )
#'
#' # Test Random Intercepts with both full and cross-lagged only
#' results_ri_full <- run_mc_sims(
#'   estimator = "RI",
#'   param_grid = data.frame(include_lagged_dv = TRUE),
#'   trials = 10,
#'   waves = 4,
#'   sample_size = 1000
#' )
#'
#' results_ri_cross <- run_mc_sims(
#'   estimator = "RI",
#'   param_grid = data.frame(include_lagged_dv = FALSE),
#'   trials = 10,
#'   waves = 4,
#'   sample_size = 1000
#' )
#'
#' # Compare regular RICLPM vs no-lag RICLPM
#' results_regular <- run_mc_sims(
#'   estimator = "RICLPM",
#'   riclpm_type = "riclpm",
#'   trials = 10,
#'   waves = 4,
#'   sample_size = 1000
#' )
#'
#' results_nolag <- run_mc_sims(
#'   estimator = "RICLPM",
#'   riclpm_type = "riclpm_nolag",
#'   trials = 10,
#'   waves = 4,
#'   sample_size = 1000
#' )
#'
#' # Compare latent change models
#' results_dual_change <- run_mc_sims(
#'   estimator = "LCHANGE",
#'   lchange_type = "dual_change",
#'   trials = 10,
#'   waves = 4,
#'   sample_size = 1000
#' )
#'
#' results_single_change <- run_mc_sims(
#'   estimator = "LCHANGE",
#'   lchange_type = "latent_change",
#'   trials = 10,
#'   waves = 4,
#'   sample_size = 1000
#' )
#'
#' # Study bias from unmeasured confounder
#' results_confounded <- run_mc_sims(
#'   estimator = "CLPM",
#'   trials = 10,
#'   waves = 4,
#'   sample_size = 1000,
#'   data_generation = "clpmu"  # Data has confounder, but CLPM ignores it
#' )
#'
#' # Custom parameter grid
#' my_grid <- expand.grid(
#'   stability_p = c(0.2, 0.5),
#'   stability_q = c(0.3, 0.6),
#'   cross_p = c(0.0, 0.1),
#'   cross_q = c(0.0, 0.1),
#'   variance_between_x = c(0.5, 1.0)
#' )
#'
#' results <- run_mc_sims(
#'   estimator = "CLPM",
#'   param_grid = my_grid,
#'   trials = 10,
#'   waves = 4,
#'   sample_size = 1000,
#'   data_generation = "riclpm"
#' )
#'
#' # Compare estimators under confounding
#' estimators <- c("OLS", "RICLPM", "CLPM", "RI")
#' all_results <- lapply(estimators, function(est) {
#'   run_mc_sims(estimator = est, trials = 5, waves = 3, data_generation = "clpmu")
#' })
#' combined_results <- do.call(rbind, all_results)
#' }
#'
#' @import dplyr tictoc
#' @export
run_mc_sims <- function(estimator,
                        riclpm_type = "riclpm",
                        param_grid = NULL,
                        lchange_type = "dual_change",
                        trials = 10,
                        waves = 3,
                        sample_size = 1000,
                        verbose = TRUE,
                        data_generation = "riclpm") {

  # Clear any lingering lavaan objects
  if (exists(".lavaan_cache", envir = .GlobalEnv)) {
    rm(.lavaan_cache, envir = .GlobalEnv)
  }

  # Clear temporary objects that might conflict
  temp_objects <- ls(pattern = "^(tmp_|temp_|model_|fit_)", envir = .GlobalEnv)
  if (length(temp_objects) > 0) {
    rm(list = temp_objects, envir = .GlobalEnv)
  }

  library(dplyr)
  library(tictoc)

  # Valid estimators - UPDATED TO INCLUDE "RI"
  valid_estimators <- c("OLS", "RICLPM", "CLPM", "CTSEM", "FI", "LCHANGE", "RI")
  if (!estimator %in% valid_estimators) {
    stop("estimator must be one of: ", paste(valid_estimators, collapse = ", "))
  }

  # Valid RICLPM types
  valid_riclpm_types <- c("riclpm", "riclpm_nolag")
  if (!riclpm_type %in% valid_riclpm_types) {
    stop("riclpm_type must be one of: ", paste(valid_riclpm_types, collapse = ", "))
  }

  # Valid data generation processes
  valid_dgps <- c("clpm", "riclpm", "clpmu")
  data_generation <- match.arg(data_generation, valid_dgps)

  # Create default grid if user does not provide one
  if (is.null(param_grid)) {
    if (data_generation == "clpmu") {
      param_grid <- expand.grid(
        stability_q = seq(0.2, 0.7, by = 0.1),
        variance_between_x = seq(0.3, 1, by = 0.1),
        stability_p = 0.2,
        cross_p = 0.0,
        cross_q = 0.0,
        variance_q = 0.5,
        variance_p = 0.5,
        variance_between_y = 0.5,
        cov_pq = 0,
        confounder_p = 0.3,
        confounder_q = 0.3,
        confounder_variance = 1,
        confounder_stability = 0.4,
        include_confounder = TRUE
      )
    } else {
      param_grid <- expand.grid(
        stability_q = seq(0.2, 0.7, by = 0.1),
        variance_between_x = seq(0.3, 1, by = 0.1),
        stability_p = 0.2,
        cross_p = 0.0,
        cross_q = 0.0,
        variance_q = 0.5,
        variance_p = 0.5,
        variance_between_y = 0.5,
        cov_pq = 0
      )
    }
    if (verbose) cat("Using default parameters\n")
  }

  # Convert to data frame if it's a list
  if (is.list(param_grid) && !is.data.frame(param_grid)) {
    param_grid <- as.data.frame(param_grid)
  }

  # Ensure param_grid is a proper data frame
  param_grid <- as.data.frame(param_grid)

  # Set default parameters
  default_params <- list(
    stability_p = 0.2,
    stability_q = 0.5,
    cross_p = 0.0,
    cross_q = 0.0,
    variance_q = 0.5,
    variance_p = 0.5,
    variance_between_x = 0.5,
    variance_between_y = 0.5,
    cov_pq = 0
  )

  # Add confounder defaults if using confounder DGP
  if (data_generation == "clpmu") {
    confounder_defaults <- list(
      confounder_p = 0.3,
      confounder_q = 0.3,
      confounder_variance = 1,
      confounder_stability = 0.4,
      include_confounder = TRUE
    )
    default_params <- c(default_params, confounder_defaults)
  }

  # Add RI-specific defaults
  if (estimator == "RI") {
    ri_defaults <- list(
      include_lagged_dv = TRUE
    )
    default_params <- c(default_params, ri_defaults)
  }

  # Add missing parameters with default values
  for (param_name in names(default_params)) {
    if (!param_name %in% names(param_grid)) {
      param_grid[[param_name]] <- default_params[[param_name]]
      if (verbose) {
        cat("Added default value for", param_name, "=", default_params[[param_name]], "\n")
      }
    }
  }

  if (verbose) {
    cat("Running", estimator, "Monte Carlo simulation\n")
    if (estimator == "RICLPM") {
      cat("RICLPM type:", riclpm_type, "\n")
    }
    if (estimator == "LCHANGE") {
      cat("Latent change type:", lchange_type, "\n")
    }
    cat("Parameter combinations:", nrow(param_grid), "\n")
    cat("Parameters being varied:", paste(names(param_grid), collapse = ", "), "\n")
    cat("Data generation process:", data_generation, "\n")
  }

  # Start timing - disable tictoc auto-saving to prevent gzfile errors
  options(tictoc.save = FALSE)
  tic(paste("Full", estimator, "Monte Carlo Simulation"))

  # SIMPLE APPROACH: Use a basic loop instead of rowwise() + unnest()
  all_results <- list()

  for (i in 1:nrow(param_grid)) {
    if (verbose && i %% max(1, floor(nrow(param_grid)/10)) == 0) {
      cat("Processing parameter combination", i, "of", nrow(param_grid), "\n")
    }

    # FIX: Convert to regular list to avoid S4 class issues
    current_params <- as.list(param_grid[i, ])

    # Build parameters for Monte Carlo function
    base_params <- list(
      trials = trials,
      waves = waves,
      stability_p = current_params[["stability_p"]],
      stability_q = current_params[["stability_q"]],
      cross_p = current_params[["cross_p"]],
      cross_q = current_params[["cross_q"]],
      variance_q = current_params[["variance_q"]],
      variance_p = current_params[["variance_p"]],
      variance_between_x = current_params[["variance_between_x"]],
      variance_between_y = current_params[["variance_between_y"]],
      cov_pq = current_params[["cov_pq"]],
      sample_size = sample_size
    )

    # Set dgp parameter based on data generation method
    if (data_generation == "clpmu") {
      # For confounded data, we need to modify the Monte Carlo functions
      # to use simCLPMu instead of their normal data generation
      base_params$dgp <- "clpmu"

      # Add confounder parameters
      if ("confounder_p" %in% names(current_params)) {
        base_params$confounder_p <- current_params[["confounder_p"]]
      }
      if ("confounder_q" %in% names(current_params)) {
        base_params$confounder_q <- current_params[["confounder_q"]]
      }
      if ("confounder_variance" %in% names(current_params)) {
        base_params$confounder_variance <- current_params[["confounder_variance"]]
      }
      if ("confounder_stability" %in% names(current_params)) {
        base_params$confounder_stability <- current_params[["confounder_stability"]]
      }
      if ("include_confounder" %in% names(current_params)) {
        base_params$include_confounder <- current_params[["include_confounder"]]
      }
    } else {
      base_params$dgp <- data_generation
    }

    # Add RICLPM-specific parameter
    if (estimator == "RICLPM") {
      base_params$estimator <- riclpm_type
    }

    # Add LCHANGE-specific parameter
    if (estimator == "LCHANGE") {
      base_params$model_type <- lchange_type
    }

    # Add RI-specific parameter
    if (estimator == "RI") {
      if ("include_lagged_dv" %in% names(current_params)) {
        base_params$include_lagged_dv <- current_params[["include_lagged_dv"]]
      }
    }

    # Call the appropriate Monte Carlo function
    tryCatch({
      if (estimator == "OLS") {
        result <- do.call(monteCarloOLS, base_params)
      } else if (estimator == "RICLPM") {
        result <- do.call(monteCarloRICLPM, base_params)
      } else if (estimator == "CLPM") {
        result <- do.call(monteCarloCLPM, base_params)
      } else if (estimator == "CTSEM") {
        safe_mc_wrapped <- function(...) suppressMessages(suppressWarnings(monteCarloCTSEM(...)))
        result <- do.call(safe_mc_wrapped, base_params)
      } else if (estimator == "FI") {
        result <- do.call(monteCarloFI, base_params)
      } else if (estimator == "LCHANGE") {
        result <- do.call(monteCarloLChange, base_params)
      } else if (estimator == "RI") {
        result <- do.call(monteCarloRI, base_params)
      }

      # Add parameter combination info to each row of results
      if (is.data.frame(result) && nrow(result) > 0) {
        for (param_name in names(current_params)) {
          if (!param_name %in% names(result)) {
            result[[param_name]] <- current_params[[param_name]]
          }
        }

        # Add riclpm_type to results for RICLPM runs
        if (estimator == "RICLPM" && !"riclpm_type" %in% names(result)) {
          result$riclpm_type <- riclpm_type
        }

        # Add lchange_type to results for LCHANGE runs
        if (estimator == "LCHANGE" && !"lchange_type" %in% names(result)) {
          result$lchange_type <- lchange_type
        }
      }

      all_results[[i]] <- result

    }, error = function(e) {
      if (verbose) {
        cat("Error in parameter combination", i, ":", e$message, "\n")
      }

      # Create error result with parameter info
      error_result <- data.frame(
        error_occurred = TRUE,
        error_message = as.character(e$message),
        param_combination = i,
        stringsAsFactors = FALSE
      )

      # Add parameter values to error result
      for (param_name in names(current_params)) {
        error_result[[param_name]] <- current_params[[param_name]]
      }

      if (estimator == "RICLPM") {
        error_result$riclpm_type <- riclpm_type
      }

      if (estimator == "LCHANGE") {
        error_result$lchange_type <- lchange_type
      }

      all_results[[i]] <- error_result
    })
  }

  # Combine all results using simple rbind
  if (verbose) cat("Combining results...\n")

  valid_results <- all_results[!sapply(all_results, is.null)]

  if (length(valid_results) > 0) {
    final_results <- do.call(rbind, valid_results)
  } else {
    final_results <- data.frame(
      error_occurred = TRUE,
      error_message = "No valid results obtained",
      stringsAsFactors = FALSE
    )
  }

  # Count and report errors
  if ("error_occurred" %in% names(final_results)) {
    error_count <- sum(final_results$error_occurred, na.rm = TRUE)
    if (error_count > 0 && verbose) {
      cat("Warning:", error_count, "trials resulted in errors\n")
    }
  }

  toc()

  return(final_results)
}
