#' Quick Fix: Replace rowwise() + unnest() with simple loop
run_mc_sims <- function(estimator,
                        param_grid = NULL,
                        trials = 10,
                        waves = 3,
                        sample_size = 1000,
                        verbose = TRUE,
                        data_generation = "riclpm") {

  # CLEAR ENVIRONMENT AT START
  closeAllConnections()
  invisible(gc())

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

  # Valid estimators
  valid_estimators <- c("OLS", "RICLPM", "CLPM", "CTSEM")
  if (!estimator %in% valid_estimators) {
    stop("estimator must be one of: ", paste(valid_estimators, collapse = ", "))
  }

  valid_dgps <- c("clpm", "riclpm", "clpm_confounder")
  data_generation <- match.arg(data_generation, valid_dgps)

  # Create default grid if user does not provide one
  if (is.null(param_grid)) {
    if (data_generation == "clpm_confounder") {
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
        confounder_stability = 0.4
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
  if (data_generation == "clpm_confounder") {
    confounder_defaults <- list(
      confounder_p = 0.3,
      confounder_q = 0.3,
      confounder_variance = 1,
      confounder_stability = 0.4
    )
    default_params <- c(default_params, confounder_defaults)
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
    cat("Parameter combinations:", nrow(param_grid), "\n")
    cat("Parameters being varied:", paste(names(param_grid), collapse = ", "), "\n")
    cat("Data generation process:", data_generation, "\n")
  }

  # Start timing
  tic(paste("Full", estimator, "Monte Carlo Simulation"))

  # SIMPLE APPROACH: Use a basic loop instead of rowwise() + unnest()
  all_results <- list()

  for (i in 1:nrow(param_grid)) {
    if (verbose && i %% max(1, floor(nrow(param_grid)/10)) == 0) {
      cat("Processing parameter combination", i, "of", nrow(param_grid), "\n")
    }

    current_params <- param_grid[i, ]

    # Build parameters for Monte Carlo function
    base_params <- list(
      trials = trials,
      waves = waves,
      stability.p = current_params$stability_p,
      stability.q = current_params$stability_q,
      cross.p = current_params$cross_p,
      cross.q = current_params$cross_q,
      variance.q = current_params$variance_q,
      variance.p = current_params$variance_p,
      variance.between.x = current_params$variance_between_x,
      variance.between.y = current_params$variance_between_y,
      cov.pq = current_params$cov_pq,
      sample_size = sample_size
    )

    # Add estimator-specific parameters
    if (estimator == "OLS") {
      base_params$data_generation <- data_generation
    } else {
      base_params$dgp <- data_generation
    }

    # Add confounder parameters if they exist
    if (data_generation == "clpm_confounder") {
      if ("confounder_p" %in% names(current_params)) {
        base_params$confounder.p <- current_params$confounder_p
      }
      if ("confounder_q" %in% names(current_params)) {
        base_params$confounder.q <- current_params$confounder_q
      }
      if ("confounder_variance" %in% names(current_params)) {
        base_params$confounder.variance <- current_params$confounder_variance
      }
      if ("confounder_stability" %in% names(current_params)) {
        base_params$confounder.stability <- current_params$confounder_stability
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
      }

      # Add parameter combination info to each row of results
      if (is.data.frame(result) && nrow(result) > 0) {
        for (param_name in names(current_params)) {
          if (!param_name %in% names(result)) {
            result[[param_name]] <- current_params[[param_name]]
          }
        }
      }

      all_results[[i]] <- result

    }, error = function(e) {
      if (verbose) {
        cat("Error in parameter combination", i, ":", e$message, "\n")
      }
      all_results[[i]] <- data.frame(
        error_occurred = TRUE,
        error_message = as.character(e$message),
        param_combination = i,
        stringsAsFactors = FALSE
      )
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

  # CLEAR ENVIRONMENT AT END
  closeAllConnections()
  invisible(gc())

  return(final_results)
}
