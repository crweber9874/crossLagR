#' @title monteCarloLChange
#' @description Monte Carlo simulation for Latent Change Score Models
#' This function combines the estimation and simulation of latent change models over a fixed number of trials.
#' The output is a data frame that includes the estimated coefficients across these trials
#' @param trials The number of trials for the Monte Carlo simulation.
#' @param waves The number of waves (time points) in the model.
#' @param dgp Specify the data generation method: "riclpm", "clpm", "clpmu"
#' @param model_type Type of latent change model: "latent_change" (single variable) or "dual_change" (bivariate)
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
#' @param confounder_type Character string specifying confounder type for clpmu dgp.
#'   Must be one of "time_variant" (default) or "time_invariant".
#' @param verbose Whether to print progress and error messages
#'
#' @import tidyr dplyr lavaan
#' @importFrom dplyr %>% select mutate row_number group_by ungroup
#'
#' @return A data frame containing the results of the simulation.
#'
#' @export
monteCarloLChange <- function(
    trials = 10,
    waves = 5,
    dgp = "riclpm",
    model_type = "dual_change",
    stability_q = 0.25,
    stability_p = 0.25,
    variance_between_y = 0.5,
    variance_between_x = 0.5,
    cross_p = 0,
    cross_q = 0,
    variance_p = 1,
    variance_q = 1,
    sample_size = 2500,
    # Confounder parameters
    confounder_p = 0.3,
    confounder_q = 0.3,
    confounder_variance = 1,
    confounder_stability = 0.4,
    include_confounder = TRUE,
    confounder_type = "time_variant",  # NEW PARAMETER
    verbose = FALSE,
    ...
) {
  library(lavaan)

  # Validate inputs
  valid_dgps <- c("riclpm", "clpm", "clpmu")
  valid_model_types <- c("latent_change", "dual_change")

  if (!dgp %in% valid_dgps) {
    stop("dgp must be one of: ", paste(valid_dgps, collapse = ", "))
  }

  if (!model_type %in% valid_model_types) {
    stop("model_type must be one of: ", paste(valid_model_types, collapse = ", "))
  }

  # Validate confounder_type
  if (!confounder_type %in% c("time_variant", "time_invariant")) {
    stop("confounder_type must be either 'time_variant' or 'time_invariant'")
  }

  # Safe coefficient extraction function
  get_coef_safe <- function(coef_name, coeffs) {
    if (coef_name %in% names(coeffs)) {
      return(as.numeric(coeffs[coef_name]))
    } else {
      return(NA_real_)
    }
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
      variance_q = variance_q
    ) %>% as.data.frame()

    results <- list()

    for (i in 1:nrow(simulation_parameters)) {
      for (j in 1:trials) {
        params <- simulation_parameters[i, ]

        # Base result row with parameter info
        base_result <- data.frame(
          variance_between_x = params$variance_between_x,
          variance_between_y = params$variance_between_y,
          stability_p = params$stability_p,
          stability_q = params$stability_q,
          cross_p = params$cross_p,
          cross_q = params$cross_q,
          variance_p = params$variance_p,
          variance_q = params$variance_q,
          trial = j,
          param_combo = i,
          estimator = "lchange",
          model_type = model_type,
          dgp = "riclpm"
        )

        # Try to generate data and fit model
        model_result <- tryCatch({
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

          # Generate model syntax
          model_syntax <- latentChange(
            waves = waves,
            model_type = model_type
          )

          # Fit model
          model <- lavaan::lavaan(model_syntax, data = dat)

          # Check convergence
          if (!lavInspect(model, "converged")) {
            stop("Model did not converge")
          }

          # Extract coefficients based on model type
          coeffs <- coef(model)

          if (model_type == "dual_change") {
            # For dual change model, extract key parameters
            change_x_to_y <- get_coef_safe("change.y", coeffs)  # x change → y change
            change_y_to_x <- get_coef_safe("change.x", coeffs)  # y change → x change
            coupling_x_to_y <- get_coef_safe("cl_x", coeffs)    # x level → y change
            coupling_y_to_x <- get_coef_safe("cl_y", coeffs)    # y level → x change
            ar_x <- get_coef_safe("ar_x", coeffs)
            ar_y <- get_coef_safe("ar_y", coeffs)
            phi_x <- get_coef_safe("phi_x", coeffs)
            phi_y <- get_coef_safe("phi_y", coeffs)

            # Map to standard naming convention for comparison with other estimators
            xlag_x <- ar_x        # Not exactly autoregressive, but similar concept
            ylag_y <- ar_y        # Not exactly autoregressive, but similar concept
            xlag_y <- change_x_to_y  # Change-to-change cross-lagged
            ylag_x <- change_y_to_x  # Change-to-change cross-lagged

          } else {
            # For single variable model, fewer parameters
            xlag_x <- NA_real_
            ylag_y <- NA_real_
            xlag_y <- NA_real_
            ylag_x <- NA_real_
            change_x_to_y <- NA_real_
            change_y_to_x <- NA_real_
            coupling_x_to_y <- NA_real_
            coupling_y_to_x <- NA_real_
            ar_x <- NA_real_
            ar_y <- NA_real_
            phi_x <- NA_real_
            phi_y <- NA_real_
          }

          list(
            success = TRUE,
            xlag_x = xlag_x,
            ylag_x = ylag_x,
            xlag_y = xlag_y,
            ylag_y = ylag_y,
            change_x_to_y = change_x_to_y,
            change_y_to_x = change_y_to_x,
            coupling_x_to_y = coupling_x_to_y,
            coupling_y_to_x = coupling_y_to_x,
            ar_x = ar_x,
            ar_y = ar_y,
            phi_x = phi_x,
            phi_y = phi_y,
            converged = TRUE,
            error_message = NA_character_,
            error_type = NA_character_
          )

        }, error = function(e) {
          error_msg <- as.character(e$message)
          error_type <- dplyr::case_when(
            grepl("lav_start_check_cov", error_msg) ~ "covariance_start_values",
            grepl("converge", error_msg) ~ "convergence_failure",
            grepl("singular", error_msg) ~ "singular_matrix",
            grepl("identification", error_msg) ~ "identification_problem",
            TRUE ~ "other_error"
          )

          if (verbose) {
            cat("Error in trial", j, "param combo", i, ":", error_msg, "\n")
          }

          list(
            success = FALSE,
            xlag_x = NA_real_,
            ylag_x = NA_real_,
            xlag_y = NA_real_,
            ylag_y = NA_real_,
            change_x_to_y = NA_real_,
            change_y_to_x = NA_real_,
            coupling_x_to_y = NA_real_,
            coupling_y_to_x = NA_real_,
            ar_x = NA_real_,
            ar_y = NA_real_,
            phi_x = NA_real_,
            phi_y = NA_real_,
            converged = FALSE,
            error_message = error_msg,
            error_type = error_type
          )
        })

        # Combine base result with model results
        final_result <- cbind(base_result, model_result[c(
          "xlag_x", "ylag_x", "xlag_y", "ylag_y",
          "change_x_to_y", "change_y_to_x",
          "coupling_x_to_y", "coupling_y_to_x",
          "ar_x", "ar_y", "phi_x", "phi_y",
          "converged", "error_message", "error_type"
        )])

        results[[length(results) + 1]] <- final_result
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
      variance_q = variance_q
    ) %>% as.data.frame()

    results <- list()

    for (i in 1:nrow(simulation_parameters)) {
      for (j in 1:trials) {
        params <- simulation_parameters[i, ]

        # Base result row
        base_result <- data.frame(
          stability_p = params$stability_p,
          stability_q = params$stability_q,
          cross_p = params$cross_p,
          cross_q = params$cross_q,
          variance_p = params$variance_p,
          variance_q = params$variance_q,
          trial = j,
          param_combo = i,
          estimator = "lchange",
          model_type = model_type,
          dgp = "clpm"
        )

        # Try to generate data and fit model
        model_result <- tryCatch({
          # Generate data
          dat <- simCLPM(
            waves = waves,
            stability_p = params$stability_p,
            stability_q = params$stability_q,
            cross_p = params$cross_p,
            cross_q = params$cross_q,
            variance_p = params$variance_p,
            variance_q = params$variance_q,
            sample.nobs = sample_size
          )$data

          # Generate model syntax
          model_syntax <- latentChange(
            waves = waves,
            model_type = model_type
          )

          # Fit model
          model <- lavaan::lavaan(model_syntax, data = dat)

          # Check convergence
          if (!lavInspect(model, "converged")) {
            stop("Model did not converge")
          }

          # Extract coefficients based on model type
          coeffs <- coef(model)

          if (model_type == "dual_change") {
            # For dual change model, extract key parameters
            change_x_to_y <- get_coef_safe("change.y", coeffs)  # x change → y change
            change_y_to_x <- get_coef_safe("change.x", coeffs)  # y change → x change
            coupling_x_to_y <- get_coef_safe("cl_x", coeffs)    # x level → y change
            coupling_y_to_x <- get_coef_safe("cl_y", coeffs)    # y level → x change
            ar_x <- get_coef_safe("ar_x", coeffs)
            ar_y <- get_coef_safe("ar_y", coeffs)
            phi_x <- get_coef_safe("phi_x", coeffs)
            phi_y <- get_coef_safe("phi_y", coeffs)

            # Map to standard naming convention for comparison with other estimators
            xlag_x <- ar_x        # Not exactly autoregressive, but similar concept
            ylag_y <- ar_y        # Not exactly autoregressive, but similar concept
            xlag_y <- change_x_to_y  # Change-to-change cross-lagged
            ylag_x <- change_y_to_x  # Change-to-change cross-lagged

          } else {
            # For single variable model, fewer parameters
            xlag_x <- NA_real_
            ylag_y <- NA_real_
            xlag_y <- NA_real_
            ylag_x <- NA_real_
            change_x_to_y <- NA_real_
            change_y_to_x <- NA_real_
            coupling_x_to_y <- NA_real_
            coupling_y_to_x <- NA_real_
            ar_x <- NA_real_
            ar_y <- NA_real_
            phi_x <- NA_real_
            phi_y <- NA_real_
          }

          list(
            success = TRUE,
            xlag_x = xlag_x,
            ylag_x = ylag_x,
            xlag_y = xlag_y,
            ylag_y = ylag_y,
            change_x_to_y = change_x_to_y,
            change_y_to_x = change_y_to_x,
            coupling_x_to_y = coupling_x_to_y,
            coupling_y_to_x = coupling_y_to_x,
            ar_x = ar_x,
            ar_y = ar_y,
            phi_x = phi_x,
            phi_y = phi_y,
            converged = TRUE,
            error_message = NA_character_,
            error_type = NA_character_
          )

        }, error = function(e) {
          error_msg <- as.character(e$message)
          error_type <- dplyr::case_when(
            grepl("lav_start_check_cov", error_msg) ~ "covariance_start_values",
            grepl("converge", error_msg) ~ "convergence_failure",
            grepl("singular", error_msg) ~ "singular_matrix",
            grepl("identification", error_msg) ~ "identification_problem",
            TRUE ~ "other_error"
          )

          if (verbose) {
            cat("Error in trial", j, "param combo", i, ":", error_msg, "\n")
          }

          list(
            success = FALSE,
            xlag_x = NA_real_,
            ylag_x = NA_real_,
            xlag_y = NA_real_,
            ylag_y = NA_real_,
            change_x_to_y = NA_real_,
            change_y_to_x = NA_real_,
            coupling_x_to_y = NA_real_,
            coupling_y_to_x = NA_real_,
            ar_x = NA_real_,
            ar_y = NA_real_,
            phi_x = NA_real_,
            phi_y = NA_real_,
            converged = FALSE,
            error_message = error_msg,
            error_type = error_type
          )
        })

        # Combine base result with model results
        final_result <- cbind(base_result, model_result[c(
          "xlag_x", "ylag_x", "xlag_y", "ylag_y",
          "change_x_to_y", "change_y_to_x",
          "coupling_x_to_y", "coupling_y_to_x",
          "ar_x", "ar_y", "phi_x", "phi_y",
          "converged", "error_message", "error_type"
        )])

        results[[length(results) + 1]] <- final_result
      }
    }

    results_df <- do.call(rbind, results)
  }

  else if(dgp == "clpmu") {  # COMPLETE SECTION FOR CONFOUNDER MODEL
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

    results <- list()

    for (i in 1:nrow(simulation_parameters)) {
      for (j in 1:trials) {
        params <- simulation_parameters[i, ]

        # Base result row
        base_result <- data.frame(
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
          trial = j,
          param_combo = i,
          estimator = "lchange",
          model_type = model_type,
          dgp = "clpmu"
        )

        # Add confounder_stability only for time_variant confounders
        if (confounder_type == "time_variant") {
          base_result$confounder_stability <- params$confounder_stability
        } else {
          base_result$confounder_stability <- NA
        }

        # Try to generate data and fit model
        model_result <- tryCatch({
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

          # Generate model syntax
          model_syntax <- latentChange(
            waves = waves,
            model_type = model_type
          )

          # Fit model
          model <- lavaan::lavaan(model_syntax, data = dat)

          # Check convergence
          if (!lavInspect(model, "converged")) {
            stop("Model did not converge")
          }

          # Extract coefficients based on model type
          coeffs <- coef(model)

          if (model_type == "dual_change") {
            # For dual change model, extract key parameters
            change_x_to_y <- get_coef_safe("change.y", coeffs)  # x change → y change
            change_y_to_x <- get_coef_safe("change.x", coeffs)  # y change → x change
            coupling_x_to_y <- get_coef_safe("cl_x", coeffs)    # x level → y change
            coupling_y_to_x <- get_coef_safe("cl_y", coeffs)    # y level → x change
            ar_x <- get_coef_safe("ar_x", coeffs)
            ar_y <- get_coef_safe("ar_y", coeffs)
            phi_x <- get_coef_safe("phi_x", coeffs)
            phi_y <- get_coef_safe("phi_y", coeffs)

            # Map to standard naming convention for comparison with other estimators
            xlag_x <- ar_x        # Not exactly autoregressive, but similar concept
            ylag_y <- ar_y        # Not exactly autoregressive, but similar concept
            xlag_y <- change_x_to_y  # Change-to-change cross-lagged
            ylag_x <- change_y_to_x  # Change-to-change cross-lagged

          } else {
            # For single variable model, fewer parameters
            xlag_x <- NA_real_
            ylag_y <- NA_real_
            xlag_y <- NA_real_
            ylag_x <- NA_real_
            change_x_to_y <- NA_real_
            change_y_to_x <- NA_real_
            coupling_x_to_y <- NA_real_
            coupling_y_to_x <- NA_real_
            ar_x <- NA_real_
            ar_y <- NA_real_
            phi_x <- NA_real_
            phi_y <- NA_real_
          }

          list(
            success = TRUE,
            xlag_x = xlag_x,
            ylag_x = ylag_x,
            xlag_y = xlag_y,
            ylag_y = ylag_y,
            change_x_to_y = change_x_to_y,
            change_y_to_x = change_y_to_x,
            coupling_x_to_y = coupling_x_to_y,
            coupling_y_to_x = coupling_y_to_x,
            ar_x = ar_x,
            ar_y = ar_y,
            phi_x = phi_x,
            phi_y = phi_y,
            converged = TRUE,
            error_message = NA_character_,
            error_type = NA_character_
          )

        }, error = function(e) {
          error_msg <- as.character(e$message)
          error_type <- dplyr::case_when(
            grepl("lav_start_check_cov", error_msg) ~ "covariance_start_values",
            grepl("converge", error_msg) ~ "convergence_failure",
            grepl("singular", error_msg) ~ "singular_matrix",
            grepl("identification", error_msg) ~ "identification_problem",
            TRUE ~ "other_error"
          )

          if (verbose) {
            cat("Error in trial", j, "param combo", i, ":", error_msg, "\n")
          }

          list(
            success = FALSE,
            xlag_x = NA_real_,
            ylag_x = NA_real_,
            xlag_y = NA_real_,
            ylag_y = NA_real_,
            change_x_to_y = NA_real_,
            change_y_to_x = NA_real_,
            coupling_x_to_y = NA_real_,
            coupling_y_to_x = NA_real_,
            ar_x = NA_real_,
            ar_y = NA_real_,
            phi_x = NA_real_,
            phi_y = NA_real_,
            converged = FALSE,
            error_message = error_msg,
            error_type = error_type
          )
        })

        # Combine base result with model results
        final_result <- cbind(base_result, model_result[c(
          "xlag_x", "ylag_x", "xlag_y", "ylag_y",
          "change_x_to_y", "change_y_to_x",
          "coupling_x_to_y", "coupling_y_to_x",
          "ar_x", "ar_y", "phi_x", "phi_y",
          "converged", "error_message", "error_type"
        )])

        results[[length(results) + 1]] <- final_result
      }
    }

    results_df <- do.call(rbind, results)
  }

  # Add success indicator
  results_df$success <- results_df$converged & is.na(results_df$error_message)

  return(results_df)
}
