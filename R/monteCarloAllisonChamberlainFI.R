#' @title monteCarloAllisonChamberlainFI
#' @description Monte Carlo simulation for the Allison-Chamberlain Fixed-Effects SEM estimator.
#'
#' This function combines data simulation and estimation using the Allison-Chamberlain
#' fixed-effects SEM approach over a fixed number of trials. It generates data from
#' either a CLPM, RI-CLPM, or CLPM with unmeasured confounder DGP, fits the
#' Allison-Chamberlain model, and extracts cross-lagged and autoregressive coefficients.
#'
#' @param trials Integer. Number of Monte Carlo trials. Default is 10.
#' @param waves Integer. Number of waves (time points). Must be >= 3. Default is 5.
#' @param stability_p Numeric. Autoregressive effect for x variable. Default is 0.25.
#' @param stability_q Numeric. Autoregressive effect for y variable. Default is 0.25.
#' @param cross_p Numeric. Cross-lagged effect of y on x. Default is 0.
#' @param cross_q Numeric. Cross-lagged effect of x on y. Default is 0.
#' @param variance_p Numeric. Within-person variance for x. Default is 1.
#' @param variance_q Numeric. Within-person variance for y. Default is 1.
#' @param sample_size Integer. Sample size for each simulation. Default is 2500.
#' @param dgp Character. Data generation process: \code{"clpm"}, \code{"riclpm"}, or
#'   \code{"clpmu"}. Default is \code{"clpm"}.
#' @param variance_between_x Numeric. Between-person variance for x (RI-CLPM DGP). Default is 0.5.
#' @param variance_between_y Numeric. Between-person variance for y (RI-CLPM DGP). Default is 0.5.
#' @param confounder_p Numeric. Effect of confounder on x (CLPMU DGP). Default is 0.3.
#' @param confounder_q Numeric. Effect of confounder on y (CLPMU DGP). Default is 0.3.
#' @param confounder_variance Numeric. Variance of confounder (CLPMU DGP). Default is 1.
#' @param confounder_stability Numeric. AR effect of confounder (CLPMU DGP). Default is 0.4.
#' @param include_confounder Logical. Include confounder in CLPMU model. Default is TRUE.
#' @param confounder_type Character. \code{"time_variant"} or \code{"time_invariant"}. Default is \code{"time_variant"}.
#' @param cov_pq Numeric. Covariance between p and q innovations. Default is 0.1.
#' @param verbose Logical. Print progress messages. Default is FALSE.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame containing estimated coefficients (\code{xlag_x}, \code{ylag_x},
#'   \code{xlag_y}, \code{ylag_y}), convergence status, and error/warning information
#'   for each trial and parameter combination.
#'
#' @details
#' The function follows the same structure as other Monte Carlo functions in the package
#' (\code{monteCarloCLPM}, \code{monteCarloRICLPM}). It generates data using the specified
#' DGP, fits the Allison-Chamberlain model via \code{lavaan::sem()}, and extracts the
#' labeled regression coefficients.
#'
#' @examples
#' \dontrun{
#' # Basic Monte Carlo with CLPM DGP
#' results <- monteCarloAllisonChamberlainFI(trials = 5, waves = 4, sample_size = 500)
#' colMeans(results[, c("xlag_x", "ylag_x", "xlag_y", "ylag_y")], na.rm = TRUE)
#'
#' # Test with RI-CLPM DGP (data has fixed effects, estimator should capture them)
#' results_ri <- monteCarloAllisonChamberlainFI(
#'     trials = 10, waves = 5, dgp = "riclpm",
#'     variance_between_x = 1, variance_between_y = 1
#' )
#' }
#'
#' @import lavaan
#' @importFrom dplyr %>% case_when
#' @export
monteCarloAllisonChamberlainFI <- function(
  trials = 10,
  waves = 5,
  stability_p = 0.25,
  stability_q = 0.25,
  cross_p = 0,
  cross_q = 0,
  variance_p = 1,
  variance_q = 1,
  sample_size = 2500,
  dgp = "clpm",
  variance_between_x = 0.5,
  variance_between_y = 0.5,
  # Confounder parameters
  confounder_p = 0.3,
  confounder_q = 0.3,
  confounder_variance = 1,
  confounder_stability = 0.4,
  include_confounder = TRUE,
  confounder_type = "time_variant",
  cov_pq = 0.1,
  verbose = FALSE,
  ...
) {
    library(lavaan)

    # Validate inputs
    valid_dgps <- c("riclpm", "clpm", "clpmu")
    if (!dgp %in% valid_dgps) {
        stop("dgp must be one of: ", paste(valid_dgps, collapse = ", "))
    }
    if (waves < 3) {
        stop("Allison-Chamberlain FI requires at least 3 waves.")
    }
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

    # Generate the Allison-Chamberlain model syntax
    model_syntax <- estimateAllisonChamberlainFI(waves = waves)

    # ==================== CLPM DGP ====================
    if (dgp == "clpm") {
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

                base_result <- data.frame(
                    stability_p = params$stability_p,
                    stability_q = params$stability_q,
                    cross_p = params$cross_p,
                    cross_q = params$cross_q,
                    variance_p = params$variance_p,
                    variance_q = params$variance_q,
                    trial = j,
                    param_combo = i,
                    estimator = "allison_chamberlain_fi",
                    dgp = "clpm"
                )

                model_result <- tryCatch(
                    {
                        warnings_collected <- c()

                        withCallingHandlers(
                            {
                                dat <- simCLPM(
                                    waves = waves,
                                    beta_x = params$stability_p,
                                    beta_y = params$stability_q,
                                    omega_yx = params$cross_p,
                                    omega_xy = params$cross_q,
                                    var_x = params$variance_p,
                                    var_y = params$variance_q,
                                    cov_xy = cov_pq,
                                    sample_size = sample_size
                                )$data

                                # Fit Allison-Chamberlain model
                                model <- lavaan::sem(model_syntax, data = dat)

                                if (!lavInspect(model, "converged")) {
                                    stop("Model did not converge")
                                }

                                # Extract coefficients
                                param_table <- lavaan::parameterEstimates(model)
                                beta_x <- param_table[param_table$label == "ar_x", "est"]
                                beta_y <- param_table[param_table$label == "ar_y", "est"]
                                omega_xy <- param_table[param_table$label == "cl_xy", "est"]
                                omega_yx <- param_table[param_table$label == "cl_yx", "est"]
                            },
                            warning = function(w) {
                                warnings_collected <<- c(warnings_collected, w$message)
                                invokeRestart("muffleWarning")
                            }
                        )

                        list(
                            success = TRUE,
                            xlag_x = as.numeric(beta_x[1]),
                            ylag_x = as.numeric(omega_yx[1]),
                            xlag_y = as.numeric(omega_xy[1]),
                            ylag_y = as.numeric(beta_y[1]),
                            converged = TRUE,
                            error_message = NA_character_,
                            error_type = NA_character_,
                            warning_messages = if (length(warnings_collected) > 0) paste(warnings_collected, collapse = "; ") else NA_character_,
                            warning_count = length(warnings_collected)
                        )
                    },
                    error = function(e) {
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
                            converged = FALSE,
                            error_message = error_msg,
                            error_type = error_type,
                            warning_messages = NA_character_,
                            warning_count = 0
                        )
                    }
                )

                final_result <- cbind(base_result, model_result[c(
                    "xlag_x", "ylag_x", "xlag_y", "ylag_y",
                    "converged", "error_message", "error_type", "warning_messages", "warning_count"
                )])

                results[[length(results) + 1]] <- final_result
            }
        }

        results_df <- do.call(rbind, results)
    }

    # ==================== RICLPM DGP ====================
    else if (dgp == "riclpm") {
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
                    estimator = "allison_chamberlain_fi",
                    dgp = "riclpm"
                )

                model_result <- tryCatch(
                    {
                        warnings_collected <- c()

                        withCallingHandlers(
                            {
                                dat <- simRICLPM(
                                    waves = waves,
                                    beta_x = params$stability_p,
                                    beta_y = params$stability_q,
                                    omega_yx = params$cross_p,
                                    omega_xy = params$cross_q,
                                    var_p = params$variance_p,
                                    var_q = params$variance_q,
                                    cov_pq = cov_pq,
                                    var_BX = params$variance_between_x,
                                    var_BY = params$variance_between_y,
                                    sample.nobs = sample_size
                                )$data

                                model <- lavaan::sem(model_syntax, data = dat)

                                if (!lavInspect(model, "converged")) {
                                    stop("Model did not converge")
                                }

                                param_table <- lavaan::parameterEstimates(model)
                                beta_x <- param_table[param_table$label == "ar_x", "est"]
                                beta_y <- param_table[param_table$label == "ar_y", "est"]
                                omega_xy <- param_table[param_table$label == "cl_xy", "est"]
                                omega_yx <- param_table[param_table$label == "cl_yx", "est"]
                            },
                            warning = function(w) {
                                warnings_collected <<- c(warnings_collected, w$message)
                                invokeRestart("muffleWarning")
                            }
                        )

                        list(
                            success = TRUE,
                            xlag_x = as.numeric(beta_x[1]),
                            ylag_x = as.numeric(omega_yx[1]),
                            xlag_y = as.numeric(omega_xy[1]),
                            ylag_y = as.numeric(beta_y[1]),
                            converged = TRUE,
                            error_message = NA_character_,
                            error_type = NA_character_,
                            warning_messages = if (length(warnings_collected) > 0) paste(warnings_collected, collapse = "; ") else NA_character_,
                            warning_count = length(warnings_collected)
                        )
                    },
                    error = function(e) {
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
                            converged = FALSE,
                            error_message = error_msg,
                            error_type = error_type,
                            warning_messages = NA_character_,
                            warning_count = 0
                        )
                    }
                )

                final_result <- cbind(base_result, model_result[c(
                    "xlag_x", "ylag_x", "xlag_y", "ylag_y",
                    "converged", "error_message", "error_type", "warning_messages", "warning_count"
                )])

                results[[length(results) + 1]] <- final_result
            }
        }

        results_df <- do.call(rbind, results)
    }

    # ==================== CLPMU DGP ====================
    else if (dgp == "clpmu") {
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
        } else {
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
            ) %>% as.data.frame()
        }

        results <- list()

        for (i in 1:nrow(simulation_parameters)) {
            for (j in 1:trials) {
                params <- simulation_parameters[i, ]

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
                    estimator = "allison_chamberlain_fi",
                    dgp = "clpmu"
                )

                if (confounder_type == "time_variant") {
                    base_result$confounder_stability <- params$confounder_stability
                } else {
                    base_result$confounder_stability <- NA
                }

                model_result <- tryCatch(
                    {
                        warnings_collected <- c()

                        withCallingHandlers(
                            {
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
                                        include_confounder = params$include_confounder,
                                        confounder_p = params$confounder_p,
                                        confounder_q = params$confounder_q,
                                        confounder_variance = params$confounder_variance,
                                        confounder_stability = params$confounder_stability,
                                        sample.nobs = sample_size
                                    )$data
                                } else {
                                    dat <- simCLPM_timeInvariantU(
                                        waves = waves,
                                        stability_p = params$stability_p,
                                        stability_q = params$stability_q,
                                        cross_p = params$cross_p,
                                        cross_q = params$cross_q,
                                        variance_p = params$variance_p,
                                        variance_q = params$variance_q,
                                        cov_pq = cov_pq,
                                        include_confounder = params$include_confounder,
                                        confounder_p = params$confounder_p,
                                        confounder_q = params$confounder_q,
                                        confounder_variance = params$confounder_variance,
                                        sample.nobs = sample_size
                                    )$data
                                }

                                model <- lavaan::sem(model_syntax, data = dat)

                                if (!lavInspect(model, "converged")) {
                                    stop("Model did not converge")
                                }

                                param_table <- lavaan::parameterEstimates(model)
                                beta_x <- param_table[param_table$label == "ar_x", "est"]
                                beta_y <- param_table[param_table$label == "ar_y", "est"]
                                omega_xy <- param_table[param_table$label == "cl_xy", "est"]
                                omega_yx <- param_table[param_table$label == "cl_yx", "est"]
                            },
                            warning = function(w) {
                                warnings_collected <<- c(warnings_collected, w$message)
                                invokeRestart("muffleWarning")
                            }
                        )

                        list(
                            success = TRUE,
                            xlag_x = as.numeric(beta_x[1]),
                            ylag_x = as.numeric(omega_yx[1]),
                            xlag_y = as.numeric(omega_xy[1]),
                            ylag_y = as.numeric(beta_y[1]),
                            converged = TRUE,
                            error_message = NA_character_,
                            error_type = NA_character_,
                            warning_messages = if (length(warnings_collected) > 0) paste(warnings_collected, collapse = "; ") else NA_character_,
                            warning_count = length(warnings_collected)
                        )
                    },
                    error = function(e) {
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
                            converged = FALSE,
                            error_message = error_msg,
                            error_type = error_type,
                            warning_messages = NA_character_,
                            warning_count = 0
                        )
                    }
                )

                final_result <- cbind(base_result, model_result[c(
                    "xlag_x", "ylag_x", "xlag_y", "ylag_y",
                    "converged", "error_message", "error_type", "warning_messages", "warning_count"
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
