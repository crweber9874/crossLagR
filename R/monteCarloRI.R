#' @title monteCarloRI
#' @description Monte Carlo simulation for Random Intercepts estimation using lmer
#' Allows the user to specify different Data Generating conditions, and apply the Random Intercepts model to the data.
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
#' @param verbose Whether to print progress and warnings
#' @param control_params List of control parameters to pass to lmer(). Default uses
#'   lmerControl with optimizer = "bobyqa" and optCtrl = list(maxfun = 20000).
#' @param ... Additional arguments
#'
#' @return A data frame containing the results of the simulation.
#'
#' @details
#' Random Intercepts models use lme4::lmer() to fit mixed-effects models with:
#' \itemize{
#'   \item Random intercepts for each individual
#'   \item Fixed effects for autoregressive and cross-lagged parameters
#'   \item Assumptions that individual intercepts are drawn from a normal distribution
#' }
#'
#' This approach is intermediate between OLS (which pools individuals) and
#' Fixed Individual Effects (which estimates separate intercepts for each individual).
#' Random Intercepts models can handle unbalanced data and are more efficient when
#' the random intercept assumption is appropriate.
#'
#' @examples
#' \dontrun{
#' # Basic Random Intercepts Monte Carlo
#' results <- monteCarloRI(
#'   trials = 100,
#'   waves = 4,
#'   dgp = "riclpm",
#'   sample_size = 1000
#' )
#'
#' # Test with confounded data
#' results_confounded <- monteCarloRI(
#'   trials = 50,
#'   waves = 4,
#'   dgp = "clpmu",
#'   confounder_type = "time_variant",
#'   confounder_p = 0.3,
#'   confounder_q = 0.3,
#'   sample_size = 1000
#' )
#'
#' # Compare with time-invariant confounder
#' results_time_invariant <- monteCarloRI(
#'   trials = 50,
#'   waves = 4,
#'   dgp = "clpmu",
#'   confounder_type = "time_invariant",
#'   sample_size = 1000
#' )
#' }
#'
#' @import dplyr lme4 stats
#' @export
monteCarloRI <- function(
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
    confounder_p = 0.3,
    confounder_q = 0.3,
    confounder_variance = 1,
    confounder_stability = 0.4,
    include_confounder = TRUE,
    confounder_type = "time_variant",
    verbose = FALSE,
    control_params = lme4::lmerControl(optimizer = "bobyqa",
                                                        optCtrl = list(maxfun = 20000)),

    ...
) {
  valid_dgps <- c("riclpm", "clpm", "clpmu")
  if (!dgp %in% valid_dgps) {
    stop("The data-generation should be random intercept cross lagged panel regression (riclpm),
         the cross lagged regression (clpm), or the cross lagged panel regression model with confounder
         (clpmu): ", paste(valid_dgps, collapse = ", "))
  }

  if (!confounder_type %in% c("time_variant", "time_invariant")) {
    stop("confounder_type 'time_variant' or 'time_invariant'")
  }

  ## Build appropriate simulation paradigm based on specificed dgp
  ## The reason for this is that the types of parameters that may vary
  ## differ based on DGP. for instance, clpm has no confounders,
  ## But clpmu does
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
      )
      simulation_parameters <- as.data.frame(simulation_parameters)
    }
  }

  # Initialize results list
  results <- list()

  for (i in 1:nrow(simulation_parameters)) {
    for (j in 1:trials) {
      params <- simulation_parameters[i, ]

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
        } else {
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
      # Catch errors during estimation
      # lavaan produces a lot of these for the riclpm
      # its not clear why?
      tryCatch({
        long_data <- reshape_long_sim_cr(raw_data)
        long_data <- as.data.frame(long_data)
        long_data <- na.omit(long_data)

        ri_results <- estimateRI(
          data = long_data,
          control_params = control_params
        )
      # output results -- these are all constant across every data generation type
        result_row <- data.frame(
          xlag_x = ri_results$xlag_x,
          ylag_x = ri_results$ylag_x,
          xlag_y = ri_results$xlag_y,
          ylag_y = ri_results$ylag_y,
          converged = ri_results$converged,
          n_obs = ri_results$n_obs,
          n_groups = ri_results$n_groups,
          trial = j,
          param_combo = i,
          estimator = "ri_lmer",
          dgp = dgp,
          stringsAsFactors = FALSE
        )

        if (!is.null(ri_results$random_effects)) {
          result_row$random_int_var_x <- ri_results$random_effects$x_model$random_intercept_var
          result_row$residual_var_x <- ri_results$random_effects$x_model$residual_var
          result_row$random_int_var_y <- ri_results$random_effects$y_model$random_intercept_var
          result_row$residual_var_y <- ri_results$random_effects$y_model$residual_var
        }

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

          if (confounder_type == "time_variant") {
            result_row$confounder_stability <- params$confounder_stability
          } else {
            result_row$confounder_stability <- NA
          }
        }

        results[[length(results) + 1]] <- result_row

      }, error = function(e) {
        if (verbose) {
          cat("Error in trial", j, "param combo", i, ":", e$message, "\n")
        }

        # Create error result with parameter info
        result_row <- data.frame(
          xlag_x = NA_real_,
          ylag_x = NA_real_,
          xlag_y = NA_real_,
          ylag_y = NA_real_,
          converged = FALSE,
          n_obs = NA_integer_,
          n_groups = NA_integer_,
          trial = j,
          param_combo = i,
          estimator = "ri_lmer",
          dgp = dgp,
          error_message = as.character(e$message),
          stringsAsFactors = FALSE
        )

        # Add dgp-specific parameters for error case too
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

          if (confounder_type == "time_variant") {
            result_row$confounder_stability <- params$confounder_stability
          } else {
            result_row$confounder_stability <- NA
          }
        }

        results[[length(results) + 1]] <- result_row
      })
    }
  }

  # Combine all results
  results_df <- do.call(rbind, results)
  row.names(results_df) <- NULL

  return(results_df)
}
