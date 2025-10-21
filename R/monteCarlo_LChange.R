#' @title monteCarloLChange
#' @description Monte Carlo simulation for Latent Change Score Models
#' This function combines the estimation and simulation of latent change models over a fixed number of trials.
#' The output is a data frame that includes the estimated coefficients across these trials
#' @param trials The number of trials for the Monte Carlo simulation.
#' @param waves The number of waves (time points) in the model. Must specify three or more waves for identification.
#' @param dgp Specify the data generation method: "riclpm", "clpm", "clpmu", "lchange". Are data generated under one of these regimes?
#' @param model_type Type of latent change model: "latent_change" or "dual_change". Use dual change with multiple variables
#' @param stability_p Autoregressive effect for p (x) variable. These parameters are only relevant when specifying an RICLPM or CLPM. They are **not** part of the dual change model
#' @param stability_q Autoregressive effect for q (y) variable. These parameters are only relevant when specifying an RICLPM or CLPM. They are **not** part of the dual change model
#' @param cross_p Cross-lagged effect of y on x. Only relevant to the CLPM and RICLPM specification.
#' @param cross_q Cross-lagged effect of x on y. Only relevant to the CLPM and RICLPM specification.
#' @param variance_p Within-person variance for p variable. Only relevant to the RICLPM and CLPM specification.
#' @param variance_q Within-person variance for q variable. Only relevant ot the RICLPM and CLPM specification
#' @param variance_between_x Between-person variance for x variable. Only relevant for RICLPM.
#' @param variance_between_y Between-person variance for y variable. Ony relevant for the RICLPM.
#' @param sample_size Sample size for simulation
#' @param confounder_p Effect of confounder on x variables (for clpmu)
#' @param confounder_q Effect of confounder on y variables (for clpmu)
#' @param confounder_variance Variance of the confounder (for clpmu)
#' @param confounder_stability Autoregressive effect of confounder (for clpmu)
#' @param include_confounder Whether to include confounder in clpmu model
#' @param confounder_type Character string specifying confounder type for clpmu dgp.
#'   Must be one of "time_variant" (default) or "time_invariant".
#' @param verbose Whether to print progress and error messages.
#' @param cov_pq Covariance between p and q (for CLPM/CLPMU)
#' @param ar_x Proportional effect for X in latent change DGP
#' @param ar_y Proportional effect for Y in latent change DGP
#' @param cl_x Coupling effect in latent change DGP
#' @param cl_y Coupling effect in latent change DGP
#' @param change_x Change-to-change effect in latent change DGP
#' @param change_y Change-to-change effect in latent change DGP
#' @param phi_x Change autoregression in latent change DGP
#' @param phi_y Change autoregression in latent change DGP
#' @param initial_var_x Initial variance for X in latent change DGP
#' @param initial_var_y Initial variance for Y in latent change DGP
#' @param constant_change_var_x Constant change variance for X in latent change DGP
#' @param constant_change_var_y Constant change variance for Y in latent change DGP
#' @param residual_variance_x Measurement error for X in latent change DGP
#' @param residual_variance_y Measurement error for Y in latent change DGP
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
    dgp = "lchange",
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
    cov_pq = 0.1,
    confounder_p = 0.3,
    confounder_q = 0.3,
    confounder_variance = 1,
    confounder_stability = 0.4,
    include_confounder = TRUE,
    confounder_type = "time_variant",
    ar_x = -0.15,
    ar_y = -0.20,
    cl_x = -0.10,
    cl_y = -0.10,
    change_x = 0.08,
    change_y = 0.08,
    phi_x = 0.15,
    phi_y = 0.15,
    initial_var_x = 1,
    initial_var_y = 1,
    constant_change_var_x = 0.5,
    constant_change_var_y = 0.5,
    residual_variance_x = 0.5,
    residual_variance_y = 0.5,
    verbose = FALSE,
    ...
) {
  library(lavaan)

  # DGP input, validation
  valid_dgps <- c("riclpm", "clpm", "clpmu", "lchange") # RICLPM, CLPM, CLMP with confounder, latent change
  valid_model_types <- c("latent_change", "dual_change")

  if (!dgp %in% valid_dgps) {
    stop("dgp must be one of: ", paste(valid_dgps, collapse = ", "))
  }

  if (!model_type %in% valid_model_types) {
    stop("model_type must be one of: ", paste(valid_model_types, collapse = ", "))
  }

  if (!confounder_type %in% c("time_variant", "time_invariant")) {
    stop("confounder_type must be either 'time_variant' or 'time_invariant'")
  }

  get_coef_safe <- function(coef_name, coeffs) {
    if (coef_name %in% names(coeffs)) {
      return(as.numeric(coeffs[coef_name]))
    } else {
      return(NA_real_)
    }
  }
# RICLPM
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

    results <- list()

    for (i in 1:nrow(simulation_parameters)) {
      for (j in 1:trials) {
        params <- simulation_parameters[i, ]

        if (verbose && j %% max(1, floor(trials/10)) == 0) {
          cat("DGP:", dgp, "| Parameter combination", i, "| Trial", j, "of", trials, "\n")
        }

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

        model_result <- tryCatch({
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

          model_syntax <- estimateLChange(
            waves = waves,
            model_type = model_type
          )

          model <- lavaan::lavaan(
            model_syntax,
            data = dat,
            optim.method = "nlminb",
            control = list(
              iter.max = 10000,
              eval.max = 10000
            )
          )

          if (!lavInspect(model, "converged")) {
            stop("Model did not converge")
          }

          coeffs <- coef(model)

          if (model_type == "dual_change") {
            change_x_to_y <- get_coef_safe("change_y", coeffs)
            change_y_to_x <- get_coef_safe("change_x", coeffs)
            coupling_x_to_y <- get_coef_safe("cl_x", coeffs)
            coupling_y_to_x <- get_coef_safe("cl_y", coeffs)
            ar_x <- get_coef_safe("ar_x", coeffs)
            ar_y <- get_coef_safe("ar_y", coeffs)
            phi_x <- get_coef_safe("phi_x", coeffs)
            phi_y <- get_coef_safe("phi_y", coeffs)

            xlag_x <- ar_x
            ylag_y <- ar_y
            xlag_y <- change_x_to_y
            ylag_x <- change_y_to_x
          } else {
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

  else if (dgp == "clpm") {
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
          estimator = "lchange",
          model_type = model_type,
          dgp = "clpm"
        )

        model_result <- tryCatch({
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

          model_syntax <- estimateLChange(
            waves = waves,
            model_type = model_type
          )

          model <- lavaan::lavaan(
            model_syntax,
            data = dat,
            optim.method = "nlminb",
            control = list(iter.max = 10000, eval.max = 10000)
          )

          if (!lavInspect(model, "converged")) {
            stop("Model did not converge")
          }

          coeffs <- coef(model)

          if (model_type == "dual_change") {
            change_x_to_y <- get_coef_safe("change_y", coeffs)
            change_y_to_x <- get_coef_safe("change_x", coeffs)
            coupling_x_to_y <- get_coef_safe("cl_x", coeffs)
            coupling_y_to_x <- get_coef_safe("cl_y", coeffs)
            ar_x <- get_coef_safe("ar_x", coeffs)
            ar_y <- get_coef_safe("ar_y", coeffs)
            phi_x <- get_coef_safe("phi_x", coeffs)
            phi_y <- get_coef_safe("phi_y", coeffs)

            xlag_x <- ar_x
            ylag_y <- ar_y
            xlag_y <- change_x_to_y
            ylag_x <- change_y_to_x
          } else {
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
            cat("Error in trial", j, "parameter combination", i, ":", error_msg, "\n")
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
# confounder specification
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
          estimator = "lchange",
          model_type = model_type,
          dgp = "clpmu"
        )

        if (confounder_type == "time_variant") {
          base_result$confounder_stability <- params$confounder_stability
        } else {
          base_result$confounder_stability <- NA
        }

        model_result <- tryCatch({
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
              confounder_p = params$confounder_p,
              confounder_q = params$confounder_q,
              confounder_variance = params$confounder_variance,
              sample.nobs = sample_size
            )$data
          }

          model_syntax <- estimateLChange(
            waves = waves,
            model_type = model_type
          )

          model <- lavaan::lavaan(
            model_syntax,
            data = dat,
            optim.method = "nlminb",
            control = list(iter.max = 10000, eval.max = 10000)
          )

          if (!lavInspect(model, "converged")) {
            stop("Model did not converge")
          }

          coeffs <- coef(model)

          if (model_type == "dual_change") {
            change_x_to_y <- get_coef_safe("change_y", coeffs)
            change_y_to_x <- get_coef_safe("change_x", coeffs)
            coupling_x_to_y <- get_coef_safe("cl_x", coeffs)
            coupling_y_to_x <- get_coef_safe("cl_y", coeffs)
            ar_x <- get_coef_safe("ar_x", coeffs)
            ar_y <- get_coef_safe("ar_y", coeffs)
            phi_x <- get_coef_safe("phi_x", coeffs)
            phi_y <- get_coef_safe("phi_y", coeffs)

            xlag_x <- ar_x
            ylag_y <- ar_y
            xlag_y <- change_x_to_y
            ylag_x <- change_y_to_x
          } else {
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
            cat("Error in trial", j, "parameter combination", i, ":", error_msg, "\n")
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

  else if (dgp == "lchange") {
    # NEW: Generate from latent change model
    simulation_parameters <- expand.grid(
      ar_x = ar_x,
      ar_y = ar_y,
      cl_x = cl_x,
      cl_y = cl_y,
      change_x = change_x,
      change_y = change_y,
      phi_x = phi_x,
      phi_y = phi_y,
      initial_var_x = initial_var_x,
      initial_var_y = initial_var_y,
      constant_change_var_x = constant_change_var_x,
      constant_change_var_y = constant_change_var_y,
      residual_variance_x = residual_variance_x,
      residual_variance_y = residual_variance_y
    ) %>% as.data.frame()

    results <- list()

    for (i in 1:nrow(simulation_parameters)) {
      for (j in 1:trials) {
        params <- simulation_parameters[i, ]

        base_result <- data.frame(
          ar_x_true = params$ar_x,
          ar_y_true = params$ar_y,
          cl_x_true = params$cl_x,
          cl_y_true = params$cl_y,
          change_x_true = params$change_x,
          change_y_true = params$change_y,
          phi_x_true = params$phi_x,
          phi_y_true = params$phi_y,
          initial_var_x = params$initial_var_x,
          initial_var_y = params$initial_var_y,
          constant_change_var_x = params$constant_change_var_x,
          constant_change_var_y = params$constant_change_var_y,
          residual_variance_x = params$residual_variance_x,
          residual_variance_y = params$residual_variance_y,
          trial = j,
          param_combo = i,
          estimator = "lchange",
          model_type = model_type,
          dgp = "lchange"
        )

        model_result <- tryCatch({
          # Generate from latent change model
          dat <- simLChange(
            waves = waves,
            model_type = model_type,
            ar_x = params$ar_x,
            ar_y = params$ar_y,
            cl_x = params$cl_x,
            cl_y = params$cl_y,
            change_x = params$change_x,
            change_y = params$change_y,
            phi_x = params$phi_x,
            phi_y = params$phi_y,
            initial_var_x = params$initial_var_x,
            initial_var_y = params$initial_var_y,
            constant_change_var_x = params$constant_change_var_x,
            constant_change_var_y = params$constant_change_var_y,
            residual_variance_x = params$residual_variance_x,
            residual_variance_y = params$residual_variance_y,
            sample.nobs = sample_size
          )$data

          model_syntax <- estimateLChange(
            waves = waves,
            model_type = model_type
          )

          model <- lavaan::lavaan(
            model_syntax,
            data = dat,
            optim.method = "nlminb",
            control = list(iter.max = 10000, eval.max = 10000)
          )

          if (!lavInspect(model, "converged")) {
            stop("Model did not converge")
          }

          coeffs <- coef(model)

          if (model_type == "dual_change") {
            change_x_to_y <- get_coef_safe("change_y", coeffs)
            change_y_to_x <- get_coef_safe("change_x", coeffs)
            coupling_x_to_y <- get_coef_safe("cl_x", coeffs)
            coupling_y_to_x <- get_coef_safe("cl_y", coeffs)
            ar_x <- get_coef_safe("ar_x", coeffs)
            ar_y <- get_coef_safe("ar_y", coeffs)
            phi_x <- get_coef_safe("phi_x", coeffs)
            phi_y <- get_coef_safe("phi_y", coeffs)

            xlag_x <- ar_x
            ylag_y <- ar_y
            xlag_y <- change_x_to_y
            ylag_x <- change_y_to_x
          } else {
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
            cat("Error in trial", j, "parameter combination", i, ":", error_msg, "\n")
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

        final_result <- cbind(base_result, model_result[c(
          "xlag_x", "ylag_x",
          "xlag_y", "ylag_y",
          "change_x_to_y", "change_y_to_x",
          "coupling_x_to_y", "coupling_y_to_x",
          "ar_x", "ar_y",
          "phi_x", "phi_y",
          "converged", "error_message",
          "error_type"
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
