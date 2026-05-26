#' @title monteCarloTraitState
#' @description Monte Carlo simulation that varies the trait-to-state variance ratio
#'   and fits four competing panel models: Bollen & Brand, RI-CLPM, CLPM, and
#'   Latent Change Score Model.
#'
#' Data are generated from a trait-plus-state DGP via \code{sim_trait_state()}.
#' The key manipulation is the ratio of between-person (trait) variance to
#' within-person (state) variance. When trait variance is 0 the DGP is a pure CLPM;
#' as trait variance increases the data increasingly violate CLPM assumptions.
#'
#' @param trials Integer. Number of Monte Carlo replications per condition. Default is 100.
#' @param waves Integer. Number of waves. Must be >= 5 for latent change model. Default is 5.
#' @param n Integer. Sample size per replication. Default is 1000.
#' @param beta_x Numeric. Autoregressive effect for X within-person process. Default is 0.3.
#' @param beta_y Numeric. Autoregressive effect for Y within-person process. Default is 0.3.
#' @param omega_xy Numeric. Cross-lagged effect from X to Y (within-person). Default is 0.1.
#' @param omega_yx Numeric. Cross-lagged effect from Y to X (within-person). Default is 0.1.
#' @param var_p Numeric. Within-person innovation variance for X. Default is 1.
#' @param var_q Numeric. Within-person innovation variance for Y. Default is 1.
#' @param cov_pq Numeric. Within-person innovation covariance. Default is 0.1.
#' @param trait_ratio Numeric vector. Ratio(s) of trait variance to state variance.
#'   For each value r, var_BX = var_BY = r * var_p. Default is c(0, 0.5, 1, 2).
#' @param cov_BXBY Numeric or NULL. Trait covariance. If NULL, set to 0.3 * var_BX. Default is NULL.
#' @param verbose Logical. Print progress. Default is TRUE.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame with one row per trial x condition x estimator, containing:
#'   \describe{
#'     \item{trait_ratio}{The trait-to-state variance ratio used for simulation}
#'     \item{var_BX, var_BY}{The actual trait variances used}
#'     \item{estimator}{Which model was fit: "clpm", "riclpm", "bollen_brand", "lchange"}
#'     \item{beta_x_est, beta_y_est}{Estimated autoregressive / proportional effects}
#'     \item{omega_xy_est, omega_yx_est}{Estimated cross-lagged / coupling effects}
#'     \item{converged}{Logical. Did the model converge?}
#'     \item{trial}{Trial number}
#'   }
#'
#' @details
#' For a single sample (trials = 1), this can be used as a one-shot demonstration.
#' For Monte Carlo (trials > 1), bias and RMSE can be computed by comparing estimates
#' to the true DGP values.
#'
#' The four estimators extract different types of parameters:
#' \itemize{
#'   \item \strong{CLPM}: autoregressive (beta) and cross-lagged (omega) on observed composites
#'   \item \strong{RI-CLPM}: autoregressive and cross-lagged on within-person latent factors
#'   \item \strong{Bollen & Brand}: autoregressive (rho) and cross-lagged (b1, b2) with latent fixed effects
#'   \item \strong{Latent Change}: proportional (beta) and coupling (omega) effects on change scores
#' }
#'
#' @examples
#' \dontrun{
#' # Single-sample demonstration
#' demo <- monteCarloTraitState(trials = 1, trait_ratio = c(0, 1, 2), n = 500)
#' demo[, c("trait_ratio", "estimator", "beta_x_est", "omega_xy_est")]
#'
#' # Full Monte Carlo
#' mc <- monteCarloTraitState(trials = 100, trait_ratio = c(0, 0.5, 1, 2), n = 1000)
#'
#' # Compute bias
#' library(dplyr)
#' mc %>%
#'   filter(converged) %>%
#'   group_by(trait_ratio, estimator) %>%
#'   summarise(
#'     bias_omega_xy = mean(omega_xy_est - 0.1, na.rm = TRUE),
#'     bias_omega_yx = mean(omega_yx_est - 0.1, na.rm = TRUE),
#'     rmse_omega_xy = sqrt(mean((omega_xy_est - 0.1)^2, na.rm = TRUE))
#'   )
#' }
#'
#' @import lavaan
#' @importFrom dplyr case_when
#' @export
monteCarloTraitState <- function(
    trials = 100,
    waves = 5,
    n = 1000,
    beta_x = 0.3,
    beta_y = 0.3,
    omega_xy = 0.1,
    omega_yx = 0.1,
    var_p = 1,
    var_q = 1,
    cov_pq = 0.1,
    trait_ratio = c(0, 0.5, 1, 2),
    cov_BXBY = NULL,
    verbose = TRUE,
    ...
) {
    library(lavaan)

    # Input validation
    if (waves < 5) {
        stop("Latent change score model requires at least 5 waves.")
    }
    if (n < 50) {
        stop("Sample size must be at least 50.")
    }
    if (any(trait_ratio < 0)) {
        stop("trait_ratio values must be non-negative.")
    }

    # Pre-build model syntaxes (these don't change across trials)
    clpm_syntax <- estimateCLPM(waves = waves)
    riclpm_syntax <- estimateRICLPM(waves = waves)
    bollen_syntax <- estimateBollen_and_Brand(
        waves = waves,
        x_effect = "lagged",
        x_autoregression = TRUE,
        y_effect_on_x = TRUE
    )
    lchange_syntax <- estimateLChange(
        waves = waves,
        variable_type = "bivariate",
        estimate_constant_change = FALSE,
        estimate_change_to_change = FALSE
    )

    # Helper: safely fit a lavaan model and extract parameters
    fit_and_extract <- function(syntax, dat, estimator_name) {
        tryCatch({
            warnings_collected <- c()

            result <- withCallingHandlers(
                {
                    fit <- lavaan::sem(syntax, data = dat, meanstructure = TRUE)

                    if (!lavInspect(fit, "converged")) {
                        stop("Model did not converge")
                    }

                    pt <- lavaan::parameterEstimates(fit)

                    # Extract parameters based on estimator
                    if (estimator_name == "clpm") {
                        list(
                            beta_x_est = pt$est[pt$label == "ar_x"][1],
                            beta_y_est = pt$est[pt$label == "ar_y"][1],
                            omega_xy_est = pt$est[pt$label == "cl_xy"][1],
                            omega_yx_est = pt$est[pt$label == "cl_yx"][1]
                        )
                    } else if (estimator_name == "riclpm") {
                        list(
                            beta_x_est = pt$est[pt$label == "ar_x"][1],
                            beta_y_est = pt$est[pt$label == "ar_y"][1],
                            omega_xy_est = pt$est[pt$label == "cl_xy"][1],
                            omega_yx_est = pt$est[pt$label == "cl_yx"][1]
                        )
                    } else if (estimator_name == "bollen_brand") {
                        list(
                            beta_x_est = pt$est[pt$label == "ar_x"][1],
                            beta_y_est = pt$est[pt$label == "ar_y"][1],
                            omega_xy_est = pt$est[pt$label == "cl_xy"][1],
                            omega_yx_est = pt$est[pt$label == "cl_yx"][1]
                        )
                    } else if (estimator_name == "lchange") {
                        list(
                            beta_x_est = pt$est[pt$label == "ar_x"][1],
                            beta_y_est = pt$est[pt$label == "ar_y"][1],
                            omega_xy_est = pt$est[pt$label == "cl_xy"][1],
                            omega_yx_est = pt$est[pt$label == "cl_yx"][1]
                        )
                    }
                },
                warning = function(w) {
                    warnings_collected <<- c(warnings_collected, w$message)
                    invokeRestart("muffleWarning")
                }
            )

            list(
                beta_x_est = as.numeric(result$beta_x_est),
                beta_y_est = as.numeric(result$beta_y_est),
                omega_xy_est = as.numeric(result$omega_xy_est),
                omega_yx_est = as.numeric(result$omega_yx_est),
                converged = TRUE,
                error_message = NA_character_,
                warning_count = length(warnings_collected)
            )
        },
        error = function(e) {
            list(
                beta_x_est = NA_real_,
                beta_y_est = NA_real_,
                omega_xy_est = NA_real_,
                omega_yx_est = NA_real_,
                converged = FALSE,
                error_message = as.character(e$message),
                warning_count = 0
            )
        })
    }

    # Main simulation loop
    results <- list()
    total_conditions <- length(trait_ratio)
    total_iterations <- total_conditions * trials

    for (r_idx in seq_along(trait_ratio)) {
        r <- trait_ratio[r_idx]
        var_BX <- r * var_p
        var_BY <- r * var_q
        trait_cov <- if (is.null(cov_BXBY)) 0.3 * sqrt(var_BX * var_BY) else cov_BXBY

        # Ensure covariance is valid (|cov| <= sqrt(var1 * var2))
        if (var_BX > 0 && var_BY > 0) {
            max_cov <- sqrt(var_BX * var_BY)
            trait_cov <- max(min(trait_cov, max_cov * 0.95), -max_cov * 0.95)
        } else {
            trait_cov <- 0
        }

        if (verbose) {
            cat(sprintf("Condition %d/%d: trait_ratio = %.2f (var_BX = %.2f, var_BY = %.2f)\n",
                        r_idx, total_conditions, r, var_BX, var_BY))
        }

        for (j in 1:trials) {
            if (verbose && j %% 10 == 0) {
                cat(sprintf("  Trial %d/%d\n", j, trials))
            }

            # Simulate data
            dat <- sim_trait_state(
                waves = waves,
                n = n,
                beta_x = beta_x,
                beta_y = beta_y,
                omega_xy = omega_xy,
                omega_yx = omega_yx,
                var_p = var_p,
                var_q = var_q,
                cov_pq = cov_pq,
                var_BX = var_BX,
                var_BY = var_BY,
                cov_BXBY = trait_cov
            )

            # Fit all four models to the same dataset
            estimators <- c("clpm", "riclpm", "bollen_brand", "lchange")
            syntaxes <- list(clpm_syntax, riclpm_syntax, bollen_syntax, lchange_syntax)

            for (e_idx in seq_along(estimators)) {
                est_name <- estimators[e_idx]
                est_result <- fit_and_extract(syntaxes[[e_idx]], dat, est_name)

                row <- data.frame(
                    trait_ratio = r,
                    var_BX = var_BX,
                    var_BY = var_BY,
                    true_beta_x = beta_x,
                    true_beta_y = beta_y,
                    true_omega_xy = omega_xy,
                    true_omega_yx = omega_yx,
                    estimator = est_name,
                    beta_x_est = est_result$beta_x_est,
                    beta_y_est = est_result$beta_y_est,
                    omega_xy_est = est_result$omega_xy_est,
                    omega_yx_est = est_result$omega_yx_est,
                    converged = est_result$converged,
                    error_message = est_result$error_message,
                    warning_count = est_result$warning_count,
                    trial = j,
                    waves = waves,
                    n = n,
                    stringsAsFactors = FALSE
                )

                results[[length(results) + 1]] <- row
            }
        }
    }

    results_df <- do.call(rbind, results)
    rownames(results_df) <- NULL

    if (verbose) {
        cat("\n=== Summary ===\n")
        cat(sprintf("Total trials: %d conditions x %d trials = %d datasets\n",
                    total_conditions, trials, total_conditions * trials))
        cat(sprintf("Models fit per dataset: 4 (CLPM, RI-CLPM, Bollen & Brand, LChange)\n"))
        cat(sprintf("Total model fits: %d\n", nrow(results_df)))
        conv_rate <- tapply(results_df$converged, results_df$estimator, mean, na.rm = TRUE)
        cat("Convergence rates by estimator:\n")
        print(round(conv_rate, 3))
    }

    return(results_df)
}
