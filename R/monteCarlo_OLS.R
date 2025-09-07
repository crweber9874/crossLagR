monteCarloOLS <- function(
    trials = 10,
    waves = 5,
    variance_between_x = 0.5,
    variance_between_y = 0.5,
    stability_q = 0.25,
    stability_p = 0.25,
    cross_p = 0,
    cross_q = 0,
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
    ...
) {
  if(dgp == "riclpm") {
    # Existing riclpm code...

  } else if(dgp == "clpm") {
    # Existing clpm code...

  } else if(dgp == "clpmu") {  # ADD THIS NEW SECTION
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
    ) %>% as.data.frame()

    results <- list()

    for (i in 1:nrow(simulation_parameters)) {
      for (j in 1:trials) {
        params <- simulation_parameters[i, ]

        # Use simCLPMu instead of simRICLPM
        dat <- simCLPMu(
          waves = waves,
          stability_p = params$stability_p,
          stability_q = params$stability_q,
          cross_p = params$cross_p,
          cross_q = params$cross_q,
          variance_p = params$variance_p,
          variance_q = params$variance_q,
          cov_pq = 0.1,  # Add this parameter
          include_confounder = params$include_confounder,
          confounder_p = params$confounder_p,
          confounder_q = params$confounder_q,
          confounder_variance = params$confounder_variance,
          confounder_stability = params$confounder_stability,
          sample.nobs = sample_size
        )$data %>%
          reshape_long_sim_cr() %>% as.data.frame() %>% na.omit()

        # Same OLS fitting as before
        model_fit_y <- lm(y ~ xlag + ylag, dat)
        model_fit_x <- lm(x ~ xlag + ylag, dat)

        xvalue_x <- model_fit_x$coefficients[["xlag"]]
        yvalue_x <- model_fit_x$coefficients[["ylag"]]
        xvalue_y <- model_fit_y$coefficients[["xlag"]]
        yvalue_y <- model_fit_y$coefficients[["ylag"]]

        # Store results
        results[[length(results) + 1]] <- data.frame(
          variance_between_x = params$variance_between_x,
          variance_between_y = params$variance_between_y,
          stability_p = params$stability_p,
          stability_q = params$stability_q,
          cross_p = params$cross_p,
          cross_q = params$cross_q,
          variance_p = params$variance_p,
          variance_q = params$variance_q,
          confounder_p = params$confounder_p,
          confounder_q = params$confounder_q,
          include_confounder = params$include_confounder,
          xlag_x = xvalue_x,
          ylag_x = yvalue_x,
          xlag_y = xvalue_y,
          ylag_y = yvalue_y,
          trial = j,
          estimator = "ols",
          dgp = "clpmu"
        )
      }
    }
    results_df <- do.call(rbind, results)  # This was missing!
  }

  return(results_df)  # Make sure this is outside all if/else blocks
}
