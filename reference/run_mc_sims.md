# run_mc_sims

Run Monte Carlo simulations across different estimators and parameter
combinations

This function provides a unified interface for running Monte Carlo
simulations across multiple estimators (OLS, RICLPM, CLPM, CTSEM, FI,
LCHANGE, RI) with varying parameter combinations. It handles parameter
grid expansion, error catching, and result compilation automatically.

## Usage

``` r
run_mc_sims(
  estimator,
  riclpm_type = "riclpm",
  param_grid = NULL,
  lchange_type = "dual_change",
  trials = 10,
  waves = 3,
  sample_size = 1000,
  verbose = TRUE,
  data_generation = "riclpm"
)
```

## Arguments

- estimator:

  Character string specifying which estimator to use. Must be one of:

  - "OLS" - Ordinary Least Squares estimation

  - "RICLPM" - Random Intercept Cross-Lagged Panel Model

  - "CLPM" - Cross-Lagged Panel Model

  - "CTSEM" - Continuous Time Structural Equation Model

  - "FI" - Fixed Individual effects (First Differences)

  - "LCHANGE" - Latent Change Score Model

  - "RI" - Random Intercepts using lmer

- riclpm_type:

  Character string specifying which RICLPM variant to use (only relevant
  when estimator = "RICLPM"). Must be one of:

  - "riclpm" - Regular RICLPM with autoregressive and cross-lagged
    effects (default)

  - "riclpm_nolag" - RICLPM without autoregressive effects (cross-lagged
    only)

- param_grid:

  Data frame or NULL. If NULL, uses default parameter grid. If provided,
  should contain columns for parameters to vary across simulations.
  Missing parameters will be filled with defaults.

- lchange_type:

  Character string specifying latent change model type (only relevant
  when estimator = "LCHANGE"). Must be one of:

  - "dual_change" - Bivariate latent change model (default)

  - "latent_change" - Single variable latent change model

- trials:

  Integer. Number of simulation trials to run for each parameter
  combination. Default is 10.

- waves:

  Integer. Number of time waves in the simulated data. Default is 3.

- sample_size:

  Integer. Sample size for each simulated dataset. Default is 1000.

- verbose:

  Logical. If TRUE, prints progress information during simulation.
  Default is TRUE.

- data_generation:

  Character string specifying the data generation process. Must be one
  of:

  - "riclpm" - Random Intercept Cross-Lagged Panel Model data generation

  - "clpm" - Cross-Lagged Panel Model data generation

  - "clpmu" - CLPM with unmeasured confounder (uses simCLPMu)

## Value

A data frame containing simulation results with the following columns:

- Parameter values used in simulation

- Estimated coefficients (varies by estimator)

- trial - Trial number

- estimator - Which estimator was used

- dgp - Data generation process used

- error_occurred - TRUE if an error occurred in that trial

- error_message - Error message if error_occurred is TRUE

## Details

The function performs the following steps:

1.  Validates input parameters

2.  Creates or validates the parameter grid

3.  Adds default values for missing parameters

4.  Loops through each parameter combination

5.  Calls the appropriate Monte Carlo function

6.  Handles errors and compiles results

7.  Returns a combined data frame of all results

Default parameter values:

- stability_p = 0.2 (autoregressive effect for X)

- stability_q = 0.5 (autoregressive effect for Y)

- cross_p = 0.0 (cross-lagged effect Y -\> X)

- cross_q = 0.0 (cross-lagged effect X -\> Y)

- variance_q = 0.5 (within-person variance for Y)

- variance_p = 0.5 (within-person variance for X)

- variance_between_x = 0.5 (between-person variance for X)

- variance_between_y = 0.5 (between-person variance for Y)

- cov_pq = 0 (within-time covariance between X and Y)

When data_generation = "clpmu", additional confounder parameters:

- confounder_p = 0.3 (effect of confounder on X)

- confounder_q = 0.3 (effect of confounder on Y)

- confounder_variance = 1 (variance of confounder)

- confounder_stability = 0.4 (autoregressive effect of confounder)

- include_confounder = TRUE (whether to include confounder)

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage with default parameters
results <- run_mc_sims(
  estimator = "RICLPM",
  trials = 5,
  waves = 3,
  sample_size = 500
)

# Test Random Intercepts with both full and cross-lagged only
results_ri_full <- run_mc_sims(
  estimator = "RI",
  param_grid = data.frame(include_lagged_dv = TRUE),
  trials = 10,
  waves = 4,
  sample_size = 1000
)

results_ri_cross <- run_mc_sims(
  estimator = "RI",
  param_grid = data.frame(include_lagged_dv = FALSE),
  trials = 10,
  waves = 4,
  sample_size = 1000
)

# Compare regular RICLPM vs no-lag RICLPM
results_regular <- run_mc_sims(
  estimator = "RICLPM",
  riclpm_type = "riclpm",
  trials = 10,
  waves = 4,
  sample_size = 1000
)

results_nolag <- run_mc_sims(
  estimator = "RICLPM",
  riclpm_type = "riclpm_nolag",
  trials = 10,
  waves = 4,
  sample_size = 1000
)

# Compare latent change models
results_dual_change <- run_mc_sims(
  estimator = "LCHANGE",
  lchange_type = "dual_change",
  trials = 10,
  waves = 4,
  sample_size = 1000
)

results_single_change <- run_mc_sims(
  estimator = "LCHANGE",
  lchange_type = "latent_change",
  trials = 10,
  waves = 4,
  sample_size = 1000
)

# Study bias from unmeasured confounder
results_confounded <- run_mc_sims(
  estimator = "CLPM",
  trials = 10,
  waves = 4,
  sample_size = 1000,
  data_generation = "clpmu"  # Data has confounder, but CLPM ignores it
)

# Custom parameter grid
my_grid <- expand.grid(
  stability_p = c(0.2, 0.5),
  stability_q = c(0.3, 0.6),
  cross_p = c(0.0, 0.1),
  cross_q = c(0.0, 0.1),
  variance_between_x = c(0.5, 1.0)
)

results <- run_mc_sims(
  estimator = "CLPM",
  param_grid = my_grid,
  trials = 10,
  waves = 4,
  sample_size = 1000,
  data_generation = "riclpm"
)

# Compare estimators under confounding
estimators <- c("OLS", "RICLPM", "CLPM", "RI")
all_results <- lapply(estimators, function(est) {
  run_mc_sims(estimator = est, trials = 5, waves = 3, data_generation = "clpmu")
})
combined_results <- do.call(rbind, all_results)
} # }
```
