# monteCarloConfounder

Monte Carlo simulation for Cross-Lagged Panel Model with unmeasured
confounder. Tests how different estimators perform when data is
generated with an unmeasured confounder that affects both X and Y
variables across time.

## Usage

``` r
monteCarloConfounder(
  trials = 10,
  waves = 5,
  estimator = "RICLPM",
  stability_p = 0.2,
  stability_q = 0.2,
  cross_p = 0.1,
  cross_q = 0.1,
  variance_p = 1,
  variance_q = 1,
  cov_pq = 0.1,
  confounder_p = 0.3,
  confounder_q = 0.3,
  confounder_variance = 1,
  confounder_stability = 0.4,
  sample_size = 2500,
  include_confounder = TRUE,
  ...
)
```

## Arguments

- trials:

  The number of trials for the Monte Carlo simulation.

- waves:

  The number of waves (time points) in the model.

- estimator:

  The estimation method to use. Options: "OLS", "RICLPM", "CLPM".

- stability_p:

  The stability parameter for the X variable (autoregressive effect).

- stability_q:

  The stability parameter for the Y variable (autoregressive effect).

- cross_p:

  The cross-lagged effect of Y on X at the next time point.

- cross_q:

  The cross-lagged effect of X on Y at the next time point.

- variance_p:

  The variance of the p latent variable.

- variance_q:

  The variance of the q latent variable.

- cov_pq:

  The covariance between X and Y within the same time point.

- confounder_p:

  The effect of the confounder on X variables.

- confounder_q:

  The effect of the confounder on Y variables.

- confounder_variance:

  The variance of the confounder.

- confounder_stability:

  The stability parameter for the confounder (autoregressive effect).

- sample_size:

  The sample size for each simulation.

- include_confounder:

  Logical. Whether to include the confounder in data generation.

- ...:

  Additional arguments.

## Value

A data frame containing the results of the simulation with columns for:
estimated coefficients, true parameters, trial number, estimator used,
and bias measures.

## Details

This function generates data using simCLPMu (which includes an
unmeasured confounder) and then fits models that may or may not account
for this confounder. This allows testing of bias introduced by
unmeasured confounders.

The function tests the following scenario: 1. Data is generated with an
unmeasured confounder affecting both X and Y 2. Models are fitted that
ignore this confounder 3. Bias in parameter estimates is measured

## Examples

``` r
if (FALSE) { # \dontrun{
# Test bias from unmeasured confounder
results <- monteCarloConfounder(
  trials = 100,
  waves = 4,
  estimator = "RICLPM",
  confounder_p = 0.3,
  confounder_q = 0.3,
  sample_size = 1000
)

# Compare multiple estimators
estimators <- c("OLS", "RICLPM", "CLPM")
results_list <- lapply(estimators, function(est) {
  monteCarloConfounder(trials = 50, estimator = est)
})
} # }
```
