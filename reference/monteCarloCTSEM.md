# monteCarloCTSEM

A somewhat strange function. It's common to estimate the cross-lagged
panel model with two OLS regression models. allows the user to specify
different Data Generating conditions, and apply the OLS model to the
data. The output is a data frame that includes the estimated
coefficients across these trials

## Usage

``` r
monteCarloCTSEM(
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
  cores = 10,
  confounder_p = 0.3,
  confounder_q = 0.3,
  confounder_variance = 1,
  confounder_stability = 0.4,
  include_confounder = TRUE,
  ...
)
```

## Arguments

- trials:

  The number of trials for the Monte Carlo simulation.

- waves:

  The number of waves (time points) in the model.

- model:

  Specify whether an "ols" or "clpm" model.

- data_generation:

  Specify the data generation method: "clpm", "ri-clpm"

- proportion_change_x:

  Proportion of change in x.

- proportion_change_y:

  Proportion of change in y.

## Value

A data frame containing the results of the simulation.
