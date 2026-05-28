# simpleSim

Monte Carlo simulation comparing OLS and CLPM estimators. Allows the
user to specify different data-generating conditions (CLPM or RI-CLPM),
and apply either an OLS or CLPM model to the data.

## Usage

``` r
monteCarlo_OLS(
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
