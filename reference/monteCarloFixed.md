# monteCarloFixed

Monte Carlo simulation for Fixed Individual Effects estimation Allows
the user to specify different Data Generating conditions, and apply the
FI model to the data. The output is a data frame that includes the
estimated coefficients across these trials

## Usage

``` r
monteCarloFixed(
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
  ...
)
```

## Arguments

- trials:

  The number of trials for the Monte Carlo simulation.

- waves:

  The number of waves (time points) in the model.

- variance_between_x:

  Between-person variance for x variable (for RICLPM data)

- variance_between_y:

  Between-person variance for y variable (for RICLPM data)

- stability_q:

  Autoregressive effect for q (y) variable

- stability_p:

  Autoregressive effect for p (x) variable

- cross_p:

  Cross-lagged effect of y on x

- cross_q:

  Cross-lagged effect of x on y

- variance_p:

  Within-person variance for p variable

- variance_q:

  Within-person variance for q variable

- sample_size:

  Sample size for simulation

- dgp:

  Specify the data generation method: "riclpm", "clpm", "clpmu"

- confounder_p:

  Effect of confounder on x variables (for clpmu)

- confounder_q:

  Effect of confounder on y variables (for clpmu)

- confounder_variance:

  Variance of the confounder (for clpmu)

- confounder_stability:

  Autoregressive effect of confounder (for clpmu)

- include_confounder:

  Whether to include confounder in clpmu model

- confounder_type:

  Character string specifying confounder type for clpmu dgp. Must be one
  of "time_variant" (default) or "time_invariant".

- ...:

  Additional arguments

## Value

A data frame containing the results of the simulation.
