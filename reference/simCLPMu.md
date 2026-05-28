# simCLPMu

Simulate data from a cross-lagged panel model (CLPM) with a time-variant
confounder.

## Usage

``` r
simCLPMu(
  waves = 10,
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
  ...
)
```

## Arguments

- waves:

  The number of waves (time points) in the model.

- stability_p:

  The stability parameter for the x variable (autoregressive effect).

- stability_q:

  The stability parameter for the y variable (autoregressive effect).

- cross_p:

  The cross-lagged effect of y on x at the next time point.

- cross_q:

  The cross-lagged effect of x on y at the next time point.

- variance_p:

  The variance of the p latent variable.

- variance_q:

  The variance of the q latent variable.

- cov_pq:

  The covariance between x and y within the same time point.

- confounder_p:

  The effect of the confounder u on x variables.

- confounder_q:

  The effect of the confounder u on y variables.

- confounder_variance:

  The variance of the confounder u.

- confounder_stability:

  The stability parameter for the confounder (autoregressive effect).

- ...:

  Additional arguments to pass to the \`lavaan::simulateData\` function.

## Value

A list containing two elements: \* \`model\`: The Lavaan model syntax
used for data simulation. \* \`data\`: The simulated data in a data
frame format.
