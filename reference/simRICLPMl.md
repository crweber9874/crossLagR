# simulateRICLPMl

Simulate data from a multiple indicator Random Intercept Cross-Lagged
Regression (RI-CLPM) model, where the first indicator of each latent
variable has a fixed factor loading of 1, and subsequent indicators have
factor loadings based on the provided parameters. Includes an optional
third variable that can interact with X and Y, also following a cross
lagged structure.

## Usage

``` r
simRICLPMl(
  waves = 10,
  num_indicators_x = 2,
  num_indicators_y = 2,
  num_indicators_z = NULL,
  stability_p = 0.2,
  stability_q = 0.2,
  stability_r = NULL,
  cross_p = 0.1,
  cross_q = 0.1,
  cross_zx = NULL,
  cross_zy = NULL,
  factor_loading_x = 1,
  factor_loading_y = 1,
  factor_loading_z = 1,
  residual_variance_x = 0.5,
  residual_variance_y = 0.5,
  residual_variance_z = 0.5,
  variance_p = 1,
  variance_q = 1,
  variance_r = NULL,
  cov_pq = 0.1,
  cov_pr = NULL,
  cov_qr = NULL,
  variance_between_x = 1,
  variance_between_y = 1,
  variance_between_z = NULL,
  cov_between_xy = 0.5,
  cov_between_xz = NULL,
  cov_between_yz = NULL,
  sample_size = 1000,
  ...
)
```

## Arguments

- waves:

  The number of waves (time points) in the model.

- num_indicators_x:

  The number of indicators for the X variable.

- num_indicators_y:

  The number of indicators for the Y variable.

- num_indicators_z:

  The number of indicators for the Z variable (if specified by user).

- stability_p:

  The stability parameter for the latent X variable. This is the
  autoregressive parameter for X. This should be specified with care,
  particularly when comparing a CLPM to RI-CLPM.

- stability_q:

  The stability parameter for the latent Y variable. This is the
  autoregressive parameter. This should be specified with care,
  particularly when comparing a CLPM to RI-CLPM.

- stability_r:

  The stability parameter for the latent Z variable (autoregressive
  effect, if applicable).

- cross_p:

  The cross-lagged effect of latent Y on latent X at the next time
  point.

- cross_q:

  The cross-lagged effect of latent X on latent Y at the next time
  point.

- cross_zx:

  The cross-lagged effect of latent Z on latent X at the next time point
  (if applicable).

- cross_zy:

  The cross-lagged effect of latent Z on latent Y at the next time point
  (if applicable).

- factor_loading_x:

  The factor loadings for the indicators of X. The first loading is
  implicitly 1, and subsequent loadings can be a single value or a
  vector of length \`num_indicators_x - 1\`.

- factor_loading_y:

  The factor loadings for the indicators of Y. The first loading is
  implicitly 1, and subsequent loadings can be a single value or a
  vector of length \`num_indicators_y - 1\`.

- factor_loading_z:

  The factor loadings for the indicators of Z (if applicable). The first
  loading is implicitly 1, and subsequent loadings can be a single value
  or a vector of length \`num_indicators_z - 1\`.

- residual_variance_x:

  The residual variances for the indicators of X (can be a single value
  or a vector of length \`num_indicators_x\`).

- residual_variance_y:

  The residual variances for the indicators of Y (can be a single value
  or a vector of length \`num_indicators_y\`).

- residual_variance_z:

  The residual variances for the indicators of Z (if applicable, can be
  a single value or a vector of length \`num_indicators_z\`).

- variance_p:

  The variance for the latent X variables.

- variance_q:

  The variance for the latent Y variables.

- variance_r:

  The variance for the latent Z variables (if applicable).

- cov_pq:

  The covariance between latent X and latent Y within the same time
  point, or wave

- cov_pr:

  The covariance between latent X and latent Z within the same time
  point (if specified by user).

- cov_qr:

  The covariance between latent Y and latent Z within the same time
  point (if specified by user).

- variance_between_x:

  The variance for the random intercept of X.

- variance_between_y:

  The variance for the random intercept of Y.

- variance_between_z:

  The variance for the random intercept of Z (if specified by user).

- cov_between_xy:

  The covariance of intercept terms between X and Y.

- cov_between_xz:

  The covariance of intercept terms between X and Z (if applicable).

- cov_between_yz:

  The covariance of intercept terms between Y and Z (if applicable).

- ...:

  Additional arguments to pass to the \`lavaan::simulateData\` function.

## Value

A list containing two elements: \* \`model\`: The lavaan model syntax
used for the data simulation. \* \`data\` : The simulated data in a data
frame format.
