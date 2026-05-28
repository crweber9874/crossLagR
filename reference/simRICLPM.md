# simRICLPM

Simulate data from a Random Intercept Cross-Lagged Panel Model (RICLPM)
with proper intercept structure

This function creates synthetic data for a Random Intercept Cross-Lagged
Panel Model (RICLPM) that matches the parameter structure and
identification constraints of the estimateRICLPM function. The model
separates stable between-person differences (random intercepts) from
within-person dynamics (autoregressive and cross-lagged effects).

## Usage

``` r
simRICLPM(
  waves = 5,
  sample.nobs = 1000,
  include_z = FALSE,
  beta_x = 0.3,
  beta_y = 0.3,
  beta_z = 0.3,
  omega_xy = 0.1,
  omega_yx = 0.1,
  omega_zx = 0.1,
  omega_zy = 0.1,
  omega_xz = 0.1,
  omega_yz = 0.1,
  var_p = 1,
  var_q = 1,
  var_r = 1,
  cov_pq = 0.1,
  cov_pr = 0.1,
  cov_qr = 0.1,
  var_BX = 1.5,
  var_BY = 1.5,
  var_BZ = 1.5,
  cov_BXBY = 0.5,
  cov_BXBZ = 0.5,
  cov_BYBZ = 0.5,
  mean_BX = 0,
  mean_BY = 0,
  mean_BZ = 0,
  constrain_beta = TRUE,
  constrain_omega = TRUE,
  constrain_residual_variances = TRUE,
  constrain_residual_covariances = TRUE,
  estimate_means = FALSE,
  ...
)
```

## Arguments

- waves:

  The number of waves (time points) in the model. Must be at least 2.

- sample.nobs:

  The number of observations (participants) to simulate. Default is
  1000.

- include_z:

  Logical. If TRUE, includes a third variable Z in the model. Default is
  FALSE.

- beta_x:

  Autoregressive effect for X (X(t-1) -\> X(t)). Default is 0.3.

- beta_y:

  Autoregressive effect for Y (Y(t-1) -\> Y(t)). Default is 0.3.

- beta_z:

  Autoregressive effect for Z (Z(t-1) -\> Z(t)), if included. Default is
  0.3.

- omega_xy:

  Cross-lagged effect of X on Y (X(t-1) -\> Y(t)). Default is 0.1.

- omega_yx:

  Cross-lagged effect of Y on X (Y(t-1) -\> X(t)). Default is 0.1.

- omega_zx:

  Cross-lagged effect of Z on X (Z(t-1) -\> X(t)), if Z included.
  Default is 0.1.

- omega_zy:

  Cross-lagged effect of Z on Y (Z(t-1) -\> Y(t)), if Z included.
  Default is 0.1.

- omega_xz:

  Cross-lagged effect of X on Z (X(t-1) -\> Z(t)), if Z included.
  Default is 0.1.

- omega_yz:

  Cross-lagged effect of Y on Z (Y(t-1) -\> Z(t)), if Z included.
  Default is 0.1.

- var_p:

  Residual variance for within-person X factors. Default is 1.

- var_q:

  Residual variance for within-person Y factors. Default is 1.

- var_r:

  Residual variance for within-person Z factors, if included. Default is
  1.

- cov_pq:

  Residual covariance between within-person X and Y factors. Default is
  0.1.

- cov_pr:

  Residual covariance between within-person X and Z factors, if Z
  included. Default is 0.1.

- cov_qr:

  Residual covariance between within-person Y and Z factors, if Z
  included. Default is 0.1.

- var_BX:

  Variance of random intercept for X. Default is 1.5.

- var_BY:

  Variance of random intercept for Y. Default is 1.5.

- var_BZ:

  Variance of random intercept for Z, if included. Default is 1.5.

- cov_BXBY:

  Covariance between random intercepts for X and Y. Default is 0.5.

- cov_BXBZ:

  Covariance between random intercepts for X and Z, if Z included.
  Default is 0.5.

- cov_BYBZ:

  Covariance between random intercepts for Y and Z, if Z included.
  Default is 0.5.

- mean_BX:

  Mean of random intercept for X. Default is 0.

- mean_BY:

  Mean of random intercept for Y. Default is 0.

- mean_BZ:

  Mean of random intercept for Z, if included. Default is 0.

- constrain_beta:

  Logical. If TRUE, constrains autoregressive effects to equality across
  waves. Default is TRUE.

- constrain_omega:

  Logical. If TRUE, constrains cross-lagged effects to equality across
  waves. Default is TRUE.

- constrain_residual_variances:

  Logical. If TRUE, constrains residual variances to equality across
  waves. Default is TRUE.

- constrain_residual_covariances:

  Logical. If TRUE, constrains residual covariances to equality across
  waves. Default is TRUE.

- estimate_means:

  Logical. If TRUE, estimates means for within-person factors. Default
  is FALSE.

- ...:

  Additional arguments to pass to lavaan::simulateData.

## Value

A list containing:

- `model`: The Lavaan model syntax used for data simulation

- `data`: The simulated data in a data frame format

- `parameters`: A list of the parameters used in simulation

## Details

The Random Intercept Cross-Lagged Panel Model (RICLPM) decomposes each
observed variable into: - A stable between-person component (random
intercept) - A time-varying within-person component (within-person
factor) - Measurement error (fixed to 0 for perfect indicators)

Variable naming convention: - X variables: x1, x2, x3, ..., x\[waves\] -
Y variables: y1, y2, y3, ..., y\[waves\] - Z variables (if included):
z1, z2, z3, ..., z\[waves\]

Model identification follows standard RICLPM practices: - Random
intercept means are estimated (unless estimate_means = FALSE) -
Within-person factor means are typically fixed to 0 (representing
deviations) - Random intercepts don't correlate with first wave
within-person factors - Observed variable residuals are fixed to 0
(perfect indicators)

## Examples

``` r
# Basic RICLPM simulation
riclpm_data <- simRICLPM(waves = 4, sample.nobs = 500)
#> ✅ Successfully simulated RICLPM data with 4 waves and 500 observations.

# RICLPM with stronger cross-lagged effects
riclpm_strong <- simRICLPM(
  waves = 5,
  omega_xy = 0.3,
  omega_yx = 0.3,
  sample.nobs = 1000
)
#> ✅ Successfully simulated RICLPM data with 5 waves and 1000 observations.

# RICLPM with third variable
riclpm_3var <- simRICLPM(
  waves = 4,
  include_z = TRUE,
  sample.nobs = 800
)
#> ✅ Successfully simulated RICLPM data with 4 waves (including Z variable) and 800 observations.

# RICLPM with unconstrained parameters
riclpm_free <- simRICLPM(
  waves = 4,
  constrain_beta = FALSE,
  constrain_omega = FALSE,
  sample.nobs = 600
)
#> Warning: lavaan->lavSimulateData():  
#>    some regression coefficients are unspecified and will be set to zero
#> ✅ Successfully simulated RICLPM data with 4 waves and 600 observations.
```
