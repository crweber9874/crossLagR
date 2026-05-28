# monteCarloTraitState

Monte Carlo simulation that varies the trait-to-state variance ratio and
fits four competing panel models: Bollen & Brand, RI-CLPM, CLPM, and
Latent Change Score Model.

Data are generated from a trait-plus-state DGP via
[`sim_trait_state()`](https://crweber9874.github.io/crossLagR/reference/sim_trait_state.md).
The key manipulation is the ratio of between-person (trait) variance to
within-person (state) variance. When trait variance is 0 the DGP is a
pure CLPM; as trait variance increases the data increasingly violate
CLPM assumptions.

## Usage

``` r
monteCarloTraitState(
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
)
```

## Arguments

- trials:

  Integer. Number of Monte Carlo replications per condition. Default is
  100.

- waves:

  Integer. Number of waves. Must be \>= 5 for latent change model.
  Default is 5.

- n:

  Integer. Sample size per replication. Default is 1000.

- beta_x:

  Numeric. Autoregressive effect for X within-person process. Default is
  0.3.

- beta_y:

  Numeric. Autoregressive effect for Y within-person process. Default is
  0.3.

- omega_xy:

  Numeric. Cross-lagged effect from X to Y (within-person). Default is
  0.1.

- omega_yx:

  Numeric. Cross-lagged effect from Y to X (within-person). Default is
  0.1.

- var_p:

  Numeric. Within-person innovation variance for X. Default is 1.

- var_q:

  Numeric. Within-person innovation variance for Y. Default is 1.

- cov_pq:

  Numeric. Within-person innovation covariance. Default is 0.1.

- trait_ratio:

  Numeric vector. Ratio(s) of trait variance to state variance. For each
  value r, var_BX = var_BY = r \* var_p. Default is c(0, 0.5, 1, 2).

- cov_BXBY:

  Numeric or NULL. Trait covariance. If NULL, set to 0.3 \* var_BX.
  Default is NULL.

- verbose:

  Logical. Print progress. Default is TRUE.

- ...:

  Additional arguments (currently unused).

## Value

A data frame with one row per trial x condition x estimator, containing:

- trait_ratio:

  The trait-to-state variance ratio used for simulation

- var_BX, var_BY:

  The actual trait variances used

- estimator:

  Which model was fit: "clpm", "riclpm", "bollen_brand", "lchange"

- beta_x_est, beta_y_est:

  Estimated autoregressive / proportional effects

- omega_xy_est, omega_yx_est:

  Estimated cross-lagged / coupling effects

- converged:

  Logical. Did the model converge?

- trial:

  Trial number

## Details

For a single sample (trials = 1), this can be used as a one-shot
demonstration. For Monte Carlo (trials \> 1), bias and RMSE can be
computed by comparing estimates to the true DGP values.

The four estimators extract different types of parameters:

- **CLPM**: autoregressive (beta) and cross-lagged (omega) on observed
  composites

- **RI-CLPM**: autoregressive and cross-lagged on within-person latent
  factors

- **Bollen & Brand**: autoregressive (rho) and cross-lagged (b1, b2)
  with latent fixed effects

- **Latent Change**: proportional (beta) and coupling (omega) effects on
  change scores

## Examples

``` r
if (FALSE) { # \dontrun{
# Single-sample demonstration
demo <- monteCarloTraitState(trials = 1, trait_ratio = c(0, 1, 2), n = 500)
demo[, c("trait_ratio", "estimator", "beta_x_est", "omega_xy_est")]

# Full Monte Carlo
mc <- monteCarloTraitState(trials = 100, trait_ratio = c(0, 0.5, 1, 2), n = 1000)

# Compute bias
library(dplyr)
mc %>%
  filter(converged) %>%
  group_by(trait_ratio, estimator) %>%
  summarise(
    bias_omega_xy = mean(omega_xy_est - 0.1, na.rm = TRUE),
    bias_omega_yx = mean(omega_yx_est - 0.1, na.rm = TRUE),
    rmse_omega_xy = sqrt(mean((omega_xy_est - 0.1)^2, na.rm = TRUE))
  )
} # }
```
