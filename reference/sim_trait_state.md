# sim_trait_state

Simulate panel data from a trait-plus-state data-generating process.

## Usage

``` r
sim_trait_state(
  waves = 5,
  n = 1000,
  beta_x = 0,
  beta_y = 0,
  omega_xy = 0,
  omega_yx = 0,
  var_p = 1,
  var_q = 1,
  cov_pq = 0,
  var_BX = 0,
  var_BY = 0,
  cov_BXBY = 0,
  cor_trait_initial = 0
)
```

## Arguments

- waves:

  Integer. Number of waves (time points). Must be at least 2.

- n:

  Integer. Sample size (number of individuals).

- beta_x:

  Numeric. Autoregressive effect for X.

- beta_y:

  Numeric. Autoregressive effect for Y.

- omega_xy:

  Numeric. Cross-lagged effect from X(t-1) to Y(t).

- omega_yx:

  Numeric. Cross-lagged effect from Y(t-1) to X(t).

- var_p:

  Numeric. Innovation variance for within-person X process.

- var_q:

  Numeric. Innovation variance for within-person Y process.

- cov_pq:

  Numeric. Innovation covariance between within-person X and Y
  processes.

- var_BX:

  Numeric. Between-person trait variance for X.

- var_BY:

  Numeric. Between-person trait variance for Y.

- cov_BXBY:

  Numeric. Between-person trait covariance between X and Y.

- cor_trait_initial:

  Numeric. Correlation between traits and wave-1 within-person states.
  Default is 0 (independent, which satisfies the RI-CLPM identification
  assumption). Non-zero values violate this assumption and can reveal
  bias in models that assume independence. Applied symmetrically:
  cor(eta_x, wp_x1) = cor(eta_y, wp_y1) = cor_trait_initial.

## Value

A wide-format data frame with columns `x1, ..., x[waves]` and
`y1, ..., y[waves]`.

## Details

The function simulates two stable traits (between-person random
intercept components) and two within-person AR(1) processes with
optional cross-lagged effects. Observed scores are the sum of trait and
state components.

When `cor_trait_initial = 0` (default), traits and wave-1 states are
drawn independently. When non-zero, they are drawn jointly from a
multivariate normal, so individuals with higher trait values tend to
also have higher (or lower) initial state values. This violates the
RI-CLPM identification constraint that random intercepts are
uncorrelated with wave-1 within-person factors.

## Examples

``` r
dat <- sim_trait_state(
  waves = 5,
  n = 500,
  beta_x = 0.4,
  beta_y = 0.3,
  omega_xy = 0.1,
  omega_yx = 0.1,
  var_p = 1,
  var_q = 1,
  cov_pq = 0.2,
  var_BX = 0.6,
  var_BY = 0.6,
  cov_BXBY = 0.3
)
head(dat)
#>            x1          x2         x3         x4         x5         y1
#> 1  0.03547892 -0.48795941  0.3523056  0.2333733  1.2112398 -0.9805910
#> 2  0.19191876 -1.36091727 -0.7200436  2.0932643  1.3891585 -0.2431305
#> 3  2.35598001  1.17648541  2.2566169  0.2655163  0.2635006  0.4932567
#> 4 -0.06751661 -0.74718624  0.9432735 -1.1902808 -0.4423140  0.3616781
#> 5  2.45776350  0.26492281  0.5769836  1.8088350  0.9738890 -0.4174472
#> 6  0.63355370  0.07654207  1.7552491  2.3540461  1.9243578 -0.1459554
#>           y2          y3         y4         y5
#> 1 0.38356978  0.70071667  0.7323103 -0.7002522
#> 2 0.62666155 -0.14631416 -1.1615107 -0.1941928
#> 3 0.52370066  1.30324349 -0.3022192  0.7673734
#> 4 0.14990173  0.08921588 -0.7434546 -1.3503778
#> 5 0.72886858  1.55767179  1.0845610  0.1010938
#> 6 0.04351836 -0.09646026  0.3803094  0.2847826
```
