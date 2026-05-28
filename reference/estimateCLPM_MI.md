# estimateCLPM_MI

Generates lavaan syntax for a multi-indicator Cross-Lagged Panel Model
(CLPM) following Usami, Murayama, & Hamaker (2019).

## Usage

``` r
estimateCLPM_MI(
  waves = 5,
  n_indicators_x = 3,
  n_indicators_y = 3,
  invariance = c("scalar", "metric", "configural", "strict"),
  correlated_residuals = TRUE,
  constrain_beta = TRUE,
  constrain_omega = TRUE,
  constrain_residual_variances = TRUE,
  constrain_residual_covariances = TRUE,
  estimate_means = TRUE,
  start_values = FALSE
)
```

## Arguments

- waves:

  Integer. Number of waves (\>= 2).

- n_indicators_x:

  Integer. Indicators for X at each wave (\>= 1).

- n_indicators_y:

  Integer. Indicators for Y at each wave (\>= 1).

- invariance:

  Character. One of `"configural"`, `"metric"`, `"scalar"`, `"strict"`.
  Default `"scalar"`.

- correlated_residuals:

  Logical. Correlate same-item residuals across waves. Default `TRUE`.

- constrain_beta, constrain_omega, constrain_residual_variances,
  constrain_residual_covariances:

  Logical. Equality constraints across waves on AR, CL, dynamic residual
  variances and covariances respectively.

- estimate_means:

  Logical. Estimate latent factor means at wave 1.

- start_values:

  Logical. Provide starting values for AR/CL.

## Value

Character string of lavaan model syntax.

## Details

Variables follow naming `x{w}_{j}`, `y{w}_{j}`. The first indicator is
the marker (loading fixed to 1, intercept fixed to 0). With
`invariance = "scalar"` the remaining loadings and intercepts are
constrained equal across waves. `"strict"` additionally constrains
indicator residual variances across waves. Same-item residuals are
correlated across all wave-pairs when `correlated_residuals = TRUE` to
absorb stable item-specific variance.

## References

Usami, S., Murayama, K., & Hamaker, E. L. (2019). *Psychological
Methods*, 24(5), 637-657.

## Examples

``` r
if (FALSE) { # \dontrun{
syntax <- estimateCLPM_MI(waves = 4, n_indicators_x = 3, n_indicators_y = 3)
fit <- lavaan::lavaan(syntax, data = my_data, meanstructure = TRUE,
                      int.ov.free = TRUE, int.lv.free = TRUE)
} # }
```
