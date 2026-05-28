# estimateLCMSR_MI

Generates lavaan syntax for a multi-indicator bivariate Latent Curve
Model with Structured Residuals (LCM-SR / GCLM).

## Usage

``` r
estimateLCMSR_MI(
  waves = 5,
  n_indicators_x = 3,
  n_indicators_y = 3,
  invariance = c("scalar", "metric", "configural", "strict"),
  correlated_residuals = TRUE,
  time_scores = NULL,
  constrain_beta = TRUE,
  constrain_omega = TRUE,
  constrain_residual_variances = TRUE,
  constrain_residual_covariances = TRUE,
  estimate_quadratic = FALSE
)
```

## Arguments

- waves:

  Integer. Number of waves (\>= 3).

- n_indicators_x, n_indicators_y:

  Integer. Indicators per wave.

- invariance:

  Character. One of `"configural"`, `"metric"`, `"scalar"`, `"strict"`.
  Default `"scalar"`.

- correlated_residuals:

  Logical. Default `TRUE`.

- time_scores:

  Numeric vector length `waves`. Default `0:(waves-1)`.

- constrain_beta, constrain_omega, constrain_residual_variances,
  constrain_residual_covariances:

  Logical.

- estimate_quadratic:

  Logical. Add quadratic growth factor (requires `waves >= 4`). Default
  `FALSE`.

## Value

Character string of lavaan syntax.

## Details

Curve-of-factors layout: indicators -\> wave latent `eta_x_w`; growth
factors `I_x`, `S_x` load on `eta_x_w` (direct effects only), and `p_w`
captures the structured residual. AR/CL operate on `p_w`, `q_w`.
Same-item residuals are correlated across waves.

## References

Curran, P. J., Howard, A. L., Bainter, S. A., Lane, S. T., & McGinley,
J. S. (2014). *Annual Review of Psychology*, 65, 637-660.

Usami, S., Murayama, K., & Hamaker, E. L. (2019). *Psychological
Methods*, 24(5), 637-657.
