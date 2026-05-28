# estimateALT_MI

Generates lavaan syntax for a multi-indicator bivariate Autoregressive
Latent Trajectory (ALT) model.

## Usage

``` r
estimateALT_MI(
  waves = 5,
  n_indicators_x = 3,
  n_indicators_y = 3,
  invariance = c("scalar", "metric", "configural", "strict"),
  correlated_residuals = TRUE,
  time_scores = NULL,
  constrain_beta = TRUE,
  constrain_omega = TRUE,
  constrain_residual_variances = TRUE,
  constrain_residual_covariances = TRUE
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

## Value

Character string of lavaan syntax.

## Details

Wave-specific latent factors `p_w`, `q_w` are defined directly from
multiple indicators. Growth factors `I_x`, `S_x`, `I_y`, `S_y` are
accumulating: they load on waves 2..T only, so they have direct effects
on endogenous waves AND indirect effects through the AR/CL dynamics.
Wave 1 is exogenous. Same-item residuals are correlated across waves.

## References

Bollen, K. A., & Curran, P. J. (2006). *Latent curve models*. Wiley.

Curran, P. J., & Bollen, K. A. (2001). The best of both worlds. In
Collins & Sayer (Eds.), *New methods for the analysis of change*
(107-135). APA.

Usami, S., Murayama, K., & Hamaker, E. L. (2019). *Psychological
Methods*, 24(5), 637-657.
