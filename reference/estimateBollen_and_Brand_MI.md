# estimateBollen_and_Brand_MI

Generates lavaan syntax for a multi-indicator Bollen & Brand (2010)
dynamic panel model.

## Usage

``` r
estimateBollen_and_Brand_MI(
  waves = 5,
  n_indicators_x = 3,
  n_indicators_y = 3,
  invariance = c("scalar", "metric", "configural", "strict"),
  correlated_residuals = TRUE,
  x_effect = "lagged",
  x_autoregression = TRUE,
  y_effect_on_x = TRUE,
  constrain_coefficients = FALSE,
  constrain_residual_variances = TRUE
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

- x_effect:

  Character. `"none"`, `"concurrent"`, or `"lagged"`. Default
  `"lagged"`.

- x_autoregression:

  Logical. If TRUE, x latents have AR paths and a latent fixed effect
  I_x. Default `TRUE`.

- y_effect_on_x:

  Logical. If TRUE, reciprocal y_t-1 -\> x_t. Default `TRUE`.

- constrain_coefficients:

  Logical. Equality constraints across waves.

- constrain_residual_variances:

  Logical.

## Value

Character string of lavaan syntax.

## Details

Latent wave factors `p_w` (X) and `q_w` (Y) with multi-indicators
replace the manifest x_w/y_w of the standard Bollen-Brand model. Latent
fixed effects `I_y`, `I_x` load with 1 on the endogenous wave latents
(`q2..qT` and `p2..pT` respectively). Same-item residuals are correlated
across waves.

## References

Bollen, K. A., & Brand, J. E. (2010). *Social Forces*, 89(1), 1-34.

Dishop, C. R., & DeShon, R. P. (2018). *Psychological Methods*, 23(4),
1089-1112.

Usami, S., Murayama, K., & Hamaker, E. L. (2019). *Psychological
Methods*, 24(5), 637-657.
