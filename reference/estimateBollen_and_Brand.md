# estimateBollen_and_Brand

Generate lavaan model syntax for the Bollen and Brand (2010) dynamic
panel model, as described in Dishop & DeShon (2018).

## Usage

``` r
estimateBollen_and_Brand(
  waves = 5,
  x_effect = "lagged",
  x_autoregression = TRUE,
  y_effect_on_x = TRUE,
  constrain_coefficients = FALSE,
  constrain_residual_variances = TRUE
)
```

## Arguments

- waves:

  Integer. Number of waves (time points). Must be \>= 3.

- x_effect:

  Character. How x influences y: `"none"` (single variable
  autoregressive model), `"concurrent"` (x_t predicts y_t), or
  `"lagged"` (x_t-1 predicts y_t). Default is `"lagged"`.

- x_autoregression:

  Logical. If TRUE, x has autoregressive paths and gets its own latent
  fixed effect (I_x). When FALSE, x is treated as fully exogenous.
  Default is FALSE.

- y_effect_on_x:

  Logical. If TRUE, includes reciprocal cross-lagged paths from y to x
  (y_t-1 predicts x_t). Requires `x_autoregression = TRUE`. Default is
  FALSE.

- constrain_coefficients:

  Logical. If TRUE, constrains autoregressive and cross-lagged effects
  to equality across waves (stationarity). Default is TRUE.

- constrain_residual_variances:

  Logical. If TRUE, constrains residual variances to equality across
  waves. Default is TRUE.

## Value

A character string containing lavaan model syntax.

## Details

This function generates lavaan syntax for the Bollen & Brand dynamic
panel model following the unified framework of Usami, Murayama, &
Hamaker (2019). This model is closely related to the ALT model, using
accumulating factors with lagged regression. It conditions on the first
observation and treats it as exogenous.

Parameter labels follow the unified naming convention:

- `ar_y`, `ar_x`: Autoregressive effects. In this model these operate on
  observed scores (no within/between decomposition), similar to the ALT
  model.

- `cl_xy`: Cross-lagged effect of X on Y.

- `cl_yx`: Cross-lagged effect of Y on X (reciprocal, if enabled).

- `d_var_y`, `d_var_x`: Dynamic residual variances.

- `I_y`, `I_x`: Latent fixed effects (individual heterogeneity),
  functioning as accumulating factors with unit loadings on endogenous
  outcomes. In the unified framework these are analogous to the A factor
  in the ALT/LCS models.

## References

Bollen, K. A., & Brand, J. E. (2010). A general panel model with random
and fixed effects: A structural equations approach. *Social Forces*,
89(1), 1-34.

Dishop, C. R., & DeShon, R. P. (2018). A tutorial on Bollen and Brand's
approach to modeling dynamics. *Psychological Methods*, 23(4),
1089-1112.

Usami, S., Murayama, K., & Hamaker, E. L. (2019). A unified framework of
longitudinal models to examine reciprocal relations. *Psychological
Methods*, 24(5), 637-657.

## Examples

``` r
if (FALSE) { # \dontrun{
# Reciprocal cross-lag dynamic panel (Dishop & DeShon Figure 5)
syntax <- estimateBollen_and_Brand(waves = 5, x_effect = "lagged",
                                    x_autoregression = TRUE, y_effect_on_x = TRUE)
cat(syntax)

# Fit with lavaan
library(lavaan)
fit <- lavaan::sem(syntax, data = my_wide_data)
summary(fit, fit.measures = TRUE, standardized = TRUE)
} # }
```
