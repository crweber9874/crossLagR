# estimateAllisonChamberlainFI

Generate lavaan model syntax for the Allison-Chamberlain fixed-effects
SEM approach to cross-lagged panel models.

## Usage

``` r
estimateAllisonChamberlainFI(
  waves = 5,
  constrain_coefficients = TRUE,
  constrain_residual_variances = TRUE,
  correlate_errors_with_future = TRUE,
  correlate_alpha_with_exogenous = TRUE
)
```

## Arguments

- waves:

  Integer. Number of waves (time points). Must be \>= 3.

- constrain_coefficients:

  Logical. If TRUE, constrains autoregressive and cross-lagged effects
  to equality across waves. Default is TRUE.

- constrain_residual_variances:

  Logical. If TRUE, constrains residual variances to equality across
  waves (after wave 1). Default is TRUE.

- correlate_errors_with_future:

  Logical. If TRUE, allows residuals of each equation to be correlated
  with future values of the opposite predetermined variable, following
  the Chamberlain correction. Default is TRUE.

- correlate_alpha_with_exogenous:

  Logical. If TRUE, allows the latent fixed-effect variables (I_x, I_y)
  to be correlated with all observed x and y values. Default is TRUE.

## Value

A character string containing lavaan model syntax for the
Allison-Chamberlain fixed-effects SEM model.

## Details

This function generates lavaan syntax for the Allison-Chamberlain
fixed-effects SEM following the unified framework of Usami, Murayama, &
Hamaker (2019). This model includes latent fixed effects (analogous to
trait factors I) and relaxes strict exogeneity to predetermination.

Parameter labels follow the unified naming convention:

- `ar_x`, `ar_y`: Autoregressive effects. These operate on observed
  scores with individual heterogeneity controlled via latent fixed
  effects.

- `cl_xy`: Cross-lagged effect of X on Y.

- `cl_yx`: Cross-lagged effect of Y on X.

- `d_var_x`, `d_var_y`: Dynamic residual variances.

- `I_x`, `I_y`: Latent fixed effects (individual heterogeneity). These
  function as trait factors in the unified framework, with unit loadings
  on endogenous outcomes.

## References

Allison, P. D. (2005). *Fixed Effects Regression Methods for
Longitudinal Data Using SAS*. SAS Institute.

Chamberlain, G. (1982). Multivariate regression models for panel data.
*Journal of Econometrics*, 18(1), 5-46.

Usami, S., Murayama, K., & Hamaker, E. L. (2019). A unified framework of
longitudinal models to examine reciprocal relations. *Psychological
Methods*, 24(5), 637-657.

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate Allison-Chamberlain FI syntax with 5 waves
syntax <- estimateAllisonChamberlainFI(waves = 5)
cat(syntax)

# Fit with lavaan
library(lavaan)
fit <- lavaan::sem(syntax, data = my_data)
summary(fit, fit.measures = TRUE, standardized = TRUE)
} # }
```
