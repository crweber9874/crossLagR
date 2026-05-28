# estimateLChange

Generate lavaan model syntax for univariate or bivariate Latent Change
Score Models (LCSM).

## Usage

``` r
estimateLChange(
  waves = 10,
  variable_type = c("univariate", "bivariate"),
  constrain_indicator_variances = TRUE,
  constrain_beta = TRUE,
  constrain_omega = TRUE,
  estimate_change_to_change = FALSE,
  constrain_change_to_change = FALSE,
  constrain_change_cross_lag = FALSE,
  estimate_constant_change = TRUE
)
```

## Arguments

- waves:

  Integer. Number of waves (time points). Must be \>= 3.

- variable_type:

  Character. Either `"univariate"` or `"bivariate"`. Default is
  `"univariate"`.

- constrain_indicator_variances:

  Logical. If TRUE, constrains observed indicator residual variances
  (unique factors) to equality across waves. Default is TRUE.

- constrain_beta:

  Logical. If TRUE, constrains proportional (autoregressive) effects to
  equality across waves. Default is TRUE.

- constrain_omega:

  Logical. If TRUE, constrains cross-lagged coupling parameters to
  equality across waves (bivariate only). Default is TRUE.

- estimate_change_to_change:

  Logical. If TRUE, includes change-to-change autoregressive and
  cross-lagged effects. Requires \>= 4 waves. Default is FALSE.

- constrain_change_to_change:

  Logical. If TRUE, constrains change-to-change autoregressive effects
  to equality across waves. Default is FALSE.

- constrain_change_cross_lag:

  Logical. If TRUE, constrains change-to-change cross-lagged effects to
  equality across waves. Default is FALSE.

- estimate_constant_change:

  Logical. If TRUE, estimates a second-order constant change
  (accumulating) factor that loads on all latent change scores. Default
  is TRUE.

## Value

A character string containing lavaan model syntax for the specified
LCSM.

## Details

This function generates lavaan syntax for a Latent Change Score Model
following the unified framework of Usami, Murayama, & Hamaker (2019).
The LCS model uses the measurement equation (separating latent true
scores from unique factors) and the dynamic equation (with accumulating
factors and lagged regression).

Parameter labels follow the unified naming convention:

- `ar_x`, `ar_y`: Proportional change coefficients. In the LCS, the
  effective autoregressive parameter is \\1 + ar\\, so `ar_x` is the
  proportional self-feedback (denoted \\\beta_x\\ in Usami et al.).

- `cl_yx`, `cl_xy`: Coupling parameters (cross-lagged effects). `cl_yx`
  = effect of Y level on X change; `cl_xy` = effect of X level on Y
  change. Referred to as \\\gamma\\ in Usami et al.

- `u_var_x`, `u_var_y`: Unique factor (measurement error) variances.

- `A_x`, `A_y`: Accumulating factors (constant change factors). Unlike
  growth factors (I, S) in the LCM-SR which have only direct effects,
  these accumulating factors have both direct and indirect effects on
  scores through the lagged relations.

Change-to-change parameters (`phi_x`, `phi_y`, `phi_cl_x`, `phi_cl_y`)
are unique to the LCS model and have no direct analogue in other models
in the unified framework.

## References

Usami, S., Murayama, K., & Hamaker, E. L. (2019). A unified framework of
longitudinal models to examine reciprocal relations. *Psychological
Methods*, 24(5), 637-657.

McArdle, J. J., & Hamagami, F. (2001). Latent difference score
structural models for linear dynamic analyses with incomplete
longitudinal data.

## Examples

``` r
if (FALSE) { # \dontrun{
# Univariate LCSM with 5 waves
syntax <- estimateLChange(waves = 5)
cat(syntax)

# Bivariate LCSM with coupling parameters
syntax <- estimateLChange(waves = 5, variable_type = "bivariate")

# Fit with lavaan
library(lavaan)
fit <- lavaan(syntax, data = my_data, meanstructure = TRUE)
summary(fit)
} # }
```
