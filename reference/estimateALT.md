# estimateALT, following Usami, Murayama, & Hamaker (2019) unified frameowk.

Generates lavaan model syntax for the bivariate Autoregressive Latent
Trajectory (ALT) model.

## Usage

``` r
estimateALT(
  waves = 5,
  time_scores = NULL,
  constrain_beta = TRUE,
  constrain_omega = TRUE,
  constrain_residual_variances = TRUE,
  constrain_residual_covariances = TRUE
)
```

## Arguments

- waves:

  Integer. Number of waves (time points). Must be \>= 3.

- time_scores:

  Numeric vector of time scores for slope factor loadings. If NULL, uses
  0, 1, 2, ..., waves-1. Default is NULL.

- constrain_beta:

  Logical. If TRUE, constrains autoregressive effects to equality across
  waves. Default is TRUE.

- constrain_omega:

  Logical. If TRUE, constrains cross-lagged effects to equality across
  waves. Default is TRUE.

- constrain_residual_variances:

  Logical. If TRUE, constrains dynamic residual variances to equality
  across waves. Default is TRUE.

- constrain_residual_covariances:

  Logical. If TRUE, constrains dynamic residual covariances to equality
  across waves. Default is TRUE.

## Value

A character string containing lavaan model syntax.

## Details

The ALT model combines latent growth trajectories (intercept I, slope S)
with autoregressive and cross-lagged dynamics, following the unified
framework of Usami, Murayama, & Hamaker (2019). Unlike the LCM-SR (which
models dynamics in the residuals from growth curves), the ALT model
applies AR/CL dynamics directly to the latent scores. The growth factors
in the ALT act as *accumulating factors* – they have both direct effects
on scores AND indirect effects through the lagged dynamics.

Parameter labels follow the unified naming convention:

- `ar_x`, `ar_y`: Autoregressive effects on latent scores. These reflect
  carry-over that compounds with the growth trajectory.

- `cl_xy`, `cl_yx`: Cross-lagged effects on latent scores.

- `I_x`, `I_y`: Latent intercept factors (accumulating).

- `S_x`, `S_y`: Latent slope factors (accumulating).

- `d_var_x`, `d_var_y`: Dynamic residual variances.

- `d_cov_xy`: Dynamic residual covariance.

The ALT conditions on the first observation (wave 1 is exogenous).
Growth factor loadings start from wave 2.

## References

Bollen, K. A., & Curran, P. J. (2006). *Latent curve models: A
structural equation perspective*. Wiley.

Curran, P. J., & Bollen, K. A. (2001). The best of both worlds:
Combining autoregressive and latent curve models. In L. M. Collins & A.
G. Sayer (Eds.), *New methods for the analysis of change* (pp. 107-135).
APA.

Usami, S., Murayama, K., & Hamaker, E. L. (2019). A unified framework of
longitudinal models to examine reciprocal relations. *Psychological
Methods*, 24(5), 637-657.

## Examples

``` r
if (FALSE) { # \dontrun{
syntax <- estimateALT(waves = 5)
cat(syntax)

library(lavaan)
fit <- lavaan::sem(syntax, data = my_data)
summary(fit, fit.measures = TRUE, standardized = TRUE)
} # }
```
