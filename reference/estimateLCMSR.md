# estimateLCMSR. This is just the GCLPM. Structure the means here.

Generates lavaan model syntax for the bivariate Latent Curve Model with
Structured Residuals (LCM-SR), also known as the General Cross-Lagged
Model (GCLM).

## Usage

``` r
estimateLCMSR(
  waves = 5,
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

- estimate_quadratic:

  Logical. If TRUE, includes quadratic growth factors. Requires \>= 4
  waves. Default is FALSE.

## Value

A character string containing lavaan model syntax.

## Details

The LCM-SR (also called the GCLM or RI-CLPM with structured means)
combines latent growth curves with autoregressive and cross-lagged
dynamics in the *residuals* from those curves, following the unified
framework of Usami, Murayama, & Hamaker (2019).

Unlike the ALT model (where AR/CL dynamics operate on the scores
directly and growth factors are accumulating), the LCM-SR first removes
the growth trajectory, then models the deviations with AR/CL dynamics.
This means the growth factors (I, S) have *only direct effects* on
scores – they are true growth factors, not accumulating factors.

Parameter labels follow the unified naming convention:

- `ar_x`, `ar_y`: Autoregressive effects on within-person deviations
  from the growth trajectory.

- `cl_xy`, `cl_yx`: Cross-lagged effects on within-person deviations.

- `I_x`, `I_y`: Growth intercept factors (direct effects only).

- `S_x`, `S_y`: Growth slope factors (direct effects only).

- `d_var_x`, `d_var_y`: Dynamic residual variances.

- `d_cov_xy`: Dynamic residual covariance.

Within-person deviations from the growth curve are labeled `p` (for X)
and `q` (for Y), analogous to the RI-CLPM's within-person factors.

## References

Curran, P. J., Howard, A. L., Bainter, S. A., Lane, S. T., & McGinley,
J. S. (2014). The separation of between-person and within-person
components of individual change over time. *Annual Review of
Psychology*, 65, 637-660.

Berry, D., & Willoughby, M. T. (2017). On the practical interpretability
of cross-lagged panel models. *Structural Equation Modeling*, 24(3),
455-468.

Usami, S., Murayama, K., & Hamaker, E. L. (2019). A unified framework of
longitudinal models to examine reciprocal relations. *Psychological
Methods*, 24(5), 637-657.

## Examples

``` r
if (FALSE) { # \dontrun{
syntax <- estimateLCMSR(waves = 5)
cat(syntax)

library(lavaan)
fit <- lavaan::sem(syntax, data = my_data, meanstructure = TRUE)
summary(fit, fit.measures = TRUE, standardized = TRUE)
} # }
```
