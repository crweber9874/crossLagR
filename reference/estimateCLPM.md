# estimateCLPM. This follows the Usami, Murayama, & Hamaker (2019) unified framework for longitudinal models.

Generates lavaan model syntax for the Cross-Lagged Panel Model (CLPM).

## Usage

``` r
estimateCLPM(
  waves = 5,
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

  The number of waves (time points) in the model.

- constrain_beta:

  Logical. If TRUE, constrains autoregressive effects to equality across
  waves. Default is TRUE.

- constrain_omega:

  Logical. If TRUE, constrains cross-lagged effects to equality across
  waves. Default is TRUE.

- constrain_residual_variances:

  Logical. If TRUE, constrains residual variances to equality across
  waves. Default is TRUE.

- constrain_residual_covariances:

  Logical. If TRUE, constrains residual covariances to equality across
  waves. Default is TRUE.

- estimate_means:

  Logical. If TRUE, estimates means for the first wave. Default is TRUE.

- start_values:

  Logical. If TRUE, provides starting values for key parameters. Default
  is FALSE.

## Value

A character string containing the lavaan model syntax for the CLPM.

## Details

This function generates lavaan syntax for a Cross-Lagged Panel Model
following the unified framework of Usami, Murayama, & Hamaker (2019). In
this framework, the CLPM uses only the decomposition equation
(group-mean centering) and the dynamic equation (lagged regression),
with no unique factors, trait factors, or growth factors.

Parameter labels follow the unified naming convention:

- `ar_y`, `ar_x`: Autoregressive effects. In the CLPM these reflect
  rank-order stability, conflating between- and within-person processes.

- `cl_xy`, `cl_yx`: Cross-lagged effects (Granger-causal paths). `cl_xy`
  = effect of X on Y; `cl_yx` = effect of Y on X.

- `d_var_x`, `d_var_y`: Dynamic residual (innovation) variances.

- `d_cov_xy`: Dynamic residual covariance (within-time association).

## References

Usami, S., Murayama, K., & Hamaker, E. L. (2019). A unified framework of
longitudinal models to examine reciprocal relations. *Psychological
Methods*, 24(5), 637-657.

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate CLPM syntax with 5 waves
syntax <- estimateCLPM(waves = 5)
cat(syntax)

# Fit with lavaan
library(lavaan)
fit <- lavaan(syntax, data = my_data, meanstructure = TRUE)
summary(fit, fit.measures = TRUE)
} # }
```
