# estimateRICLPM

Generates lavaan model syntax for the Random Intercept Cross-Lagged
Panel Model (RI-CLPM).

## Usage

``` r
estimateRICLPM(
  waves,
  include_z = FALSE,
  time_invariant_vars = NULL,
  constrain_beta = TRUE,
  constrain_omega = TRUE,
  constrain_residual_variances = TRUE,
  constrain_residual_covariances = TRUE,
  start_values = FALSE,
  label_autoregressive = c("ar_x", "ar_y", "ar_z"),
  label_crosslagged = c("cl_yx", "cl_xy", "cl_xz", "cl_yz", "cl_zx", "cl_zy")
)
```

## Arguments

- waves:

  The number of waves (time points) in the model.

- include_z:

  Logical. If TRUE, includes a third variable Z in the model. Default is
  FALSE.

- time_invariant_vars:

  A vector of variable names for time-invariant variables.

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

- start_values:

  Logical. If TRUE, provides starting values for key parameters. Default
  is FALSE.

- label_autoregressive:

  Character vector for autoregressive parameter labels. Default is
  `c("ar_x", "ar_y", "ar_z")`.

- label_crosslagged:

  Character vector for cross-lagged parameter labels. Default is
  `c("cl_yx", "cl_xy", "cl_xz", "cl_yz", "cl_zx", "cl_zy")`.

## Value

A character string containing the lavaan model syntax for the RI-CLPM.

## Details

This function generates lavaan syntax for a Random Intercept
Cross-Lagged Panel Model following the unified framework of Usami,
Murayama, & Hamaker (2019). The RI-CLPM uses both the decomposition
equation (separating stable trait factors from within-person deviations)
and the dynamic equation (lagged regression on the within-person
deviations).

Parameter labels follow the unified naming convention:

- `ar_x`, `ar_y`, `ar_z`: Autoregressive effects. In the RI-CLPM these
  represent *within-person carry-over*, not rank-order stability (which
  is captured by the trait factor I).

- `cl_xy`, `cl_yx`: Cross-lagged effects operating at the within-person
  level. `cl_xy` = effect of X on Y; `cl_yx` = effect of Y on X.

- `d_var_x`, `d_var_y`, `d_var_z`: Dynamic residual (innovation)
  variances for within-person deviations.

- `d_cov_xy`, `d_cov_xz`, `d_cov_yz`: Dynamic residual covariances
  (within-time, within-person associations).

- `I_x`, `I_y`, `I_z`: Stable trait factors (random intercepts)
  representing time-invariant between-person differences. These are
  unique to models that include decomposition of trait and state (no
  direct analogue in the CLPM).

Within-person latent deviations are labeled `p` (for X), `q` (for Y),
and `r` (for Z). These structural latent variables have no direct
analogue across all models and retain model-specific names.

## References

Usami, S., Murayama, K., & Hamaker, E. L. (2019). A unified framework of
longitudinal models to examine reciprocal relations. *Psychological
Methods*, 24(5), 637-657.

Hamaker, E. L., Kuiper, R. M., & Grasman, R. P. P. P. (2015). A critique
of the cross-lagged panel model. *Psychological Methods*, 20, 102-116.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic RI-CLPM with 5 waves (assumes data has x1-x5, y1-y5)
model_syntax <- estimateRICLPM(waves = 5)

# RI-CLPM with unconstrained cross-lagged effects
model_syntax <- estimateRICLPM(waves = 5, constrain_omega = FALSE)

# RI-CLPM with third variable (assumes data has x1-x5, y1-y5, z1-z5)
model_syntax <- estimateRICLPM(waves = 5, include_z = TRUE)

# Fit the model
library(lavaan)
fitted_model <- lavaan(model_syntax, data = your_data)
summary(fitted_model, fit.measures = TRUE, standardized = TRUE)
} # }
```
