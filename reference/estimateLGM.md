# estimateLGM

Generates model syntax for Linear Latent Growth Model (LGM) with options
to constrain parameters and include multiple variables.

## Usage

``` r
estimateLGM(
  waves = 5,
  variable_type = c("univariate", "bivariate"),
  time_scores = NULL,
  constrain_residual_variances = TRUE,
  constrain_residual_covariances = TRUE,
  estimate_quadratic = FALSE,
  center_time = FALSE,
  start_values = FALSE
)
```

## Arguments

- waves:

  The number of waves (time points) in the model.

- variable_type:

  Specify whether estimating "univariate" (single variable) or
  "bivariate" (dual variable). Default is "univariate".

- time_scores:

  A numeric vector specifying the time scores for each wave. If NULL,
  uses 0, 1, 2, ..., waves-1. Default is NULL.

- constrain_residual_variances:

  Logical. If TRUE, constrains residual variances to equality across
  waves. Default is TRUE.

- constrain_residual_covariances:

  Logical. If TRUE, constrains residual covariances to equality across
  waves (bivariate only). Default is TRUE.

- estimate_quadratic:

  Logical. If TRUE, includes quadratic growth factor. Default is FALSE.

- center_time:

  Logical. If TRUE, centers time scores around their mean for
  interpretation. Default is FALSE.

- start_values:

  Logical. If TRUE, provides starting values for key parameters. Default
  is FALSE.

## Value

A character string containing the Lavaan model syntax for LGM.

## Details

This function generates the model syntax for a Linear Latent Growth
Model (LGM) with a specified number of waves. The LGM models systematic
change over time using latent intercept and slope factors.

For univariate models: - Intercept factor (I): Represents initial
level/starting point - Slope factor (S): Represents rate of linear
change over time - Quadratic factor (Q): Represents
acceleration/deceleration (if specified)

For bivariate models: - All components from univariate model for both X
and Y variables - Covariances between growth factors across variables -
Optional constraints on residual parameters

Key parameters estimated: - mean_I_x, mean_I_y: Mean intercepts (initial
levels) - mean_S_x, mean_S_y: Mean slopes (rates of change) - mean_Q_x,
mean_Q_y: Mean quadratic terms (if included) - var_I_x, var_I_y:
Intercept variances (individual differences in initial level) - var_S_x,
var_S_y: Slope variances (individual differences in rate of change) -
var_Q_x, var_Q_y: Quadratic variances (if included) - cov_IS_x,
cov_IS_y: Intercept-slope covariances - var_p, var_q: Latent variable
unique variances - cov_pq: Latent variable covariances (bivariate only)

Time scoring options: - Default: 0, 1, 2, ..., waves-1 (intercept =
initial level) - Custom: User-specified time points - Centered: Time
scores centered around their mean for easier interpretation

The model assumes linear change over time, with optional quadratic
component. Individual differences in both initial levels and rates of
change are estimated.

## Examples

``` r
# Basic univariate LGM with 5 waves
model_syntax <- estimateLGM(waves = 5)
cat(model_syntax)
#> # Latent Variables
#> p1 =~ 1*x1
#> p2 =~ 1*x2
#> p3 =~ 1*x3
#> p4 =~ 1*x4
#> p5 =~ 1*x5
#> 
#> # Latent Growth Factors
#> I =~ 1*p1 + 1*p2 + 1*p3 + 1*p4 + 1*p5
#> S =~ 0*p1 + 1*p2 + 2*p3 + 3*p4 + 4*p5
#> 
#> # Growth Factor Means
#> I ~ mean_I*1
#> S ~ mean_S*1
#> 
#> # Growth Factor Variances and Covariances
#> I ~~ var_I*I
#> S ~~ var_S*S
#> I ~~ cov_IS*S
#> 
#> # Latent Variable Means (Fixed to 0)
#> p1 ~ 0*1
#> p2 ~ 0*1
#> p3 ~ 0*1
#> p4 ~ 0*1
#> p5 ~ 0*1
#> 
#> # Observed Variable Intercepts (Fixed to 0)
#> x1 ~ 0*1
#> x2 ~ 0*1
#> x3 ~ 0*1
#> x4 ~ 0*1
#> x5 ~ 0*1
#> 
#> # Latent Variable Variances
#> p1 ~~ var_p*p1
#> p2 ~~ var_p*p2
#> p3 ~~ var_p*p3
#> p4 ~~ var_p*p4
#> p5 ~~ var_p*p5
#> 
#> # Fix Observed Variable Residuals to Zero
#> x1 ~~ 0*x1
#> x2 ~~ 0*x2
#> x3 ~~ 0*x3
#> x4 ~~ 0*x4
#> x5 ~~ 0*x5

# Bivariate LGM with custom time scores
model_syntax <- estimateLGM(
  waves = 4,
  variable_type = "bivariate",
  time_scores = c(0, 1, 3, 5)
)

# LGM with quadratic growth
model_syntax <- estimateLGM(
  waves = 6,
  estimate_quadratic = TRUE,
  center_time = TRUE
)
```
