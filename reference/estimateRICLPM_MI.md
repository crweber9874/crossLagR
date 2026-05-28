# estimateRICLPM_MI

Generates lavaan syntax for a multi-indicator Random Intercept
Cross-Lagged Panel Model (RI-CLPM).

## Usage

``` r
estimateRICLPM_MI(
  waves = 5,
  n_indicators_x = 3,
  n_indicators_y = 3,
  invariance = c("scalar", "metric", "configural", "strict"),
  correlated_residuals = TRUE,
  constrain_beta = TRUE,
  constrain_omega = TRUE,
  constrain_residual_variances = TRUE,
  constrain_residual_covariances = TRUE,
  start_values = FALSE
)
```

## Arguments

- waves:

  Integer. Number of waves (\>= 2).

- n_indicators_x:

  Integer. Indicators for X per wave (\>= 1).

- n_indicators_y:

  Integer. Indicators for Y per wave (\>= 1).

- invariance:

  Character. One of `"configural"`, `"metric"`, `"scalar"`, `"strict"`.
  Default `"scalar"`.

- correlated_residuals:

  Logical. Correlate same-item residuals across waves. Default `TRUE`.

- constrain_beta, constrain_omega, constrain_residual_variances,
  constrain_residual_covariances:

  Logical. Equality constraints across waves.

- start_values:

  Logical.

## Value

Character string of lavaan syntax.

## Details

Second-order structure: wave-specific latent factor `eta_x_w` aggregates
indicators `x{w}_1 ... x{w}_{J}`; it is then perfectly decomposed into a
between-person trait `I_x` and a within-person deviation `p_w`
(analogously for Y). All AR/CL dynamics operate on `p_w`, `q_w`.
Same-item residuals are correlated across waves (per Mulder & Hamaker,
2021).

## References

Hamaker, E. L., Kuiper, R. M., & Grasman, R. P. P. P. (2015). A critique
of the cross-lagged panel model. *Psychological Methods*, 20, 102-116.

Mulder, J. D., & Hamaker, E. L. (2021). Three extensions of the random
intercept cross-lagged panel model. *SEM*, 28(4), 638-648.

Usami, S., Murayama, K., & Hamaker, E. L. (2019). *Psychological
Methods*, 24(5), 637-657.

## Examples

``` r
if (FALSE) { # \dontrun{
syntax <- estimateRICLPM_MI(waves = 4, n_indicators_x = 3, n_indicators_y = 3)
fit <- lavaan::lavaan(syntax, data = my_data, meanstructure = TRUE)
} # }
```
