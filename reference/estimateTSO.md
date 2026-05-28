# estimateTSO

Generates lavaan model syntax for a bivariate Trait-State-Occasion (TSO)
variance decomposition.

## Usage

``` r
estimateTSO(
  waves,
  constrain_state_variances = TRUE,
  constrain_state_covariances = TRUE
)
```

## Arguments

- waves:

  Integer. Number of waves (time points). Must be \>= 3.

- constrain_state_variances:

  Logical. If TRUE, constrains state (within-person) variances to
  equality across waves. Default is TRUE.

- constrain_state_covariances:

  Logical. If TRUE, constrains within-time state covariances to equality
  across waves. Default is TRUE.

## Value

A character string containing lavaan model syntax.

## Details

This is a simple variance-decomposition model for the bivariate panel
case with a single indicator per wave per variable. Each observed
variable is decomposed into a time-invariant trait component and a
time-varying state (occasion) component: \$\$x\_{it} = T\_{xi} +
s\_{xit}\$\$ \$\$y\_{it} = T\_{yi} + s\_{yit}\$\$

The trait factors `T_x` and `T_y` capture between-person variance, while
the state residuals capture within-person variance. The variance ratio
`var(T_x) / var(s_x)` estimates the between-to-within variance ratio.

This model is the random-intercept-only special case of the RI-CLPM (no
lagged dynamics), which is what is needed for variance decomposition.
With a single indicator per occasion, the classic multi-indicator
Cole-Martin-Steiger TSO is not identified, so this simplified version is
used.

Parameter labels:

- `T_var_x`, `T_var_y`: Trait variances (between-person).

- `T_cov_xy`: Trait covariance.

- `s_var_x`, `s_var_y`: State variances (within-person).

- `s_cov_xy`: Within-time state covariance.
