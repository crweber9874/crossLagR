# estimateRICLPM_nolag

Generates lavaan model syntax for RI-CLPM without autoregressive paths
(cross-lagged only).

## Usage

``` r
estimateRICLPM_nolag(
  time_varying_x,
  time_varying_y,
  time_varying_z = NULL,
  time_invariant_vars = NULL,
  waves
)
```

## Arguments

- time_varying_x:

  A vector of variable names for X, ordered earliest to latest.

- time_varying_y:

  A vector of variable names for Y, ordered earliest to latest.

- time_varying_z:

  An optional vector of variable names for Z.

- time_invariant_vars:

  A vector of variable names for time-invariant variables.

- waves:

  The number of waves (time points) in the model.

## Value

A character string containing lavaan model syntax.

## Details

This function generates lavaan syntax for an RI-CLPM with no
autoregressive paths (cross-lagged only), following the unified
framework of Usami, Murayama, & Hamaker (2019). See
[`estimateRICLPM`](https://crweber9874.github.io/crossLagR/reference/estimateRICLPM.md)
for the unified parameter naming convention.

Parameter labels:

- `cl_yx`: Cross-lagged effect of Y on X (no AR for X/Y).

- `cl_xy`: Cross-lagged effect of X on Y.

- `I_x`, `I_y`, `I_z`: Trait factors (random intercepts).

## References

Usami, S., Murayama, K., & Hamaker, E. L. (2019). A unified framework of
longitudinal models to examine reciprocal relations. *Psychological
Methods*, 24(5), 637-657.

## Examples

``` r
if (FALSE) { # \dontrun{
syntax <- estimateRICLPM_nolag(
  time_varying_x = paste0("x", 1:5),
  time_varying_y = paste0("y", 1:5),
  waves = 5
)
cat(syntax)
} # }
```
