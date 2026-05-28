# estimateRICLPM3

Generates model syntax for a three-variable Random Intercept
Cross-Lagged Panel Model (RI-CLPM) with user-specified variable names.

## Usage

``` r
estimateRICLPM3(
  data,
  time_varying_x,
  time_varying_y,
  time_varying_z,
  time_invariant_vars = NULL,
  waves = 10
)
```

## Arguments

- data:

  A data frame containing the user's data.

- time_varying_x:

  A vector of variable names for X, ordered earliest to latest.

- time_varying_y:

  A vector of variable names for Y, ordered earliest to latest.

- time_varying_z:

  A vector of variable names for Z, ordered earliest to latest.

- time_invariant_vars:

  A vector of variable names for time-invariant predictors.

- waves:

  The number of waves (time points) in the model.

## Value

A character string containing lavaan model syntax for the three-variable
RI-CLPM.

## Details

This function generates lavaan syntax for a three-variable RI-CLPM
following the unified framework of Usami, Murayama, & Hamaker (2019).
See
[`estimateRICLPM`](https://crweber9874.github.io/crossLagR/reference/estimateRICLPM.md)
for the unified parameter naming convention.

Note: This variant uses equality constraints across all cross-lagged
paths (`cl`) and all autoregressive paths (`ar`).

## References

Usami, S., Murayama, K., & Hamaker, E. L. (2019). A unified framework of
longitudinal models to examine reciprocal relations. *Psychological
Methods*, 24(5), 637-657.

## Examples

``` r
# Not Run
# model_syntax = estimateRICLPM3(data = my_data,
#   time_varying_x = paste0("x", 1:5),
#   time_varying_y = paste0("y", 1:5),
#   time_varying_z = paste0("z", 1:5), waves = 5)
```
