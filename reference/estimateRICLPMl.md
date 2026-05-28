# estimateRICLPMl

Generates lavaan model syntax for RI-CLPM with multiple indicators per
construct, where the first indicator has a fixed factor loading of 1 and
subsequent indicators have free factor loadings.

## Usage

``` r
estimateRICLPMl(
  time_varying_x = list(c("x1_1", "x1_2"), c("x2_1", "x2_2")),
  time_varying_y = list(c("y1_1", "y1_2"), c("y2_1", "y2_2")),
  time_varying_z = NULL,
  time_invariant_vars = NULL,
  waves
)
```

## Arguments

- time_varying_x:

  A list of vectors, each containing variable names for the X indicators
  at each wave.

- time_varying_y:

  A list of vectors, each containing variable names for the Y indicators
  at each wave.

- time_varying_z:

  An optional list of vectors for Z indicators at each wave.

- time_invariant_vars:

  A vector of variable names for time-invariant variables.

- waves:

  The number of waves (time points) in the model.

## Value

A character string containing lavaan model syntax.

## Details

This function generates lavaan syntax for an RI-CLPM with multiple
indicators, following the unified framework of Usami, Murayama, &
Hamaker (2019). See
[`estimateRICLPM`](https://crweber9874.github.io/crossLagR/reference/estimateRICLPM.md)
for the unified parameter naming convention.

Parameter labels:

- `ar_x`, `ar_y`, `ar_z`: Autoregressive (within-person carry-over).

- `cl_xy`, `cl_yx`, etc.: Cross-lagged effects.

- `I_x`, `I_y`, `I_z`: Trait factors (random intercepts).

## References

Usami, S., Murayama, K., & Hamaker, E. L. (2019). A unified framework of
longitudinal models to examine reciprocal relations. *Psychological
Methods*, 24(5), 637-657.

## Examples

``` r
if (FALSE) { # \dontrun{
model_syntax <- estimateRICLPMl(
  time_varying_x = list(c("x1_1", "x1_2"), c("x2_1", "x2_2")),
  time_varying_y = list(c("y1_1", "y1_2"), c("y2_1", "y2_2")),
  waves = 2)
cat(model_syntax)
} # }
```
