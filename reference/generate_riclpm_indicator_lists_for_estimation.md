# generate_riclpm_indicator_lists_for_estimation

Generates lists of indicator variable names in the format required by
the estimateRICLPMl function.

## Usage

``` r
generate_riclpm_indicator_lists_for_estimation(
  waves,
  num_indicators,
  has_z = FALSE,
  x_prefix = "x",
  y_prefix = "y",
  z_prefix = "z",
  indicator_separator = "_"
)
```

## Arguments

- waves:

  The number of waves (time points) in the model.

- num_indicators:

  The number of indicators per construct (X, Y, and optionally Z) at
  each wave.

- has_z:

  Logical, indicating whether the model includes a Z variable. Defaults
  to FALSE.

- x_prefix:

  The prefix for the X variable indicators. Defaults to "x".

- y_prefix:

  The prefix for the Y variable indicators. Defaults to "y".

- z_prefix:

  The prefix for the Z variable indicators. Defaults to "z".

- indicator_separator:

  The separator between the prefix, indicator number, and wave number.
  Defaults to "\_".

## Value

A list containing the \`time_varying_x\`, \`time_varying_y\`, and
optionally \`time_varying_z\` lists, formatted for \`estimateRICLPMl\`.

## Examples

``` r
indicator_lists_2waves_3indicators_for_est <- generate_riclpm_indicator_lists_for_estimation(waves = 2, num_indicators = 3)
print(indicator_lists_2waves_3indicators_for_est)
#> $time_varying_x
#> $time_varying_x[[1]]
#> [1] "x1_1" "x1_2"
#> 
#> $time_varying_x[[2]]
#> [1] "x2_1" "x2_2"
#> 
#> $time_varying_x[[3]]
#> [1] "x3_1" "x3_2"
#> 
#> 
#> $time_varying_y
#> $time_varying_y[[1]]
#> [1] "y1_1" "y1_2"
#> 
#> $time_varying_y[[2]]
#> [1] "y2_1" "y2_2"
#> 
#> $time_varying_y[[3]]
#> [1] "y3_1" "y3_2"
#> 
#> 

indicator_lists_3waves_2indicators_withZ_for_est <- generate_riclpm_indicator_lists_for_estimation(waves = 3, num_indicators = 2, has_z = TRUE)
print(indicator_lists_3waves_2indicators_withZ_for_est)
#> $time_varying_x
#> $time_varying_x[[1]]
#> [1] "x1_1" "x1_2" "x1_3"
#> 
#> $time_varying_x[[2]]
#> [1] "x2_1" "x2_2" "x2_3"
#> 
#> 
#> $time_varying_y
#> $time_varying_y[[1]]
#> [1] "y1_1" "y1_2" "y1_3"
#> 
#> $time_varying_y[[2]]
#> [1] "y2_1" "y2_2" "y2_3"
#> 
#> 
#> $time_varying_z
#> $time_varying_z[[1]]
#> [1] "z1_1" "z1_2" "z1_3"
#> 
#> $time_varying_z[[2]]
#> [1] "z2_1" "z2_2" "z2_3"
#> 
#> 

indicator_lists_2waves_2indicators_custom_prefix_for_est <- generate_riclpm_indicator_lists_for_estimation(waves = 2, num_indicators = 2, x_prefix = "varX", y_prefix = "varY")
print(indicator_lists_2waves_2indicators_custom_prefix_for_est)
#> $time_varying_x
#> $time_varying_x[[1]]
#> [1] "varX1_1" "varX1_2"
#> 
#> $time_varying_x[[2]]
#> [1] "varX2_1" "varX2_2"
#> 
#> 
#> $time_varying_y
#> $time_varying_y[[1]]
#> [1] "varY1_1" "varY1_2"
#> 
#> $time_varying_y[[2]]
#> [1] "varY2_1" "varY2_2"
#> 
#> 
```
