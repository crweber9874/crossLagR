# estimateStTr

Generates model syntax for state trait model

## Usage

``` r
estimateStTr(data = data, waves = 10, ...)
```

## Arguments

- data:

  A data frame containing

- waves:

  The number of waves (time points) in the model.

- time_varying_x:

  A vector of variable names for the X variables in data frame. Ordered
  in time from earliest to latest.

- time_varying_y:

  A vector of variable names for the Y variables in data frame. Ordered
  in time from earliest to latest.

- y_vars:

  A vector of variable names for the Y variables.

## Value

A character string containing the Lavaan model syntax for RI-CLPM. \*
\`model\`: The Lavaan model syntax used for data simulation.

## Details

This function generates the model syntax for a Random Intercept
Cross-Lagged Panel Model (RICLPM) with a specified number of waves. The
model includes latent variables for each wave, stability paths, and
covariances between the latent variables.

## Examples

``` r
# Not Run
# Generate model syntax for a RICLPM with 5 waves
#  model_syntax = estimateRICLPM(data = my_data, x_vars = c("x1", "x2", "x3", "x4", "x5"), y_vars = c("y1", "y2", "y3", "y4", "y5"), waves = 5)
#  cat(model_syntax)
```
