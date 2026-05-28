# riclpm_graph

This function draws a tidySEM graph of a RI-CLPM model

Create a simulated dataset of a cross-lagged model with a specified
number of waves and structural parameters.

## Usage

``` r
riclpm_graph(model = fit, waves = c(3, 4, 5, 6))
```

## Arguments

- waves:

  Number of waves in the fitted model.

- fit:

  A fitted lavaan model object.

## Value

RI-CLPM SEM plot
