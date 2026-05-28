# simWithinBetween

Parametric bootstrap of between- and within-person variance components
implied by a
[`withinBetween`](https://crweber9874.github.io/crossLagR/reference/withinBetween.md)
fit.

## Usage

``` r
simWithinBetween(wb, n_units, n_waves, n_sims = 500, seed = NULL)
```

## Arguments

- wb:

  A data frame returned by
  [`withinBetween()`](https://crweber9874.github.io/crossLagR/reference/withinBetween.md).

- n_units:

  Integer. Number of units (persons) per simulated panel. Defaults to a
  sensible value matching typical panels.

- n_waves:

  Integer. Number of waves per simulated panel.

- n_sims:

  Integer. Number of bootstrap replicates. Default 500.

- seed:

  Optional integer seed for reproducibility.

## Value

A data frame with columns `variable`, `sim`, `between_var`,
`within_var`.

## Details

For each variable, simulates `n_sims` panels from \\y\_{it} = \alpha_i +
\varepsilon\_{it}\\ with \\\alpha_i \sim N(0, \sigma^2_B)\\ and
\\\varepsilon\_{it} \sim N(0, \sigma^2_W)\\, then recomputes the
empirical between- and within-person variance components in each
replicate. The resulting distribution describes the sampling variability
of the decomposition under the implied model.
