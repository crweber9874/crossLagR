# plotWithinBetween

Violin plot comparing observed between/within variance components to a
parametric bootstrap.

## Usage

``` r
plotWithinBetween(observed, simulated)
```

## Arguments

- observed:

  A data frame returned by
  [`withinBetween`](https://crweber9874.github.io/crossLagR/reference/withinBetween.md).

- simulated:

  A data frame returned by
  [`simWithinBetween`](https://crweber9874.github.io/crossLagR/reference/simWithinBetween.md).

## Value

A ggplot object with one facet per variable, two violins (between,
within) showing the bootstrap distribution, and a point overlay marking
the observed value.
