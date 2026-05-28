# withinBetween

Decompose panel variables into between-person and within-person variance
components.

## Usage

``` r
withinBetween(data, id, time, vars, method = c("dplyr", "tso"))
```

## Arguments

- data:

  A long-format panel data frame.

- id:

  Character. Name of the unit (person) identifier column.

- time:

  Character. Name of the time/wave column.

- vars:

  Character vector of variable names to decompose.

- method:

  Either `"dplyr"` (default) for an empirical decomposition via
  group-mean centering, or `"tso"` to fit a bivariate
  Trait-State-Occasion model with lavaan and read the structural
  variance estimates. The `"tso"` option requires exactly two variables.

## Value

A data frame with one row per variable and columns `variable`,
`between_var`, `within_var`, `total_var`, `icc`, `method`.

## Details

The `"dplyr"` method estimates between-person variance as the variance
of person-specific means and within-person variance as the variance of
deviations from those means. The `"tso"` method fits the bivariate
Trait-State-Occasion model from
[`estimateTSO`](https://crweber9874.github.io/crossLagR/reference/estimateTSO.md)
and reads `T_var_x`, `T_var_y` (between) and `s_var_x`, `s_var_y`
(within) from the lavaan parameter table.

## Examples

``` r
if (FALSE) { # \dontrun{
withinBetween(upenn_long, id = "identifier", time = "wave",
              vars = c("party_identification", "political_legitimacy"))

withinBetween(upenn_long, id = "identifier", time = "wave",
              vars = c("party_identification", "political_legitimacy"),
              method = "tso")
} # }
```
