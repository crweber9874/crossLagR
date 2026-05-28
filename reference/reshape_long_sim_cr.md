# reshape_long_sim_cr

Reshape Data from Wide to Long Format for Cross-Lagged Analysis

This function takes a wide-format data frame containing variables named
'x1', 'x2', ..., 'xn', and 'y1', 'y2', ..., 'yn', and reshapes it into a
long format suitable for cross-lagged analysis. The resulting data frame
includes columns for individual ID, wave, x, y, xlag (lagged x), and
ylag (lagged y).

## Usage

``` r
reshape_long_sim_cr(data = dat)
```

## Arguments

- data:

  A wide-format data frame containing variables named 'x1', 'x2', ...,
  'xn', and 'y1', 'y2', ..., 'yn'.

## Value

A long-format data frame with columns 'id', 'wave', 'x', 'y', 'xlag',
and 'ylag'.

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- simCLPM(waves = 3)
head(dat)
long_dat <- reshape_long_sim_cr(dat)
head(long_dat)
} # }
```
