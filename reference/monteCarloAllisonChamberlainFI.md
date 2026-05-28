# monteCarloAllisonChamberlainFI

Monte Carlo simulation for the Allison-Chamberlain Fixed-Effects SEM
estimator.

This function combines data simulation and estimation using the
Allison-Chamberlain fixed-effects SEM approach over a fixed number of
trials. It generates data from either a CLPM, RI-CLPM, or CLPM with
unmeasured confounder DGP, fits the Allison-Chamberlain model, and
extracts cross-lagged and autoregressive coefficients.

## Usage

``` r
monteCarloAllisonChamberlainFI(
  trials = 10,
  waves = 5,
  stability_p = 0.25,
  stability_q = 0.25,
  cross_p = 0,
  cross_q = 0,
  variance_p = 1,
  variance_q = 1,
  sample_size = 2500,
  dgp = "clpm",
  variance_between_x = 0.5,
  variance_between_y = 0.5,
  confounder_p = 0.3,
  confounder_q = 0.3,
  confounder_variance = 1,
  confounder_stability = 0.4,
  include_confounder = TRUE,
  confounder_type = "time_variant",
  cov_pq = 0.1,
  verbose = FALSE,
  ...
)
```

## Arguments

- trials:

  Integer. Number of Monte Carlo trials. Default is 10.

- waves:

  Integer. Number of waves (time points). Must be \>= 3. Default is 5.

- stability_p:

  Numeric. Autoregressive effect for x variable. Default is 0.25.

- stability_q:

  Numeric. Autoregressive effect for y variable. Default is 0.25.

- cross_p:

  Numeric. Cross-lagged effect of y on x. Default is 0.

- cross_q:

  Numeric. Cross-lagged effect of x on y. Default is 0.

- variance_p:

  Numeric. Within-person variance for x. Default is 1.

- variance_q:

  Numeric. Within-person variance for y. Default is 1.

- sample_size:

  Integer. Sample size for each simulation. Default is 2500.

- dgp:

  Character. Data generation process: `"clpm"`, `"riclpm"`, or
  `"clpmu"`. Default is `"clpm"`.

- variance_between_x:

  Numeric. Between-person variance for x (RI-CLPM DGP). Default is 0.5.

- variance_between_y:

  Numeric. Between-person variance for y (RI-CLPM DGP). Default is 0.5.

- confounder_p:

  Numeric. Effect of confounder on x (CLPMU DGP). Default is 0.3.

- confounder_q:

  Numeric. Effect of confounder on y (CLPMU DGP). Default is 0.3.

- confounder_variance:

  Numeric. Variance of confounder (CLPMU DGP). Default is 1.

- confounder_stability:

  Numeric. AR effect of confounder (CLPMU DGP). Default is 0.4.

- include_confounder:

  Logical. Include confounder in CLPMU model. Default is TRUE.

- confounder_type:

  Character. `"time_variant"` or `"time_invariant"`. Default is
  `"time_variant"`.

- cov_pq:

  Numeric. Covariance between p and q innovations. Default is 0.1.

- verbose:

  Logical. Print progress messages. Default is FALSE.

- ...:

  Additional arguments (currently unused).

## Value

A data frame containing estimated coefficients (`xlag_x`, `ylag_x`,
`xlag_y`, `ylag_y`), convergence status, and error/warning information
for each trial and parameter combination.

## Details

The function follows the same structure as other Monte Carlo functions
in the package (`monteCarloCLPM`, `monteCarloRICLPM`). It generates data
using the specified DGP, fits the Allison-Chamberlain model via
[`lavaan::sem()`](https://rdrr.io/pkg/lavaan/man/sem.html), and extracts
the labeled regression coefficients.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic Monte Carlo with CLPM DGP
results <- monteCarloAllisonChamberlainFI(trials = 5, waves = 4, sample_size = 500)
colMeans(results[, c("xlag_x", "ylag_x", "xlag_y", "ylag_y")], na.rm = TRUE)

# Test with RI-CLPM DGP (data has fixed effects, estimator should capture them)
results_ri <- monteCarloAllisonChamberlainFI(
    trials = 10, waves = 5, dgp = "riclpm",
    variance_between_x = 1, variance_between_y = 1
)
} # }
```
