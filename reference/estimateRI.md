# estimateRI

Estimate Random Intercepts model using lmer() from lme4 package

This function fits a random intercepts model using lme4::lmer(), where
individual-specific intercepts are treated as random effects. This
approach is intermediate between OLS (which pools all individuals) and
Fixed Individual Effects (which estimates separate intercepts for each
individual).

## Usage

``` r
estimateRI(
  data,
  y_outcome = "y",
  x_outcome = "x",
  y_lag = "ylag",
  x_lag = "xlag",
  id_var = "id",
  return_models = FALSE,
  control_params = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun =
    20000))
)
```

## Arguments

- data:

  Long-format data frame with columns for variables and lagged
  variables. Must include: id, x, y, xlag, ylag (created by
  reshape_long_sim_cr function).

- y_outcome:

  Character string specifying the Y outcome variable name. Default is
  "y".

- x_outcome:

  Character string specifying the X outcome variable name. Default is
  "x".

- y_lag:

  Character string specifying the lagged Y predictor name. Default is
  "ylag".

- x_lag:

  Character string specifying the lagged X predictor name. Default is
  "xlag".

- id_var:

  Character string specifying the ID variable name. Default is "id".

- return_models:

  Logical. If TRUE, returns the full lmer objects. If FALSE, returns
  just the coefficients. Default is FALSE.

- control_params:

  List of control parameters to pass to lmer(). Default uses lmerControl
  with optimizer = "bobyqa" and optCtrl = list(maxfun = 20000).

## Value

A list containing:

- `xlag_x`: Autoregressive effect of X (X(t-1) -\> X(t))

- `ylag_x`: Cross-lagged effect of Y on X (Y(t-1) -\> X(t))

- `xlag_y`: Cross-lagged effect of X on Y (X(t-1) -\> Y(t))

- `ylag_y`: Autoregressive effect of Y (Y(t-1) -\> Y(t))

- `converged`: Whether both models converged successfully

- `n_obs`: Number of observations used in analysis

- `n_groups`: Number of groups (individuals) in the analysis

- `models`: If return_models=TRUE, includes the fitted lmer objects

- `random_effects`: Random intercept variances and residual variances

## Details

Random Intercepts estimation works by:

1.  Fitting mixed-effects models with random intercepts for each
    individual

2.  Allowing for correlation within individuals while estimating
    population-level effects

3.  Assuming individual intercepts are drawn from a normal distribution

The approach handles unbalanced data well and provides more efficient
estimates than Fixed Individual Effects when the random intercept
assumption is appropriate.

Model specifications:

- Y equation: `y ~ xlag + ylag + (1|id)`

- X equation: `x ~ xlag + ylag + (1|id)`

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate some data
sim_data <- simRICLPM(waves = 4, sample.nobs = 500)$data
long_data <- reshape_long_sim_cr(sim_data)

# Estimate Random Intercepts model
ri_results <- estimateRI(long_data)
print(ri_results)

# Get full model objects
ri_full <- estimateRI(long_data, return_models = TRUE)
summary(ri_full$models$x_model)
summary(ri_full$models$y_model)

# Check random effects
ri_full$random_effects
} # }
```
