# monteCarloLChange

Monte Carlo simulation for Latent Change Score Models This function
combines the estimation and simulation of latent change models over a
fixed number of trials. The output is a data frame that includes the
estimated coefficients across these trials

## Usage

``` r
monteCarloLChange(
  trials = 10,
  waves = 5,
  dgp = "lchange",
  model_type = "dual_change",
  stability_q = 0.25,
  stability_p = 0.25,
  variance_between_y = 0.5,
  variance_between_x = 0.5,
  cross_p = 0,
  cross_q = 0,
  variance_p = 1,
  variance_q = 1,
  sample_size = 2500,
  cov_pq = 0.1,
  confounder_p = 0.3,
  confounder_q = 0.3,
  confounder_variance = 1,
  confounder_stability = 0.4,
  include_confounder = TRUE,
  confounder_type = "time_variant",
  ar_x = -0.15,
  ar_y = -0.2,
  cl_x = -0.1,
  cl_y = -0.1,
  change_x = 0.08,
  change_y = 0.08,
  phi_x = 0.15,
  phi_y = 0.15,
  initial_var_x = 1,
  initial_var_y = 1,
  constant_change_var_x = 0.5,
  constant_change_var_y = 0.5,
  residual_variance_x = 0.5,
  residual_variance_y = 0.5,
  verbose = FALSE,
  ...
)
```

## Arguments

- trials:

  The number of trials for the Monte Carlo simulation.

- waves:

  The number of waves (time points) in the model. Must specify three or
  more waves for identification.

- dgp:

  Specify the data generation method: "riclpm", "clpm", "clpmu",
  "lchange". Are data generated under one of these regimes?

- model_type:

  Type of latent change model: "latent_change" or "dual_change". Use
  dual change with multiple variables

- stability_q:

  Autoregressive effect for q (y) variable. These parameters are only
  relevant when specifying an RICLPM or CLPM. They are \*\*not\*\* part
  of the dual change model

- stability_p:

  Autoregressive effect for p (x) variable. These parameters are only
  relevant when specifying an RICLPM or CLPM. They are \*\*not\*\* part
  of the dual change model

- variance_between_y:

  Between-person variance for y variable. Ony relevant for the RICLPM.

- variance_between_x:

  Between-person variance for x variable. Only relevant for RICLPM.

- cross_p:

  Cross-lagged effect of y on x. Only relevant to the CLPM and RICLPM
  specification.

- cross_q:

  Cross-lagged effect of x on y. Only relevant to the CLPM and RICLPM
  specification.

- variance_p:

  Within-person variance for p variable. Only relevant to the RICLPM and
  CLPM specification.

- variance_q:

  Within-person variance for q variable. Only relevant ot the RICLPM and
  CLPM specification

- sample_size:

  Sample size for simulation

- cov_pq:

  Covariance between p and q (for CLPM/CLPMU)

- confounder_p:

  Effect of confounder on x variables (for clpmu)

- confounder_q:

  Effect of confounder on y variables (for clpmu)

- confounder_variance:

  Variance of the confounder (for clpmu)

- confounder_stability:

  Autoregressive effect of confounder (for clpmu)

- include_confounder:

  Whether to include confounder in clpmu model

- confounder_type:

  Character string specifying confounder type for clpmu dgp. Must be one
  of "time_variant" (default) or "time_invariant".

- ar_x:

  Proportional effect for X in latent change DGP

- ar_y:

  Proportional effect for Y in latent change DGP

- cl_x:

  Coupling effect in latent change DGP

- cl_y:

  Coupling effect in latent change DGP

- change_x:

  Change-to-change effect in latent change DGP

- change_y:

  Change-to-change effect in latent change DGP

- phi_x:

  Change autoregression in latent change DGP

- phi_y:

  Change autoregression in latent change DGP

- initial_var_x:

  Initial variance for X in latent change DGP

- initial_var_y:

  Initial variance for Y in latent change DGP

- constant_change_var_x:

  Constant change variance for X in latent change DGP

- constant_change_var_y:

  Constant change variance for Y in latent change DGP

- residual_variance_x:

  Measurement error for X in latent change DGP

- residual_variance_y:

  Measurement error for Y in latent change DGP

- verbose:

  Whether to print progress and error messages.

## Value

A data frame containing the results of the simulation.
