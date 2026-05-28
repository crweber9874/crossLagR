# compCTSEM

Create a continuous time model specification using the ctsem package.

Wraps `ctsem::ctModel` to define a continuous-time structural equation
model. Returns the model object for subsequent fitting.

## Usage

``` r
compCTSEM(
  data = dat_long,
  id = "id",
  manifest_var_names = c("x", "y"),
  latent_var_names = c("x", "y"),
  LAMBDA = diag(2),
  type = "stanct",
  time_points = 10,
  cores = cores,
  chains = chains,
  ...
)
```

## Arguments

- data:

  A long-format data frame.

- id:

  Character. Name of the ID column. Default is `"id"`.

- manifest_var_names:

  Character vector of manifest variable names. Default is `c("x", "y")`.

- latent_var_names:

  Character vector of latent variable names. Default is `c("x", "y")`.

- LAMBDA:

  Factor loading matrix. Default is `diag(2)`.

- type:

  Character. Type of ctsem model. Default is `"stanct"`.

- time_points:

  Integer. Number of time points. Default is 10.

- cores:

  Integer. Number of CPU cores for parallel computation.

- chains:

  Integer. Number of MCMC chains.

- ...:

  Additional arguments passed to `ctsem::ctModel`.

## Value

A list; returns model and drift matrix \#'

## Examples

``` r
if (FALSE) { # \dontrun{
model <- compCTSEM(data = my_long_data, time_points = 5)
} # }
```
