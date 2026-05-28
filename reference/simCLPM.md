# simCLPM

Simulate data from a Cross-Lagged Panel Model (CLPM) with options for
parameter constraints and starting values.

## Usage

``` r
simCLPM(
  waves = 5,
  beta_y = 0.4,
  beta_x = 0.4,
  omega_xy = 0.2,
  omega_yx = 0.2,
  var_y = 1,
  var_x = 1,
  cov_xy = 0.3,
  var_y1 = NULL,
  var_x1 = NULL,
  cov_xy1 = NULL,
  mean_y1 = 0,
  mean_x1 = 0,
  sample_size = 1000,
  seed = NULL,
  ...
)
```

## Arguments

- waves:

  The number of waves (time points) in the model.

- beta_y:

  The autoregressive effect for the y variable (stability parameter).

- beta_x:

  The autoregressive effect for the x variable (stability parameter).

- omega_xy:

  The cross-lagged effect of x on y at the next time point.

- omega_yx:

  The cross-lagged effect of y on x at the next time point.

- var_y:

  The residual variance for the y latent variables.

- var_x:

  The residual variance for the x latent variables.

- cov_xy:

  The covariance between x and y within the same time point.

- var_y1:

  The variance for the first wave y variable. If NULL, uses var_y.

- var_x1:

  The variance for the first wave x variable. If NULL, uses var_x.

- cov_xy1:

  The covariance for the first wave between x and y. If NULL, uses
  cov_xy.

- mean_y1:

  The mean for the first wave y variable. Default is 0.

- mean_x1:

  The mean for the first wave x variable. Default is 0.

- sample_size:

  The number of observations to simulate. Default is 1000.

- seed:

  Random seed for reproducibility. Default is NULL.

- ...:

  Additional arguments to pass to the \`lavaan::simulateData\` function.

## Value

A list containing two elements: \* \`model\`: The Lavaan model syntax
used for data simulation. \* \`data\`: The simulated data in a data
frame format.

## Details

This function simulates data from a Cross-Lagged Panel Model (CLPM)
where: - p variables represent latent y constructs - q variables
represent latent x constructs - beta parameters represent autoregressive
(stability) effects - omega parameters represent cross-lagged
(reciprocal) effects

The CLPM assumes all relationships are due to within-person processes
and does not separate stable between-person differences from dynamic
effects (unlike the RI-CLPM).

Parameter interpretation: - beta_y: How much y at time t predicts y at
time t+1 - beta_x: How much x at time t predicts x at time t+1 -
omega_xy: How much x at time t predicts y at time t+1 - omega_yx: How
much y at time t predicts x at time t+1

## Examples

``` r
# Basic CLPM simulation
sim_result <- simCLPM(
  waves = 5,
  beta_y = 0.4,
  beta_x = 0.3,
  omega_xy = 0.2,
  omega_yx = 0.1,
  sample_size = 500
)

# View the generated model
cat(sim_result$model)
#>     p1 =~ 1*y1
#>     q1 =~ 1*x1
#>     p2 =~ 1*y2
#>     q2 =~ 1*x2
#>     p3 =~ 1*y3
#>     q3 =~ 1*x3
#>     p4 =~ 1*y4
#>     q4 =~ 1*x4
#>     p5 =~ 1*y5
#>     q5 =~ 1*x5
#>     y1 ~ 1
#>     x1 ~ 1
#>     y2 ~ 1
#>     x2 ~ 1
#>     y3 ~ 1
#>     x3 ~ 1
#>     y4 ~ 1
#>     x4 ~ 1
#>     y5 ~ 1
#>     x5 ~ 1
#>     p1 ~ 0*1
#>     q1 ~ 0*1
#>     p2 ~ 0*1
#>     q2 ~ 0*1
#>     p3 ~ 0*1
#>     q3 ~ 0*1
#>     p4 ~ 0*1
#>     q4 ~ 0*1
#>     p5 ~ 0*1
#>     q5 ~ 0*1
#>     p2 ~ 0.4*p1 + 0.2*q1
#>     q2 ~ 0.3*q1 + 0.1*p1
#>     p3 ~ 0.4*p2 + 0.2*q2
#>     q3 ~ 0.3*q2 + 0.1*p2
#>     p4 ~ 0.4*p3 + 0.2*q3
#>     q4 ~ 0.3*q3 + 0.1*p3
#>     p5 ~ 0.4*p4 + 0.2*q4
#>     q5 ~ 0.3*q4 + 0.1*p4
#>     p1 ~~ 1*p1
#>     q1 ~~ 1*q1
#>     p1 ~~ 0.3*q1
#>     p2 ~~ 1*p2
#>     q2 ~~ 1*q2
#>     p2 ~~ 0.3*q2
#>     p3 ~~ 1*p3
#>     q3 ~~ 1*q3
#>     p3 ~~ 0.3*q3
#>     p4 ~~ 1*p4
#>     q4 ~~ 1*q4
#>     p4 ~~ 0.3*q4
#>     p5 ~~ 1*p5
#>     q5 ~~ 1*q5
#>     p5 ~~ 0.3*q5
#>     y1 ~~ 0*y1
#>     x1 ~~ 0*x1
#>     y2 ~~ 0*y2
#>     x2 ~~ 0*x2
#>     y3 ~~ 0*y3
#>     x3 ~~ 0*x3
#>     y4 ~~ 0*y4
#>     x4 ~~ 0*x4
#>     y5 ~~ 0*y5
#>     x5 ~~ 0*x5

# View the data
head(sim_result$data)
#>           y1         x1         y2         x2          y3         x3         y4
#> 1  0.1766046  0.8087061  1.7378893 -0.2601594  1.11140003  1.3681487 -0.8035027
#> 2  0.5099755  1.5115138 -0.2219419 -0.4018028 -0.72733658 -1.8541145 -0.6905359
#> 3  0.3793157  1.0991914  0.9293786  3.2643445  0.12129024  2.1901878  0.3523382
#> 4 -2.2606166 -1.4035611 -1.2635553  1.0038483 -1.75338508 -1.0694481 -1.4845073
#> 5 -0.8338965  1.0441023  0.8016214 -0.1334976 -0.03687628 -1.1050282 -2.1377252
#> 6 -0.6032593 -0.9682397 -1.4419061  0.7830082  0.06761841  0.5424185 -1.0280294
#>            x4         y5         x5
#> 1  0.80663950 -1.2011982 -0.6210664
#> 2 -0.98481246 -1.8985272 -0.6588583
#> 3 -0.66466800  0.3430064 -0.5593006
#> 4 -0.93899214 -0.9527972  1.2524949
#> 5 -2.01977664  0.6367791 -1.4151973
#> 6  0.03972786 -0.2319081  0.8289402
```
