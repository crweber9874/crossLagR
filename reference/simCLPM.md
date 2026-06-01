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
#>           y1         x1          y2          x2         y3         x3
#> 1  0.4826720  0.4176978  1.80522025 -0.13900572  0.8236712  1.7126768
#> 2  1.5928859  0.1280702  0.01628455  0.02685647 -1.7453624 -0.6351248
#> 3 -0.6404779  2.4020017  0.70503707  2.86066936  1.0799810  1.0422462
#> 4 -2.0884470 -1.6235117 -1.22568018  1.07199986 -1.9152387 -0.8756436
#> 5 -0.9950844  1.2500236  0.76616209 -0.19730226  0.1146537 -1.2864711
#> 6 -0.8305783 -0.6778342 -1.49191334  0.69302617  0.2813173  0.2865343
#>           y4         x4         y5         x5
#> 1 -0.9582472  0.6455208 -0.8326336 -0.8243582
#> 2 -1.2380441 -1.5548735 -0.5944928 -1.3781340
#> 3  0.8679352 -0.1278327 -0.8850231  0.1180525
#> 4 -1.5715544 -1.0296249 -0.7454716  1.1381389
#> 5 -2.0562303 -1.9349248  0.4426776 -1.3081353
#> 6 -0.9130992  0.1593922 -0.5056444  0.9799269
```
