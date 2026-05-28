# ISCAP panel: party identification ↔ ideology

This vignette walks through a minimal cross-lagged analysis of the ISCAP
panel \[@upenn_iscap_2022\]: do shifts in party identification predict
shifts in ideological self-placement, or the reverse? See the [user
guide](https://crweber9874.github.io/crossLagR/book/06_iscap.html) for
the full treatment.

## Setup

``` r

library(crossLagR)
library(lavaan)
library(dplyr)
library(tidyr)

data(upenn)
```

The ISCAP panel includes three waves we’ll use here (2012, 2016, 2018).
The relevant columns are `party_identification_*` (7-point, 7 = Strong
Republican) and `ideology_*` (7-point, 7 = Extremely Conservative).

## Wide data

``` r

wide <- upenn |>
  transmute(
    x1 = party_identification_2012,
    x2 = party_identification_2016,
    x3 = party_identification_2018,
    y1 = ideology_2012,
    y2 = ideology_2016,
    y3 = ideology_2018
  ) |>
  drop_na()

n_obs  <- nrow(wide)
n_miss <- sum(is.na(upenn[, c("party_identification_2012",
                              "party_identification_2016",
                              "party_identification_2018",
                              "ideology_2012",
                              "ideology_2016",
                              "ideology_2018")]))

cat("Complete cases used:", n_obs, " (", n_miss, "missing cells in source)\n")
#> Complete cases used: 1082  ( 4671 missing cells in source)
```

## CLPM

[`estimateCLPM()`](https://crweber9874.github.io/crossLagR/reference/estimateCLPM.md)
returns lavaan model syntax with the conventional `beta_*`
(autoregressive) and `omega_*` (cross-lag) labels.

``` r

clpm_fit <- lavaan(estimateCLPM(waves = 3), data = wide)
summary(clpm_fit, standardized = TRUE)
#> lavaan 0.6-21 ended normally after 44 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                        19
#>   Number of equality constraints                     7
#> 
#>   Number of observations                          1082
#> 
#> Model Test User Model:
#>                                                       
#>   Test statistic                               698.327
#>   Degrees of freedom                                15
#>   P-value (Chi-square)                           0.000
#> 
#> Parameter Estimates:
#> 
#>   Standard errors                             Standard
#>   Information                                 Expected
#>   Information saturated (h1) model          Structured
#> 
#> Latent Variables:
#>                    Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#>   p1 =~                                                                 
#>     y1                1.000                               1.531    1.000
#>   q1 =~                                                                 
#>     x1                1.000                               2.154    1.000
#>   p2 =~                                                                 
#>     y2                1.000                               1.704    1.000
#>   q2 =~                                                                 
#>     x2                1.000                               2.176    1.000
#>   p3 =~                                                                 
#>     y3                1.000                               1.853    1.000
#>   q3 =~                                                                 
#>     x3                1.000                               2.224    1.000
#> 
#> Regressions:
#>                    Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#>   p2 ~                                                                  
#>     p1      (ar_y)    0.916    0.010   94.344    0.000    0.823    0.823
#>     q1      (cl_x)    0.090    0.010    9.134    0.000    0.114    0.114
#>   q2 ~                                                                  
#>     q1      (ar_x)    0.736    0.015   48.361    0.000    0.728    0.728
#>     p1      (cl_y)    0.250    0.015   16.761    0.000    0.176    0.176
#>   p3 ~                                                                  
#>     p2      (ar_y)    0.916    0.010   94.344    0.000    0.843    0.843
#>     q2      (cl_x)    0.090    0.010    9.134    0.000    0.106    0.106
#>   q3 ~                                                                  
#>     q2      (ar_x)    0.736    0.015   48.361    0.000    0.720    0.720
#>     p2      (cl_y)    0.250    0.015   16.761    0.000    0.192    0.192
#> 
#> Covariances:
#>                    Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#>   p1 ~~                                                                 
#>     q1      (d__1)    2.223    0.121   18.390    0.000    0.674    0.674
#>  .p2 ~~                                                                 
#>    .q2      (d_c_)    0.158    0.018    8.793    0.000    0.193    0.193
#>  .p3 ~~                                                                 
#>    .q3      (d_c_)    0.158    0.018    8.793    0.000    0.193    0.193
#> 
#> Intercepts:
#>                    Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#>     p1                4.201    0.047   90.262    0.000    2.744    2.744
#>     q1                3.834    0.065   58.556    0.000    1.780    1.780
#>    .p2                0.000                               0.000    0.000
#>    .q2                0.000                               0.000    0.000
#>    .p3                0.000                               0.000    0.000
#>    .q3                0.000                               0.000    0.000
#>    .y1                0.000                               0.000    0.000
#>    .x1                0.000                               0.000    0.000
#>    .y2                0.000                               0.000    0.000
#>    .x2                0.000                               0.000    0.000
#>    .y3                0.000                               0.000    0.000
#>    .x3                0.000                               0.000    0.000
#> 
#> Variances:
#>                    Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#>     p1   (d_vr_y1)    2.343    0.101   23.259    0.000    1.000    1.000
#>     q1   (d_vr_x1)    4.638    0.199   23.259    0.000    1.000    1.000
#>    .p2    (d_vr_y)    0.532    0.016   32.894    0.000    0.183    0.183
#>    .q2    (d_vr_x)    1.259    0.038   32.894    0.000    0.266    0.266
#>    .p3    (d_vr_y)    0.532    0.016   32.894    0.000    0.155    0.155
#>    .q3    (d_vr_x)    1.259    0.038   32.894    0.000    0.255    0.255
#>    .y1                0.000                               0.000    0.000
#>    .x1                0.000                               0.000    0.000
#>    .y2                0.000                               0.000    0.000
#>    .x2                0.000                               0.000    0.000
#>    .y3                0.000                               0.000    0.000
#>    .x3                0.000                               0.000    0.000
```

## RI-CLPM

[`estimateRICLPM()`](https://crweber9874.github.io/crossLagR/reference/estimateRICLPM.md)
adds random intercepts $`I_x, I_y`$ that absorb stable between-person
differences. The AR and CL coefficients are now estimated on the
within-person deviations.

``` r

riclpm_fit <- lavaan(estimateRICLPM(waves = 3), data = wide)
summary(riclpm_fit, standardized = TRUE)
#> lavaan 0.6-21 ended normally after 44 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                        22
#>   Number of equality constraints                    10
#> 
#>   Number of observations                          1082
#> 
#> Model Test User Model:
#>                                                       
#>   Test statistic                                59.581
#>   Degrees of freedom                                15
#>   P-value (Chi-square)                           0.000
#> 
#> Parameter Estimates:
#> 
#>   Standard errors                             Standard
#>   Information                                 Expected
#>   Information saturated (h1) model          Structured
#> 
#> Latent Variables:
#>                    Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#>   I_x =~                                                                
#>     x1                1.000                               2.046    0.939
#>     x2                1.000                               2.046    0.936
#>     x3                1.000                               2.046    0.936
#>   I_y =~                                                                
#>     y1                1.000                               1.412    0.922
#>     y2                1.000                               1.412    0.917
#>     y3                1.000                               1.412    0.917
#>   p1 =~                                                                 
#>     x1                1.000                               0.747    0.343
#>   q1 =~                                                                 
#>     y1                1.000                               0.595    0.388
#>   p2 =~                                                                 
#>     x2                1.000                               0.768    0.352
#>   q2 =~                                                                 
#>     y2                1.000                               0.614    0.399
#>   p3 =~                                                                 
#>     x3                1.000                               0.770    0.352
#>   q3 =~                                                                 
#>     y3                1.000                               0.615    0.399
#> 
#> Regressions:
#>                    Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#>   p2 ~                                                                  
#>     p1      (ar_x)   -0.245    0.035   -6.964    0.000   -0.238   -0.238
#>     q1      (cl_y)    0.050    0.054    0.926    0.355    0.039    0.039
#>   q2 ~                                                                  
#>     q1      (ar_y)    0.247    0.057    4.320    0.000    0.239    0.239
#>     p1      (cl_x)    0.022    0.025    0.876    0.381    0.027    0.027
#>   p3 ~                                                                  
#>     p2      (ar_x)   -0.245    0.035   -6.964    0.000   -0.245   -0.245
#>     q2      (cl_y)    0.050    0.054    0.926    0.355    0.040    0.040
#>   q3 ~                                                                  
#>     q2      (ar_y)    0.247    0.057    4.320    0.000    0.246    0.246
#>     p2      (cl_x)    0.022    0.025    0.876    0.381    0.028    0.028
#> 
#> Covariances:
#>                    Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#>   p1 ~~                                                                 
#>     q1      (d_c_)    0.066    0.014    4.742    0.000    0.149    0.149
#>  .p2 ~~                                                                 
#>    .q2      (d_c_)    0.066    0.014    4.742    0.000    0.149    0.149
#>  .p3 ~~                                                                 
#>    .q3      (d_c_)    0.066    0.014    4.742    0.000    0.149    0.149
#>   I_x ~~                                                                
#>     I_y               2.237    0.116   19.260    0.000    0.775    0.775
#> 
#> Intercepts:
#>                    Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#>     p1                0.000                               0.000    0.000
#>     q1                0.000                               0.000    0.000
#>    .p2                0.000                               0.000    0.000
#>    .q2                0.000                               0.000    0.000
#>    .p3                0.000                               0.000    0.000
#>    .q3                0.000                               0.000    0.000
#>     I_x               3.876    0.063   61.332    0.000    1.895    1.895
#>     I_y               4.250    0.045   95.144    0.000    3.011    3.011
#>    .x1                0.000                               0.000    0.000
#>    .y1                0.000                               0.000    0.000
#>    .x2                0.000                               0.000    0.000
#>    .y2                0.000                               0.000    0.000
#>    .x3                0.000                               0.000    0.000
#>    .y3                0.000                               0.000    0.000
#> 
#> Variances:
#>                    Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#>     p1    (d_vr_x)    0.557    0.021   26.715    0.000    1.000    1.000
#>     q1    (d_vr_y)    0.354    0.018   19.809    0.000    1.000    1.000
#>    .p2    (d_vr_x)    0.557    0.021   26.715    0.000    0.944    0.944
#>    .q2    (d_vr_y)    0.354    0.018   19.809    0.000    0.940    0.940
#>    .p3    (d_vr_x)    0.557    0.021   26.715    0.000    0.941    0.941
#>    .q3    (d_vr_y)    0.354    0.018   19.809    0.000    0.937    0.937
#>     I_x               4.185    0.186   22.489    0.000    1.000    1.000
#>     I_y               1.993    0.095   20.946    0.000    1.000    1.000
#>    .x1                0.000                               0.000    0.000
#>    .y1                0.000                               0.000    0.000
#>    .x2                0.000                               0.000    0.000
#>    .y2                0.000                               0.000    0.000
#>    .x3                0.000                               0.000    0.000
#>    .y3                0.000                               0.000    0.000
```

## Reading the contrast

A pattern that recurs in the book:

- The CLPM AR coefficients are inflated by stable between-person trait
  variance and the CL coefficients tend to be conservatively biased.
- The RI-CLPM AR coefficients are smaller (within-person carry-over
  only). CL coefficients show whether **change** in one variable
  forecasts **change** in the other.

For the ISCAP data the cross-lags are modest in either specification,
with party identification appearing to lead ideology slightly more than
the reverse. See [the book
chapter](https://crweber9874.github.io/crossLagR/book/06_iscap.html) for
the full discussion including a system-legitimacy block, latent-change
specifications, and Bayesian fits via `brms`.

## References
