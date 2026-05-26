# crossLagR Variable and Parameter Conventions

## Latent Variable Naming

The package uses short letters for latent within-person factors across all CLPM and RICLPM functions. These appear in lavaan model syntax and in simulated data output.

| Symbol | Represents | Notes |
|--------|-----------|-------|
| `p`    | Within-person latent factor for X (first variable) | `p1`, `p2`, ..., `p[waves]` |
| `q`    | Within-person latent factor for Y (second variable) | `q1`, `q2`, ..., `q[waves]` |
| `r`    | Within-person latent factor for Z (third variable, optional) | `r1`, `r2`, ..., `r[waves]` |

These latent factors are defined via measurement equations like `p1 =~ 1*x1`, meaning each observed variable is a perfect indicator of its latent factor.

## Observed Variable Naming

| Symbol | Represents | Used in |
|--------|-----------|---------|
| `x1, x2, ..., x[waves]` | Observed X variable at each wave | All estimation and simulation functions |
| `y1, y2, ..., y[waves]` | Observed Y variable at each wave | All estimation and simulation functions |
| `z1, z2, ..., z[waves]` | Observed Z variable at each wave (3-variable models) | `estimateRICLPM(include_z=TRUE)`, `estimateRICLPM3()` |

**Your data must use these exact column names** (or be renamed to match) before passing to `lavaan()`.

## Random Intercepts (RICLPM only)

| Symbol | Represents |
|--------|-----------|
| `BX`   | Between-person random intercept for X |
| `BY`   | Between-person random intercept for Y |
| `BZ`   | Between-person random intercept for Z |

These are latent variables with unit loadings on all observed waves: `BX =~ 1*x1 + 1*x2 + ... + 1*x[waves]`. They capture stable individual differences that are constant across time.

## Latent Fixed Effects (Bollen & Brand only)

| Symbol | Represents |
|--------|-----------|
| `eta_y` | Unobserved heterogeneity for Y |
| `eta_x` | Unobserved heterogeneity for X |

These are analogous to unit fixed effects in econometrics. They have unit loadings on endogenous outcomes (waves 2+) and are allowed to correlate freely with exogenous variables.

## Parameter Labels

### Autoregressive (Stability) Parameters

| Label | Meaning | Direction |
|-------|---------|-----------|
| `beta_x` or `rho_x` | Autoregressive effect for X | X(t-1) -> X(t) |
| `beta_y` or `rho_y` | Autoregressive effect for Y | Y(t-1) -> Y(t) |
| `beta_z` or `rho_z` | Autoregressive effect for Z | Z(t-1) -> Z(t) |

- `beta_*` is used in CLPM and RICLPM functions
- `rho_*` is used in Bollen & Brand functions
- When constrained, these are equal across all waves
- When unconstrained, wave-specific labels are used (e.g., `beta_x2`, `beta_x3`)

### Cross-Lagged Parameters

| Label | Meaning | Direction |
|-------|---------|-----------|
| `omega_xy` | Cross-lagged effect of X on Y | X(t-1) -> Y(t) |
| `omega_yx` | Cross-lagged effect of Y on X | Y(t-1) -> X(t) |
| `omega_xz` | Cross-lagged effect of X on Z | X(t-1) -> Z(t) |
| `omega_yz` | Cross-lagged effect of Y on Z | Y(t-1) -> Z(t) |
| `omega_zx` | Cross-lagged effect of Z on X | Z(t-1) -> X(t) |
| `omega_zy` | Cross-lagged effect of Z on Y | Z(t-1) -> Y(t) |
| `b1`       | Cross-lagged X -> Y (Bollen & Brand) | X(t-1) -> Y(t) |
| `b2`       | Cross-lagged Y -> X (Bollen & Brand) | Y(t-1) -> X(t) |

**Reading convention**: `omega_xy` means "the effect of X on Y" (X predicts Y). The subscript order is `[from][to]`.

### Variance and Covariance Parameters

| Label | Meaning |
|-------|---------|
| `var_p`, `var_q`, `var_r` | Residual variances for within-person latent factors (RICLPM) |
| `var_x`, `var_y` | Residual variances for observed variables (CLPM, Bollen & Brand) |
| `var_x1`, `var_y1` | First-wave variances (freely estimated, not constrained with later waves) |
| `cov_pq` | Residual covariance between p and q within the same wave |
| `cov_xy`, `cov_xy1` | Covariance between x and y (general / first wave) |
| `var_BX`, `var_BY` | Random intercept variances (simulation functions) |
| `cov_BXBY` | Random intercept covariance (simulation functions) |

## Common Function Parameters

### Shared Across All Estimation Functions

| Parameter | Type | Default | Meaning |
|-----------|------|---------|---------|
| `waves` | integer | varies | Number of time points. Must be >= 2 (CLPM) or >= 3 (Bollen & Brand) |

### Constraint Parameters (Estimation Functions)

| Parameter | Type | Default | Effect when TRUE |
|-----------|------|---------|-----------------|
| `constrain_beta` | logical | TRUE | AR effects equal across waves |
| `constrain_omega` | logical | TRUE | Cross-lagged effects equal across waves |
| `constrain_coefficients` | logical | TRUE | Both AR and CL constrained (Bollen & Brand) |
| `constrain_residual_variances` | logical | TRUE | Residual variances equal across waves (wave 2+) |
| `constrain_residual_covariances` | logical | TRUE | Residual covariances equal across waves |
| `estimate_means` | logical | FALSE | Include mean structure in the model |
| `start_values` | logical | FALSE | Provide lavaan starting values |

### Simulation Function Parameters

| Parameter | Type | Default | Meaning |
|-----------|------|---------|---------|
| `sample_size` or `sample.nobs` | integer | 1000 | Number of simulated participants |
| `seed` | integer | NULL | Random seed for reproducibility |
| `beta_x`, `beta_y` | numeric | 0.3 | True autoregressive effects |
| `omega_xy`, `omega_yx` | numeric | 0.1 | True cross-lagged effects |

### Monte Carlo Function Parameters

| Parameter | Type | Default | Meaning |
|-----------|------|---------|---------|
| `trials` | integer | varies | Number of simulation replications |
| `estimator` | character | varies | Which model to fit to the data |
| `dgp` or `data_generation` | character | varies | Data-generating process |
| `verbose` | logical | TRUE | Print progress messages |

## Function Families

### 1. Estimation: `estimate*()`
Generate lavaan model syntax strings. Pass the returned string to `lavaan::lavaan()` or `lavaan::sem()`.

| Function | Model | Key Feature |
|----------|-------|-------------|
| `estimateCLPM()` | Cross-Lagged Panel Model | No random intercepts; all effects are within-person |
| `estimateRICLPM()` | Random Intercept CLPM | Separates between-person (BX, BY) from within-person (p, q) |
| `estimateRICLPM3()` | 3-Variable RICLPM | Adds Z variable with custom variable names |
| `estimateRICLPMl()` | RICLPM with multiple indicators | Multiple observed indicators per latent factor |
| `estimateRICLPM_nolag()` | RICLPM without lag structure | Contemporaneous associations only |
| `estimateBollen_and_Brand()` | Dynamic Panel Model | Latent fixed effects (eta), conditions on first obs |
| `estimateAllisonChamberlainFI()` | Allison-Chamberlain FE | Chamberlain predetermined-variable correction |
| `estimateLChange()` | Latent Change Score | Models wave-to-wave change as latent variable |
| `estimateLGM()` | Latent Growth Model | Intercept and slope factors |
| `estimateSTARTS()` | STARTS Model | Stable Trait, AR Trait, State decomposition |
| `estimateHLM()` | Hierarchical Linear Model | brms-based random effects model |
| `estimateCTSEM()` | Continuous-Time SEM | Handles unequal time spacing |
| `estimateFI()` | Fixed Individual Effects | OLS-based within estimator |
| `estimateRI()` | Random Intercept (brms) | brms-based random intercept model |

### 2. Simulation: `sim*()`
Generate synthetic data under known parameters. Returns a list with `$model` (syntax) and `$data` (data frame).

| Function | Generates Data From |
|----------|-------------------|
| `simCLPM()` | CLPM |
| `simRICLPM()` | RICLPM (2 or 3 variables) |
| `simRICLPM3()` | 3-variable RICLPM |
| `simRICLPMl()` | RICLPM with multiple indicators |
| `simLChange()` | Latent Change Score |
| `simLGM()` | Latent Growth Model |
| `simCLPMu()` | CLPM with unmeasured confounder |
| `simCLPM_timeInvariant()` | CLPM with time-invariant predictors |
| `sim_trait_state()` | Trait-state model |

### 3. Monte Carlo: `monteCarlo*()`
Run repeated simulation-and-estimation to evaluate model performance.

| Function | Tests |
|----------|-------|
| `monteCarloCLPM()` | CLPM estimation under various DGPs |
| `monteCarloRICLPM()` | RICLPM estimation |
| `monteCarloCTSEM()` | Continuous-time SEM |
| `monteCarloLChange()` | Latent change score |
| `monteCarloOLS()` | OLS baseline |
| `monteCarloFixed()` | Fixed effects |
| `monteCarloRI()` | Random intercepts |
| `monteCarloConfounder()` | Confounder scenarios |
| `monteCarloAllisonChamberlainFI()` | Allison-Chamberlain FE |
| `run_mc_sims()` | Unified MC orchestrator |

### 4. Utilities

| Function | Purpose |
|----------|---------|
| `zero.one()` | Min-max normalize a vector to [0, 1] |
| `reshape_long_sim_cr()` | Reshape wide simulated data to long format with lagged variables |
| `generate_riclpm_indicator_lists_for_estimation()` | Generate indicator name lists for multi-indicator RICLPM |
| `riclpm_graph()` | Create tidySEM path diagram of RICLPM |
| `compCTSEM()` | Compile CTSEM model |
| `build_book()` | Build bookdown documentation |
| `preview_book()` | Preview documentation in browser |

## Typical Workflow

```r
# 1. Generate model syntax
syntax <- estimateRICLPM(waves = 5)

# 2. Fit the model with lavaan
fit <- lavaan::lavaan(syntax, data = my_data, missing = "ML")

# 3. Examine results
summary(fit, fit.measures = TRUE, standardized = TRUE)

# Or for Bollen & Brand:
bb_syntax <- estimateBollen_and_Brand(waves = 5, x_effect = "lagged",
                                       x_autoregression = TRUE, y_effect_on_x = TRUE)
bb_fit <- lavaan::sem(bb_syntax, data = my_data, missing = "ML")
```

## Data Requirements

All lavaan-based estimation functions expect **wide-format data** with columns named `x1, x2, ..., x[waves]` and `y1, y2, ..., y[waves]`. If your data uses different names (e.g., `party_id_2012`, `party_id_2016`), rename them before fitting:

```r
my_data <- my_data |> rename(x1 = party_id_2012, x2 = party_id_2016, ...)
```
