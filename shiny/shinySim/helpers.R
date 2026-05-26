# helpers.R — Simulation app helper functions
# Wraps crossLagR sim*/estimate* into a unified interface

# ============================================================================
# DGP Registry: maps DGP name → sim function + default parameters
# ============================================================================

dgp_registry <- list(
  CLPM = list(
    sim_fn = "simCLPM",
    params = list(
      beta_x = 0.3, beta_y = 0.3,
      omega_xy = 0.1, omega_yx = 0.1,
      var_x = 1, var_y = 1, cov_xy = 0.3
    ),
    desc = "Standard CLPM (no trait factors)"
  ),
  `RI-CLPM` = list(
    sim_fn = "simRICLPM",
    params = list(
      beta_x = 0.3, beta_y = 0.3,
      omega_xy = 0.1, omega_yx = 0.1,
      var_p = 1, var_q = 1, cov_pq = 0.1,
      var_BX = 1.5, var_BY = 1.5, cov_BXBY = 0.5
    ),
    desc = "Random Intercept CLPM (trait + state)"
  ),
  `Latent Change` = list(
    sim_fn = "simLChange",
    params = list(
      beta_x = -0.2, beta_y = -0.2,
      omega_x = -0.1, omega_y = -0.1,
      initial_var_x = 1, initial_var_y = 1,
      indicator_variance_x = 0.5, indicator_variance_y = 0.5,
      latent_variance_x = 1, latent_variance_y = 1,
      cov_initial_xy = 0.3
    ),
    desc = "Bivariate Latent Change Score Model"
  )
)

# ============================================================================
# Estimator Registry: maps estimator name → estimate function + extraction
# ============================================================================

estimator_registry <- list(
  CLPM = list(
    est_fn = "estimateCLPM",
    args = list(constrain_beta = FALSE, constrain_omega = FALSE,
                constrain_residual_variances = FALSE,
                constrain_residual_covariances = FALSE),
    extract = "cl_labels"
  ),
  `RI-CLPM` = list(
    est_fn = "estimateRICLPM",
    args = list(constrain_beta = FALSE, constrain_omega = FALSE,
                constrain_residual_variances = FALSE,
                constrain_residual_covariances = FALSE),
    extract = "cl_labels"
  ),
  `Latent Change` = list(
    est_fn = "estimateLChange",
    args = list(variable_type = "bivariate",
                constrain_beta = TRUE, constrain_omega = TRUE,
                constrain_indicator_variances = TRUE,
                estimate_constant_change = TRUE),
    extract = "lchange_coupling"
  ),
  `AC-FI` = list(
    est_fn = "estimateAllisonChamberlainFI",
    args = list(constrain_coefficients = TRUE,
                correlate_alpha_with_exogenous = FALSE),
    extract = "cl_labels"
  ),
  ALT = list(
    est_fn = "estimateALT",
    args = list(constrain_coefficients = TRUE),
    extract = "cl_labels"
  ),
  `Bollen-Brand` = list(
    est_fn = "estimateBollenBrand",
    args = list(constrain_coefficients = TRUE),
    extract = "cl_labels"
  )
)

# ============================================================================
# Pairwise caveat text: DGP × Estimator educational guidance
# ============================================================================

get_caveat_text <- function(dgp, est) {
  if (dgp == est) return(NULL)
  key <- paste0(dgp, " -> ", est)
  caveats <- list(
    # ---- RI-CLPM DGP, wrong estimator ----
    `RI-CLPM -> CLPM` = paste0(
      "You specified an RI-CLPM data-generating process but estimated a standard CLPM. ",
      "The CLPM conflates time-varying (within-person) and time-invariant (trait) effects ",
      "into a single set of coefficients. When stable between-person differences exist, ",
      "the CLPM autoregressive estimates will be biased upward because they absorb trait ",
      "stability, and cross-lagged estimates may be biased in either direction depending ",
      "on the covariance between traits. ",
      "Try this: Increase the ICC (proportion of between-person variance) using the slider ",
      "and watch how the CLPM overestimates the autoregressive parameters."
    ),
    `RI-CLPM -> Latent Change` = paste0(
      "You specified an RI-CLPM DGP but estimated a Latent Change Score model. ",
      "The LCS parameterizes dynamics as proportional effects (levels predicting change) ",
      "and coupling effects, which are conceptually different from CLPM-style cross-lags. ",
      "When the true DGP is an RI-CLPM, the LCS estimates will not map directly onto the ",
      "true parameters. The LCS proportional effect absorbs a mix of autoregressive and ",
      "trait components. Coupling parameters may partially recover cross-lagged effects, ",
      "but the correspondence is inexact."
    ),
    `RI-CLPM -> AC-FI` = paste0(
      "You specified an RI-CLPM DGP and estimated an Allison-Chamberlain Fixed Effects model. ",
      "The AC-FI removes stable between-person confounds (traits) through differencing, ",
      "which should eliminate trait bias. However, AC-FI suffers from Nickell bias: the ",
      "demeaned lag is mechanically correlated with the demeaned error, introducing a downward ",
      "bias on autoregressive parameters that is O(1/T). With short panels (T = 3-5), expect ",
      "the AR estimates to be substantially attenuated. Cross-lagged estimates are also affected. ",
      "Try this: Increase the number of waves and watch the AC-FI estimates converge to truth."
    ),
    `RI-CLPM -> ALT` = paste0(
      "You specified an RI-CLPM DGP but estimated an Autoregressive Latent Trajectory model. ",
      "The ALT models growth trajectories alongside autoregressive dynamics. When the true DGP ",
      "has random intercepts but no growth, the ALT may absorb trait variance into its slope factor, ",
      "leading to biased cross-lagged estimates. The two models partition variance differently."
    ),
    `RI-CLPM -> Bollen-Brand` = paste0(
      "You specified an RI-CLPM DGP but estimated a Bollen-Brand (dynamic panel) model. ",
      "The Bollen-Brand uses predetermined variable constraints for identification, which differs ",
      "from the RI-CLPM's random intercept approach. When traits covary with initial conditions, ",
      "the Bollen-Brand estimates may be inconsistent."
    ),
    # ---- CLPM DGP, wrong estimator ----
    `CLPM -> RI-CLPM` = paste0(
      "You specified a CLPM DGP but estimated an RI-CLPM. Since the true DGP has no ",
      "between-person trait variance, the random intercept variances should estimate near zero. ",
      "The within-person (state) parameters should be close to the truth. This is a good ",
      "sensitivity check: if the RI-CLPM recovers the CLPM parameters and the trait variance ",
      "is near zero, it confirms the CLPM is appropriate. ",
      "Note: With small samples, the RI-CLPM may have convergence issues when fitting ",
      "unnecessary random intercepts."
    ),
    `CLPM -> Latent Change` = paste0(
      "You specified a CLPM DGP but estimated a Latent Change Score model. ",
      "The LCS proportional effects and coupling parameters have a different causal interpretation ",
      "than CLPM cross-lags. The LCS beta (proportional) is not the same as the CLPM beta (AR): ",
      "the LCS models change from the previous level, while the CLPM directly regresses the current ",
      "level on the prior level. The mapping is approximately: LCS beta ≈ CLPM beta - 1."
    ),
    `CLPM -> AC-FI` = paste0(
      "You specified a CLPM DGP (no traits) but estimated an Allison-Chamberlain Fixed Effects model. ",
      "Since there are no between-person confounds to remove, the AC-FI differencing is unnecessary ",
      "and introduces Nickell bias without any offsetting benefit. The AR estimates will be biased ",
      "downward. The CLPM is the correct estimator here. ",
      "This demonstrates why blindly applying fixed effects is not always beneficial."
    ),
    `CLPM -> ALT` = paste0(
      "You specified a CLPM DGP but estimated an ALT model. The ALT includes growth factors that ",
      "don't exist in the DGP. With enough data, the slope variances should estimate near zero. ",
      "However, with smaller samples, the ALT may attribute noise to growth, distorting AR/CL estimates."
    ),
    `CLPM -> Bollen-Brand` = paste0(
      "You specified a CLPM DGP but estimated a Bollen-Brand model. ",
      "The Bollen-Brand uses predetermined variable constraints. When the DGP is a simple CLPM, ",
      "the Bollen-Brand should approximately recover the true parameters, but may differ due to ",
      "its different identification strategy."
    ),
    # ---- Latent Change DGP, wrong estimator ----
    `Latent Change -> CLPM` = paste0(
      "You specified a Latent Change DGP but estimated a CLPM. The CLPM AR coefficient ",
      "conflates the LCS proportional effect with the unit autoregression: CLPM beta ≈ 1 + LCS beta. ",
      "Cross-lagged parameters correspond roughly to the LCS coupling effects, but the mapping is ",
      "inexact because the CLPM does not model change scores explicitly."
    ),
    `Latent Change -> RI-CLPM` = paste0(
      "You specified a Latent Change DGP but estimated an RI-CLPM. The RI-CLPM will try to ",
      "partition variance into random intercepts and within-person dynamics. The LCS constant change ",
      "factor may be partially absorbed into the random intercept. Within-person AR and CL estimates ",
      "will not directly correspond to LCS proportional and coupling effects."
    ),
    `Latent Change -> AC-FI` = paste0(
      "You specified a Latent Change DGP but estimated an AC-FI model. The fixed effects ",
      "differencing removes person-mean levels, which partially corresponds to the LCS framework. ",
      "However, the AC-FI AR/CL interpretation differs from LCS proportional/coupling effects, and ",
      "Nickell bias will attenuate the autoregressive estimates with short panels."
    ),
    `Latent Change -> ALT` = paste0(
      "You specified a Latent Change DGP but estimated an ALT model. Both models include growth ",
      "components, but they parameterize dynamics differently. The ALT growth slope may absorb the ",
      "LCS constant change factor, while the AR/CL components capture residual dynamics."
    ),
    `Latent Change -> Bollen-Brand` = paste0(
      "You specified a Latent Change DGP but estimated a Bollen-Brand model. ",
      "The Bollen-Brand will estimate autoregressive and cross-lagged effects that are not directly ",
      "comparable to LCS proportional and coupling effects."
    )
  )
  caveats[[key]]
}

# ============================================================================
# Simulate one dataset
# ============================================================================

simulate_data <- function(dgp_name, waves, sample_size, params) {
  reg <- dgp_registry[[dgp_name]]
  fn <- match.fun(reg$sim_fn)

  # Build args
  call_args <- c(list(waves = waves), params)
  if (dgp_name == "CLPM") {
    call_args$sample_size <- sample_size
  } else if (dgp_name == "Latent Change") {
    call_args$sample.nobs <- sample_size
    call_args$variable_type <- "bivariate"
  } else {
    call_args$sample.nobs <- sample_size
  }

  result <- do.call(fn, call_args)
  result$data
}

# ============================================================================
# Fit one estimator
# ============================================================================

fit_estimator <- function(est_name, data, waves) {
  reg <- estimator_registry[[est_name]]
  fn <- match.fun(reg$est_fn)

  call_args <- c(list(waves = waves), reg$args)
  syntax <- do.call(fn, call_args)

  fit <- tryCatch(
    suppressWarnings(lavaan::lavaan(
      syntax, data = data, missing = "fiml",
      meanstructure = TRUE, int.ov.free = TRUE
    )),
    error = function(e) NULL
  )

  if (is.null(fit) || !lavaan::lavInspect(fit, "converged")) {
    return(list(converged = FALSE, params = data.frame()))
  }

  # Extract parameters based on extraction method
  ss <- lavaan::standardizedSolution(fit)
  pe <- lavaan::parameterEstimates(fit)

  if (reg$extract == "cl_labels") {
    params_df <- ss |>
      dplyr::filter(op == "~", grepl("^(cl_|ar_|beta_|omega_)", label)) |>
      dplyr::select(label, est.std, se, pvalue, ci.lower, ci.upper)
  } else if (reg$extract == "lchange_coupling") {
    params_df <- ss |>
      dplyr::filter(op == "~") |>
      dplyr::filter(
        grepl("^ld_", lhs) & grepl("^cf_", rhs)
      ) |>
      dplyr::mutate(
        label = dplyr::case_when(
          grepl("^ld_x", lhs) & grepl("^cf_x", rhs) ~ "beta_x",
          grepl("^ld_y", lhs) & grepl("^cf_y", rhs) ~ "beta_y",
          grepl("^ld_x", lhs) & grepl("^cf_y", rhs) ~ "omega_x",
          grepl("^ld_y", lhs) & grepl("^cf_x", rhs) ~ "omega_y",
          TRUE ~ paste0(lhs, "~", rhs)
        )
      ) |>
      dplyr::group_by(label) |>
      dplyr::summarise(
        est.std = mean(est.std), se = mean(se),
        pvalue = mean(pvalue),
        ci.lower = mean(ci.lower), ci.upper = mean(ci.upper),
        .groups = "drop"
      )
  }

  list(
    converged = TRUE,
    params = params_df,
    fit_measures = tryCatch(
      lavaan::fitMeasures(fit, c("cfi", "rmsea", "srmr")),
      error = function(e) c(cfi = NA, rmsea = NA, srmr = NA)
    )
  )
}

# ============================================================================
# Get true parameter values from DGP for comparison
# ============================================================================

get_true_values <- function(dgp_name, params) {
  if (dgp_name %in% c("CLPM", "RI-CLPM")) {
    data.frame(
      base_label = c("ar_x", "ar_y", "cl_xy", "cl_yx",
                      "beta_x", "beta_y", "omega_xy", "omega_yx"),
      true = c(params$beta_x, params$beta_y,
               params$omega_xy, params$omega_yx,
               params$beta_x, params$beta_y,
               params$omega_xy, params$omega_yx)
    )
  } else if (dgp_name == "Latent Change") {
    data.frame(
      base_label = c("beta_x", "beta_y", "omega_x", "omega_y"),
      true = c(params$beta_x, params$beta_y,
               params$omega_x, params$omega_y)
    )
  }
}

# Strip trailing wave number from estimated label to get the base parameter name
# e.g., "ar_y2" -> "ar_y", "cl_xy3" -> "cl_xy", "beta_x" -> "beta_x"
strip_wave <- function(label) {
  gsub("[0-9]+$", "", label)
}

# Join estimated params with true values by base label
match_true_values <- function(est_params, true_vals) {
  est_params |>
    dplyr::mutate(base_label = strip_wave(label)) |>
    dplyr::left_join(true_vals, by = "base_label") |>
    dplyr::mutate(
      bias = round(est.std - true, 4),
      dplyr::across(c(est.std, se, true), ~ round(.x, 4))
    )
}


# ============================================================================
# MC summary statistics
# ============================================================================

mc_summary <- function(mc_results, true_vals) {
  mc_results |>
    dplyr::mutate(base_label = strip_wave(label)) |>
    dplyr::group_by(estimator, base_label) |>
    dplyr::summarise(
      mean_est = mean(est.std, na.rm = TRUE),
      sd_est = sd(est.std, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::left_join(true_vals, by = "base_label") |>
    dplyr::mutate(
      bias = mean_est - true,
      rmse = sqrt(bias^2 + sd_est^2),
      rel_bias = dplyr::if_else(abs(true) > 0.001, bias / true, NA_real_)
    )
}
