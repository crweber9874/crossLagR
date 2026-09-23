test_that("simCLPMu honours sample.nobs for both confounder types", {
  expect_equal(nrow(simCLPMu(waves = 3, sample.nobs = 137)$data), 137)
  expect_equal(
    nrow(simCLPMu(waves = 3, confounder_type = "time_invariant",
                  sample.nobs = 137)$data),
    137
  )
})

test_that("confounder_type switches the generated lavaan syntax", {
  tv <- simCLPMu(waves = 3, sample.nobs = 50)$model
  ti <- simCLPMu(waves = 3, confounder_type = "time_invariant",
                 sample.nobs = 50)$model

  # Time-variant: one confounder per wave, with an AR(1) path between them.
  expect_match(tv, "u2 ~ 0\\.4 \\* u1")
  expect_false(grepl("U =~", tv))

  # Time-invariant: a single latent U loading on every observed variable.
  expect_match(ti, "U =~")
  expect_false(grepl("u2 ~ ", ti))
})

test_that("simCLPMu rejects an unknown confounder_type", {
  expect_error(simCLPMu(waves = 3, confounder_type = "nonsense"))
})

test_that("simCLPM_timeInvariantU still works but is deprecated", {
  expect_warning(out <- simCLPM_timeInvariantU(waves = 3, sample.nobs = 50),
                 class = "deprecatedWarning")
  expect_equal(nrow(out$data), 50)
  expect_match(out$model, "U =~")
})

test_that("run_mc_sims passes confounder_type through to the clpmu DGP", {
  res <- suppressWarnings(suppressMessages(run_mc_sims(
    estimator       = "CLPM",
    data_generation = "clpmu",
    param_grid      = data.frame(confounder_type = "time_invariant"),
    trials = 2, waves = 3, sample_size = 200, verbose = FALSE
  )))
  expect_equal(nrow(res), 2)
  expect_true(all(res$n_obs == 200))
  expect_false(any(res$error_occurred))
})

test_that("riclpm_type = 'riclpm_nolag' reaches estimateRICLPM_nolag", {
  res <- suppressWarnings(suppressMessages(run_mc_sims(
    estimator   = "RICLPM",
    riclpm_type = "riclpm_nolag",
    param_grid  = data.frame(stability_p = 0.2),
    trials = 2, waves = 3, sample_size = 200, verbose = FALSE
  )))
  expect_true(all(res$estimator == "RICLPM_NOLAG"))
  expect_true(all(res$riclpm_type == "riclpm_nolag"))
})

test_that("run_mc_sims no longer advertises the removed FI estimator", {
  expect_error(run_mc_sims(estimator = "FI"), "estimator must be one of")
})
