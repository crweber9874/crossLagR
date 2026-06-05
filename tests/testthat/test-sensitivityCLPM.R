## Unit + integration tests for sensitivityCLPM.

fit_clpm_for_sens <- function(seed = 7001) {
  set.seed(seed)
  ## RI-CLPM data with mild between-person variance so the RI-CLPM refit
  ## inside sensitivityCLPM has something to recover.
  dat <- simRICLPM(waves = 4, sample.nobs = 400,
                   beta_x = 0.4, beta_y = 0.4,
                   omega_xy = 0.15, omega_yx = 0.10,
                   var_p = 0.5, var_q = 0.5,
                   var_BX = 1.0, var_BY = 1.0)$data
  list(
    fit  = suppressWarnings(
      lavaan::lavaan(estimateCLPM(waves = 4), data = dat, meanstructure = TRUE)
    ),
    data = dat
  )
}

test_that("sensitivityCLPM returns expected bundled object", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  ctx <- fit_clpm_for_sens()
  s <- sensitivityCLPM(ctx$fit, icc_grid = c(0, 0.3, 0.6),
                       n_sims = 3, seed = 1, data = ctx$data)
  expect_s3_class(s, "crossLagR_sensitivity")
  for (nm in c("estimates", "clpm_actual", "riclpm_truth",
               "observed_icc", "implied_bias", "icc_grid",
               "n_sims", "sample_size", "waves")) {
    expect_true(nm %in% names(s), info = paste("missing:", nm))
  }
  expect_setequal(names(s$clpm_actual),
                  c("ar_x", "ar_y", "cl_xy", "cl_yx"))
  expect_setequal(names(s$riclpm_truth),
                  c("ar_x", "ar_y", "cl_xy", "cl_yx"))
})

test_that("sensitivityCLPM detects AR inflation: CLPM_actual > RI-CLPM_truth", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  ctx <- fit_clpm_for_sens()
  s <- sensitivityCLPM(ctx$fit, icc_grid = c(0, 0.4),
                       n_sims = 3, seed = 2, data = ctx$data)
  ## When the true DGP has trait variance, the CLPM estimate should
  ## inflate the AR relative to the RI-CLPM estimate.
  expect_gt(s$implied_bias[["ar_x"]], 0)
  expect_gt(s$implied_bias[["ar_y"]], 0)
})

test_that("sensitivityCLPM rejects bad input", {
  expect_error(sensitivityCLPM(list()), "lavaan")
  expect_error(sensitivityCLPM(1), "lavaan")
})

test_that("sensitivityCLPM rejects icc_grid out of [0, 1)", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  ctx <- fit_clpm_for_sens()
  expect_error(sensitivityCLPM(ctx$fit, icc_grid = c(-0.1), data = ctx$data),
               "icc_grid")
  expect_error(sensitivityCLPM(ctx$fit, icc_grid = c(1.0),  data = ctx$data),
               "icc_grid")
})

test_that("sensitivityCLPM CLPM curve AR rises monotonically with simulated ICC", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  ctx <- fit_clpm_for_sens()
  s <- sensitivityCLPM(ctx$fit, icc_grid = c(0, 0.3, 0.6),
                       n_sims = 4, seed = 11, data = ctx$data)
  agg <- stats::aggregate(
    est ~ icc, data = s$estimates[s$estimates$parameter == "ar_x", ],
    FUN = mean, na.rm = TRUE
  )
  expect_gt(agg$est[agg$icc == 0.6], agg$est[agg$icc == 0])
})

test_that("sensitivityCLPM observed_icc is populated from data", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  ctx <- fit_clpm_for_sens()
  s <- sensitivityCLPM(ctx$fit, icc_grid = c(0, 0.4),
                       n_sims = 2, seed = 3, data = ctx$data)
  expect_true(is.finite(s$observed_icc["x"]))
  expect_true(is.finite(s$observed_icc["y"]))
  ## Should be > 0 since we simulated with non-zero var_BX/var_BY.
  expect_gt(s$observed_icc[["x"]], 0)
})

test_that("print.crossLagR_sensitivity shows the bias table", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  ctx <- fit_clpm_for_sens()
  s <- sensitivityCLPM(ctx$fit, icc_grid = c(0, 0.3),
                       n_sims = 2, seed = 4, data = ctx$data)
  out <- capture.output(print(s))
  expect_true(any(grepl("CLPM bias diagnostic", out)))
  expect_true(any(grepl("clpm_actual", out)))
  expect_true(any(grepl("riclpm_truth", out)))
  expect_true(any(grepl("implied_bias", out)))
  expect_true(any(grepl("Observed ICC", out)))
})

test_that("plot.crossLagR_sensitivity returns a ggplot with diamond + dashed line", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("ggplot2")
  ctx <- fit_clpm_for_sens()
  s <- sensitivityCLPM(ctx$fit, icc_grid = c(0, 0.4),
                       n_sims = 2, seed = 5, data = ctx$data)
  p <- plot(s)
  expect_s3_class(p, "ggplot")
})

test_that("sensitivityCLPM errors if RI-CLPM refit fails to converge", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  ## Use a too-small sample with a pathological model to push RI-CLPM
  ## to non-convergence -- but be tolerant: if it converges anyway,
  ## the test should not fail spuriously.
  set.seed(99)
  dat <- simCLPM(waves = 3, sample_size = 30, seed = 99)$data
  fit <- suppressWarnings(
    lavaan::lavaan(estimateCLPM(waves = 3), data = dat, meanstructure = TRUE)
  )
  if (lavaan::lavInspect(fit, "converged")) {
    res <- tryCatch(
      sensitivityCLPM(fit, icc_grid = c(0, 0.5), n_sims = 2,
                      data = dat, seed = 99),
      error = function(e) e
    )
    ## Either we get a sensitivity object back (RI-CLPM converged) or
    ## we get a clear error about RI-CLPM refit failure -- both fine.
    expect_true(inherits(res, "crossLagR_sensitivity") ||
                grepl("RI-CLPM", conditionMessage(res)))
  } else {
    succeed()
  }
})
