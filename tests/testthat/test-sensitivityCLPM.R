## Unit + integration tests for sensitivityCLPM.

fit_clpm_for_sens <- function() {
  set.seed(7001)
  dat <- simCLPM(waves = 4, sample_size = 300, seed = 7001,
                 beta_x = 0.4, beta_y = 0.4,
                 omega_xy = 0.15, omega_yx = 0.10)$data
  lavaan::lavaan(estimateCLPM(waves = 4), data = dat, meanstructure = TRUE)
}

test_that("sensitivityCLPM returns expected structure", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  fit <- fit_clpm_for_sens()
  s <- sensitivityCLPM(fit, icc_grid = c(0, 0.3, 0.6), n_sims = 3, seed = 1)
  expect_s3_class(s, "crossLagR_sensitivity")
  expect_named(s, c("estimates", "original", "icc_grid",
                    "n_sims", "sample_size", "waves"))
  expect_s3_class(s$estimates, "data.frame")
  expect_true(all(c("icc", "sim", "estimator", "parameter", "est",
                    "original_est") %in% names(s$estimates)))
  expect_setequal(unique(s$estimates$parameter),
                  c("ar_x", "ar_y", "cl_xy", "cl_yx"))
  expect_setequal(unique(s$estimates$estimator), c("CLPM", "RI-CLPM"))
})

test_that("sensitivityCLPM original matches CLPM fit estimates", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  fit <- fit_clpm_for_sens()
  s <- sensitivityCLPM(fit, icc_grid = c(0), n_sims = 2, seed = 2)
  expect_equal(s$original[["ar_x"]],
               lavaan::parameterEstimates(fit)$est[
                 lavaan::parameterEstimates(fit)$label == "ar_x"][1])
})

test_that("sensitivityCLPM rejects bad input", {
  expect_error(sensitivityCLPM(list()), "lavaan")
  expect_error(sensitivityCLPM(1), "lavaan")
})

test_that("sensitivityCLPM rejects icc_grid out of [0, 1)", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  fit <- fit_clpm_for_sens()
  expect_error(sensitivityCLPM(fit, icc_grid = c(-0.1)),  "icc_grid")
  expect_error(sensitivityCLPM(fit, icc_grid = c(1.0)),   "icc_grid")
  expect_error(sensitivityCLPM(fit, icc_grid = c(1.5)),   "icc_grid")
})

test_that("sensitivityCLPM ICC = 0 case is near-baseline (within MC noise)", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  fit <- fit_clpm_for_sens()
  s <- sensitivityCLPM(fit, icc_grid = 0, n_sims = 8, seed = 3,
                       sample_size = 500)
  clpm0 <- s$estimates[s$estimates$estimator == "CLPM" &
                       s$estimates$parameter == "ar_x", ]
  expect_lt(abs(mean(clpm0$est, na.rm = TRUE) - s$original[["ar_x"]]), 0.10)
})

test_that("sensitivityCLPM CLPM AR_x increases monotonically with injected ICC", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  fit <- fit_clpm_for_sens()
  s <- sensitivityCLPM(fit, icc_grid = c(0, 0.3, 0.6), n_sims = 6,
                       seed = 11, sample_size = 600)
  agg <- stats::aggregate(
    est ~ icc, data = s$estimates[s$estimates$estimator == "CLPM" &
                                  s$estimates$parameter == "ar_x", ],
    FUN = mean, na.rm = TRUE
  )
  ## Hamaker-style bias: CLPM AR inflates with ICC.
  expect_gt(agg$est[agg$icc == 0.6], agg$est[agg$icc == 0])
})

test_that("sensitivityCLPM RI-CLPM stays closer to truth than CLPM at high ICC", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  fit <- fit_clpm_for_sens()
  s <- sensitivityCLPM(fit, icc_grid = 0.6, n_sims = 8, seed = 13,
                       sample_size = 600)
  truth <- s$original[["ar_x"]]
  ar_x  <- s$estimates[s$estimates$parameter == "ar_x", ]
  clpm_err   <- abs(mean(ar_x$est[ar_x$estimator == "CLPM"],    na.rm = TRUE) - truth)
  riclpm_err <- abs(mean(ar_x$est[ar_x$estimator == "RI-CLPM"], na.rm = TRUE) - truth)
  expect_lt(riclpm_err, clpm_err)
})

test_that("print.crossLagR_sensitivity runs without error", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  fit <- fit_clpm_for_sens()
  s <- sensitivityCLPM(fit, icc_grid = c(0, 0.3), n_sims = 2, seed = 4)
  out <- capture.output(print(s))
  expect_true(any(grepl("sensitivity analysis", out, ignore.case = TRUE)))
  expect_true(any(grepl("plot\\(\\) ", out)))
})

test_that("plot.crossLagR_sensitivity returns a ggplot", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("ggplot2")
  fit <- fit_clpm_for_sens()
  s <- sensitivityCLPM(fit, icc_grid = c(0, 0.4), n_sims = 2, seed = 5)
  p <- plot(s)
  expect_s3_class(p, "ggplot")
})
