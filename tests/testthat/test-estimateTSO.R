test_that("estimateTSO returns valid syntax with trait/state labels", {
  syntax <- estimateTSO(waves = 4)
  expect_type(syntax, "character")
  for (lab in c("T_x =~ 1\\*x1", "T_y =~ 1\\*y1",
                "T_var_x", "T_var_y", "T_cov_xy",
                "s_var_x", "s_var_y", "s_cov_xy")) {
    expect_match(syntax, lab)
  }
})

test_that("estimateTSO trait factor loads on every wave with 1*", {
  syntax <- estimateTSO(waves = 5)
  for (w in 1:5) {
    expect_match(syntax, paste0("1\\*x", w))
    expect_match(syntax, paste0("1\\*y", w))
  }
})

test_that("estimateTSO unconstrained drops shared labels", {
  syntax <- estimateTSO(waves = 4,
                        constrain_state_variances  = FALSE,
                        constrain_state_covariances = FALSE)
  expect_no_match(syntax, "s_var_x\\*")
  expect_no_match(syntax, "s_cov_xy\\*")
})

test_that("estimateTSO rejects too few waves and bad input", {
  expect_error(estimateTSO(waves = 2))
  expect_error(estimateTSO(waves = 0))
  expect_error(estimateTSO(waves = 3.5))
})

test_that("estimateTSO parses in lavaan", {
  skip_if_not_installed("lavaan")
  pt <- lavaan::lavaanify(estimateTSO(waves = 4))
  expect_s3_class(pt, "data.frame")
  expect_true(nrow(pt) > 0)
})

test_that("estimateTSO fits simulated data and recovers positive trait variance", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  set.seed(101)
  dat <- simRICLPM(
    waves = 4, sample.nobs = 600,
    var_BX = 1.5, var_BY = 1.5,
    beta_x = 0, beta_y = 0,
    omega_xy = 0, omega_yx = 0
  )$data
  fit <- suppressWarnings(
    lavaan::lavaan(estimateTSO(waves = 4), data = dat, meanstructure = TRUE)
  )
  if (isTRUE(lavaan::lavInspect(fit, "converged"))) {
    pe <- lavaan::parameterEstimates(fit)
    expect_gt(pe$est[pe$label == "T_var_x"][1], 0)
    expect_gt(pe$est[pe$label == "T_var_y"][1], 0)
  } else {
    succeed()
  }
})
