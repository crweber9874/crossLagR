test_that("estimateCLPM returns valid lavaan syntax with unified labels", {
  syntax <- estimateCLPM(waves = 4)
  expect_type(syntax, "character")
  expect_match(syntax, "ar_y")
  expect_match(syntax, "ar_x")
  expect_match(syntax, "cl_xy")
  expect_match(syntax, "cl_yx")
  expect_match(syntax, "d_var_y")
  expect_match(syntax, "d_var_x")
  expect_match(syntax, "d_cov_xy")
})

test_that("estimateCLPM unconstrained labels are wave-specific", {
  syntax <- estimateCLPM(waves = 3, constrain_beta = FALSE, constrain_omega = FALSE)
  expect_match(syntax, "ar_y2")
  expect_match(syntax, "ar_y3")
  expect_match(syntax, "cl_xy2")
  expect_match(syntax, "cl_xy3")
})

test_that("estimateCLPM rejects invalid inputs", {
  expect_error(estimateCLPM(waves = 0))
  expect_error(estimateCLPM(waves = 1))
  expect_error(estimateCLPM(waves = -1))
  expect_error(estimateCLPM(waves = 3, constrain_beta = "yes"))
})

test_that("estimateCLPM parses in lavaan", {
  skip_if_not_installed("lavaan")
  syntax <- estimateCLPM(waves = 3)
  model <- lavaan::lavaanify(syntax)
  expect_s3_class(model, "data.frame")
  expect_true(nrow(model) > 0)
})

test_that("estimateCLPM fits simulated data", {
  skip_if_not_installed("lavaan")
  set.seed(123)
  dat <- simCLPM(
    waves = 4,
    beta_x = 0.3, beta_y = 0.3,
    omega_xy = 0.1, omega_yx = 0.1,
    var_x = 1, var_y = 1,
    sample_size = 500
  )$data
  syntax <- estimateCLPM(waves = 4)
  fit <- lavaan::lavaan(syntax, data = dat)
  expect_true(lavaan::lavInspect(fit, "converged"))

  pt <- lavaan::parameterEstimates(fit)
  ar_x_rows <- pt[pt$label == "ar_x", ]
  expect_true(nrow(ar_x_rows) > 0)
  expect_true(is.finite(ar_x_rows$est[1]))
})
