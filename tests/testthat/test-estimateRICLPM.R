test_that("estimateRICLPM returns valid lavaan syntax with unified labels", {
  syntax <- estimateRICLPM(waves = 4)
  expect_type(syntax, "character")
  expect_match(syntax, "ar_x")
  expect_match(syntax, "ar_y")
  expect_match(syntax, "cl_xy")
  expect_match(syntax, "cl_yx")
  expect_match(syntax, "I_x")
  expect_match(syntax, "I_y")
})

test_that("estimateRICLPM unconstrained labels are wave-specific", {
  syntax <- estimateRICLPM(waves = 3, constrain_beta = FALSE, constrain_omega = FALSE)
  expect_match(syntax, "ar_x2")
  expect_match(syntax, "ar_x3")
  expect_match(syntax, "cl_xy2")
})

test_that("estimateRICLPM with third variable includes Z labels", {
  syntax <- estimateRICLPM(waves = 4, include_z = TRUE)
  expect_match(syntax, "ar_z")
  expect_match(syntax, "I_z")
  expect_match(syntax, "cl_xz")
  expect_match(syntax, "cl_yz")
})

test_that("estimateRICLPM parses in lavaan", {
  skip_if_not_installed("lavaan")
  syntax <- estimateRICLPM(waves = 4)
  model <- lavaan::lavaanify(syntax)
  expect_s3_class(model, "data.frame")
  expect_true(nrow(model) > 0)
})

test_that("estimateRICLPM fits simulated data", {
  skip_if_not_installed("lavaan")
  set.seed(456)
  dat <- simRICLPM(
    waves = 4,
    beta_x = 0.3, beta_y = 0.3,
    omega_xy = 0.1, omega_yx = 0.1,
    var_p = 1, var_q = 1,
    var_BX = 1.5, var_BY = 1.5,
    sample.nobs = 500
  )$data
  syntax <- estimateRICLPM(waves = 4)
  fit <- suppressWarnings(lavaan::lavaan(syntax, data = dat))
  expect_true(lavaan::lavInspect(fit, "converged"))

  pt <- lavaan::parameterEstimates(fit)
  # Constrained labels appear multiple times; check label exists and estimate is finite
  ar_x_rows <- pt[pt$label == "ar_x", ]
  cl_xy_rows <- pt[pt$label == "cl_xy", ]
  expect_true(nrow(ar_x_rows) > 0)
  expect_true(nrow(cl_xy_rows) > 0)
  expect_true(is.finite(ar_x_rows$est[1]))
  expect_true(is.finite(cl_xy_rows$est[1]))
})
