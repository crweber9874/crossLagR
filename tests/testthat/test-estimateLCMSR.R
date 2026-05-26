test_that("estimateLCMSR returns valid lavaan syntax with unified labels", {
  syntax <- estimateLCMSR(waves = 4)
  expect_type(syntax, "character")
  expect_match(syntax, "ar_x")
  expect_match(syntax, "ar_y")
  expect_match(syntax, "cl_xy")
  expect_match(syntax, "cl_yx")
  expect_match(syntax, "I_x")
  expect_match(syntax, "I_y")
  expect_match(syntax, "S_x")
  expect_match(syntax, "S_y")
  expect_match(syntax, "d_var_x")
  expect_match(syntax, "d_var_y")
  expect_match(syntax, "d_cov_xy")
})

test_that("estimateLCMSR unconstrained labels are wave-specific", {
  syntax <- estimateLCMSR(waves = 4, constrain_beta = FALSE, constrain_omega = FALSE)
  expect_match(syntax, "ar_x2")
  expect_match(syntax, "ar_x3")
  expect_match(syntax, "cl_xy2")
})

test_that("estimateLCMSR rejects invalid inputs", {
  expect_error(estimateLCMSR(waves = 2))
  expect_error(estimateLCMSR(waves = 4, time_scores = c(0, 1)))
  expect_error(estimateLCMSR(waves = 3, estimate_quadratic = TRUE))
})

test_that("estimateLCMSR parses in lavaan", {
  skip_if_not_installed("lavaan")
  syntax <- estimateLCMSR(waves = 4)
  model <- lavaan::lavaanify(syntax)
  expect_s3_class(model, "data.frame")
  expect_true(nrow(model) > 0)
})

test_that("estimateLCMSR custom time scores work", {
  syntax <- estimateLCMSR(waves = 4, time_scores = c(0, 1, 3, 6))
  expect_match(syntax, "3\\*x3")
  expect_match(syntax, "6\\*x4")
})

test_that("estimateLCMSR growth factors load on observed variables", {
  syntax <- estimateLCMSR(waves = 4)
  expect_match(syntax, "I_x =~ 1\\*x1")
  expect_match(syntax, "I_y =~ 1\\*y1")
  expect_match(syntax, "S_x =~ 0\\*x1")
  expect_match(syntax, "S_y =~ 0\\*y1")
})

test_that("estimateLCMSR quadratic factors work", {
  skip_if_not_installed("lavaan")
  syntax <- estimateLCMSR(waves = 5, estimate_quadratic = TRUE)
  expect_match(syntax, "Q_x")
  expect_match(syntax, "Q_y")
  model <- lavaan::lavaanify(syntax)
  expect_s3_class(model, "data.frame")
})
