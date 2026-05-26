test_that("estimateALT returns valid lavaan syntax with unified labels", {
  syntax <- estimateALT(waves = 4)
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

test_that("estimateALT unconstrained labels are wave-specific", {
  syntax <- estimateALT(waves = 4, constrain_beta = FALSE, constrain_omega = FALSE)
  expect_match(syntax, "ar_x2")
  expect_match(syntax, "ar_x3")
  expect_match(syntax, "cl_xy2")
})

test_that("estimateALT rejects invalid inputs", {
  expect_error(estimateALT(waves = 2))
  expect_error(estimateALT(waves = 4, time_scores = c(0, 1)))
})

test_that("estimateALT parses in lavaan", {
  skip_if_not_installed("lavaan")
  syntax <- estimateALT(waves = 4)
  model <- lavaan::lavaanify(syntax)
  expect_s3_class(model, "data.frame")
  expect_true(nrow(model) > 0)
})

test_that("estimateALT custom time scores work", {
  syntax <- estimateALT(waves = 4, time_scores = c(0, 1, 3, 6))
  expect_match(syntax, "3\\*p3")
  expect_match(syntax, "6\\*p4")
})
