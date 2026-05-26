test_that("estimateLChange returns valid lavaan syntax with unified labels", {
  syntax <- estimateLChange(waves = 4, variable_type = "bivariate")
  expect_type(syntax, "character")
  expect_match(syntax, "ar_x")
  expect_match(syntax, "ar_y")
  expect_match(syntax, "cl_xy")
  expect_match(syntax, "cl_yx")
  expect_match(syntax, "u_var_x")
  expect_match(syntax, "u_var_y")
})

test_that("estimateLChange univariate has no cross-lagged labels", {
  syntax <- estimateLChange(waves = 4, variable_type = "univariate")
  expect_match(syntax, "ar_x")
  expect_no_match(syntax, "cl_xy")
  expect_no_match(syntax, "ar_y")
})

test_that("estimateLChange with accumulating factor includes A labels", {
  syntax <- estimateLChange(
    waves = 4, variable_type = "bivariate",
    estimate_constant_change = TRUE
  )
  expect_match(syntax, "A_x")
  expect_match(syntax, "A_y")
})

test_that("estimateLChange rejects too few waves", {
  expect_error(estimateLChange(waves = 2))
})

test_that("estimateLChange parses in lavaan", {
  skip_if_not_installed("lavaan")
  syntax <- estimateLChange(waves = 4, variable_type = "bivariate")
  model <- lavaan::lavaanify(syntax)
  expect_s3_class(model, "data.frame")
  expect_true(nrow(model) > 0)
})
