test_that("estimateAllisonChamberlainFI returns valid lavaan syntax with unified labels", {
  syntax <- estimateAllisonChamberlainFI(waves = 4)
  expect_type(syntax, "character")
  expect_match(syntax, "ar_x")
  expect_match(syntax, "ar_y")
  expect_match(syntax, "cl_xy")
  expect_match(syntax, "cl_yx")
  expect_match(syntax, "I_x")
  expect_match(syntax, "I_y")
  expect_match(syntax, "d_var_x")
  expect_match(syntax, "d_var_y")
})

test_that("estimateAllisonChamberlainFI rejects too few waves", {
  expect_error(estimateAllisonChamberlainFI(waves = 2))
})

test_that("estimateAllisonChamberlainFI parses in lavaan", {
  skip_if_not_installed("lavaan")
  syntax <- estimateAllisonChamberlainFI(waves = 4)
  model <- lavaan::lavaanify(syntax)
  expect_s3_class(model, "data.frame")
  expect_true(nrow(model) > 0)
})
