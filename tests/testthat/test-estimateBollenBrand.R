test_that("estimateBollen_and_Brand returns valid lavaan syntax with unified labels", {
  syntax <- estimateBollen_and_Brand(waves = 4)
  expect_type(syntax, "character")
  expect_match(syntax, "ar_y")
  expect_match(syntax, "cl_xy")
  expect_match(syntax, "I_y")
})

test_that("estimateBollen_and_Brand reciprocal model includes X labels", {
  syntax <- estimateBollen_and_Brand(
    waves = 4, x_effect = "lagged",
    x_autoregression = TRUE, y_effect_on_x = TRUE
  )
  expect_match(syntax, "ar_x")
  expect_match(syntax, "cl_yx")
  expect_match(syntax, "I_x")
})

test_that("estimateBollen_and_Brand constrained uses single labels", {
  syntax <- estimateBollen_and_Brand(
    waves = 4, x_effect = "lagged",
    x_autoregression = TRUE, y_effect_on_x = TRUE,
    constrain_coefficients = TRUE
  )
  expect_match(syntax, "ar_y\\*y")
  expect_match(syntax, "ar_x\\*x")
  expect_match(syntax, "cl_xy\\*x")
  expect_match(syntax, "cl_yx\\*y")
})

test_that("estimateBollen_and_Brand rejects invalid inputs", {
  expect_error(estimateBollen_and_Brand(waves = 2))
  expect_error(estimateBollen_and_Brand(waves = 4, x_effect = "bad"))
  expect_error(estimateBollen_and_Brand(
    waves = 4, y_effect_on_x = TRUE, x_autoregression = FALSE
  ))
})

test_that("estimateBollen_and_Brand parses in lavaan", {
  skip_if_not_installed("lavaan")
  syntax <- estimateBollen_and_Brand(
    waves = 4, x_effect = "lagged",
    x_autoregression = TRUE, y_effect_on_x = TRUE
  )
  model <- lavaan::lavaanify(syntax)
  expect_s3_class(model, "data.frame")
  expect_true(nrow(model) > 0)
})
