test_that("estimateLGM univariate returns valid syntax with growth labels", {
  syntax <- estimateLGM(waves = 4)
  expect_type(syntax, "character")
  expect_match(syntax, "I =~ 1\\*p1")
  expect_match(syntax, "S =~ 0\\*p1")
  expect_match(syntax, "mean_I")
  expect_match(syntax, "mean_S")
  expect_match(syntax, "var_I")
  expect_match(syntax, "var_S")
})

test_that("estimateLGM bivariate adds I_x/S_x/I_y/S_y labels", {
  syntax <- estimateLGM(waves = 4, variable_type = "bivariate")
  for (lab in c("I_x", "S_x", "I_y", "S_y",
                "mean_I_x", "mean_I_y",
                "var_I_x", "var_I_y",
                "var_S_x", "var_S_y")) {
    expect_match(syntax, lab, fixed = TRUE)
  }
})

test_that("estimateLGM rejects invalid inputs", {
  expect_error(estimateLGM(waves = 2))
  expect_error(estimateLGM(waves = 0))
  expect_error(estimateLGM(waves = 4, time_scores = c(0, 1, 2)))
  expect_error(estimateLGM(waves = 3, estimate_quadratic = TRUE))
  expect_error(estimateLGM(waves = 4, center_time = "yes"))
})

test_that("estimateLGM custom time scores embed in syntax", {
  syntax <- estimateLGM(waves = 4, time_scores = c(0, 1, 3, 7))
  expect_match(syntax, "3\\*p3")
  expect_match(syntax, "7\\*p4")
})

test_that("estimateLGM quadratic factor appears with Q label", {
  syntax <- estimateLGM(waves = 5, estimate_quadratic = TRUE)
  expect_match(syntax, "Q =~")
  expect_match(syntax, "mean_Q")
  expect_match(syntax, "var_Q")
})

test_that("estimateLGM parses in lavaan (univariate + bivariate)", {
  skip_if_not_installed("lavaan")
  for (vt in c("univariate", "bivariate")) {
    syntax <- estimateLGM(waves = 4, variable_type = vt)
    pt <- lavaan::lavaanify(syntax)
    expect_s3_class(pt, "data.frame")
    expect_true(nrow(pt) > 0)
  }
})

test_that("estimateLGM center_time produces symmetric time scores", {
  syntax <- estimateLGM(waves = 5, center_time = TRUE)
  expect_match(syntax, "-2\\*p1")
  expect_match(syntax, "2\\*p5")
})
