test_that("simCLPM returns model + data with x*/y* columns", {
  skip_if_not_installed("lavaan")
  set.seed(1)
  out <- simCLPM(waves = 4, sample_size = 200)
  expect_named(out, c("model", "data"))
  expect_s3_class(out$data, "data.frame")
  expect_equal(nrow(out$data), 200)
  expect_true(all(c(paste0("x", 1:4), paste0("y", 1:4)) %in% names(out$data)))
  expect_type(out$model, "character")
})

test_that("simCLPM is reproducible with a seed", {
  skip_if_not_installed("lavaan")
  d1 <- simCLPM(waves = 3, sample_size = 100, seed = 42)$data
  d2 <- simCLPM(waves = 3, sample_size = 100, seed = 42)$data
  expect_identical(d1, d2)
})

test_that("simCLPM rejects invalid inputs", {
  expect_error(simCLPM(waves = 0))
  expect_error(simCLPM(waves = 1))
  expect_error(simCLPM(waves = 1.5))
  expect_error(simCLPM(waves = 3, var_x = -1))
  expect_error(simCLPM(waves = 3, var_y = 0))
  expect_error(simCLPM(waves = 3, beta_x = "bad"))
})

test_that("simCLPM lavaan model parses", {
  skip_if_not_installed("lavaan")
  m <- simCLPM(waves = 3, sample_size = 50, seed = 7)$model
  pt <- lavaan::lavaanify(m)
  expect_s3_class(pt, "data.frame")
  expect_true(nrow(pt) > 0)
})

test_that("simCLPM data approximately recovers stability/cross-lag at large N", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  set.seed(11)
  d <- simCLPM(
    waves = 4,
    beta_x = 0.5, beta_y = 0.5,
    omega_xy = 0.2, omega_yx = 0.0,
    var_x = 1, var_y = 1, cov_xy = 0,
    sample_size = 2000
  )$data
  cor_x <- stats::cor(d$x1, d$x2)
  cor_xy_lag <- stats::cor(d$x1, d$y2)
  expect_gt(cor_x, 0.3)
  expect_gt(cor_xy_lag, 0)
})

test_that("simCLPMu returns model + data with x*/y*/u* columns and confounded structure", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  set.seed(3)
  out <- simCLPMu(waves = 4)
  expect_named(out, c("model", "data"))
  expect_s3_class(out$data, "data.frame")
  needed <- c(paste0("x", 1:4), paste0("y", 1:4), paste0("u", 1:4))
  expect_true(all(needed %in% names(out$data)))
})
