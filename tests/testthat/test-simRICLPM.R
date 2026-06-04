test_that("simRICLPM returns expected structure with x*/y* columns", {
  skip_if_not_installed("lavaan")
  set.seed(2)
  out <- simRICLPM(waves = 4, sample.nobs = 300)
  expect_true(all(c("model", "data") %in% names(out)))
  expect_s3_class(out$data, "data.frame")
  expect_equal(nrow(out$data), 300)
  expect_true(all(c(paste0("x", 1:4), paste0("y", 1:4)) %in% names(out$data)))
})

test_that("simRICLPM include_z = TRUE adds z columns", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  set.seed(5)
  d <- simRICLPM(waves = 3, sample.nobs = 200, include_z = TRUE)$data
  expect_true(all(paste0("z", 1:3) %in% names(d)))
})

test_that("simRICLPM larger random-intercept variance increases ICC", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  ## cov_BXBY must respect |cov| <= sqrt(var_BX * var_BY) for the low-var run
  set.seed(17)
  low  <- simRICLPM(waves = 4, sample.nobs = 1000,
                    var_BX = 0.5, var_BY = 0.5, cov_BXBY = 0)$data
  high <- simRICLPM(waves = 4, sample.nobs = 1000,
                    var_BX = 3.0, var_BY = 3.0, cov_BXBY = 0)$data
  icc_low  <- stats::cor(low$x1,  low$x4)
  icc_high <- stats::cor(high$x1, high$x4)
  expect_gt(icc_high, icc_low)
})

test_that("generate_riclpm_indicator_lists_for_estimation returns expected list", {
  out <- generate_riclpm_indicator_lists_for_estimation(waves = 3, num_indicators = 2)
  expect_named(out, c("time_varying_x", "time_varying_y"))
  expect_length(out$time_varying_x, 2)
  expect_equal(out$time_varying_x[[1]], c("x1_1", "x1_2", "x1_3"))
  expect_equal(out$time_varying_y[[2]], c("y2_1", "y2_2", "y2_3"))
})

test_that("generate_riclpm_indicator_lists_for_estimation with z adds z list", {
  out <- generate_riclpm_indicator_lists_for_estimation(
    waves = 2, num_indicators = 2, has_z = TRUE,
    x_prefix = "X", y_prefix = "Y", z_prefix = "Z", indicator_separator = "."
  )
  expect_true("time_varying_z" %in% names(out))
  expect_equal(out$time_varying_z[[1]], c("Z1.1", "Z1.2"))
  expect_equal(out$time_varying_x[[2]], c("X2.1", "X2.2"))
})

test_that("generate_riclpm_indicator_lists_for_estimation rejects bad input", {
  expect_error(generate_riclpm_indicator_lists_for_estimation(waves = 0, num_indicators = 2))
  expect_error(generate_riclpm_indicator_lists_for_estimation(waves = 2, num_indicators = 0))
  expect_error(generate_riclpm_indicator_lists_for_estimation(waves = 2, num_indicators = 2,
                                                              has_z = "yes"))
  expect_error(generate_riclpm_indicator_lists_for_estimation(waves = 2, num_indicators = 2,
                                                              x_prefix = c("a", "b")))
})
