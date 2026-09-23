test_that("clpmBandStyle returns every styling key the plot needs", {
  st <- clpmBandStyle()
  expect_type(st, "list")
  expect_true(all(c("band_fill", "band_alpha", "line_colour", "ref_linetype",
                    "base_size", "ar_label", "cl_label") %in% names(st)))
  expect_equal(clpmBandStyle(band_fill = "steelblue")$band_fill, "steelblue")
})

test_that("clpmSensitivityBands returns a plot carrying its simulation", {
  skip_on_cran()
  p <- clpmSensitivityBands(0.55, 0.12, r_grid = c(0.5, 1, 2),
                            reps = 3, n_obs = 300, seed = 1)
  expect_s3_class(p, "ggplot")
  sim <- attr(p, "sim")
  expect_s3_class(sim, "data.frame")
  expect_equal(nrow(sim), 3L)
  expect_true(all(c("bw_ratio", "AR", "AR_lo", "AR_hi",
                    "CL", "CL_lo", "CL_hi", "n_obs", "n_miss") %in% names(sim)))
  expect_true(all(sim$AR_lo <= sim$AR & sim$AR <= sim$AR_hi))
})

test_that("supplying `sim` skips the simulation entirely", {
  skip_on_cran()
  sim <- clpmSensitivityBands(0.55, 0.12, r_grid = c(0.5, 1),
                              reps = 3, n_obs = 300, seed = 1, plot = FALSE)
  # no lavaan work should happen here; a bad grid would error if it did
  p <- clpmSensitivityBands(0.55, 0.12, sim = sim, r_grid = numeric(0))
  expect_s3_class(p, "ggplot")
  expect_equal(attr(p, "sim"), sim)
})

test_that("plot = FALSE returns the swept data frame", {
  skip_on_cran()
  out <- clpmSensitivityBands(0.55, 0.12, r_grid = c(1, 2),
                              reps = 2, n_obs = 200, seed = 2, plot = FALSE)
  expect_s3_class(out, "data.frame")
  expect_false(inherits(out, "ggplot"))
})

test_that("clpmSensitivityBands validates its inputs", {
  expect_error(clpmSensitivityBands(c(0.5, 0.6), 0.12), "length")
  expect_error(clpmSensitivityBands(0.55, 0.12, conf = 1.5), "conf")
})
