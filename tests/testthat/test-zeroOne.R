test_that("zero.one rescales numeric vector to [0,1]", {
  x <- c(-5, 0, 5, 10)
  out <- zero.one(x)
  expect_equal(min(out), 0)
  expect_equal(max(out), 1)
  expect_equal(out, (x - min(x)) / (max(x) - min(x)))
})

test_that("zero.one handles NA values", {
  out <- zero.one(c(NA, 1, 2, 3))
  expect_true(is.na(out[1]))
  expect_equal(min(out, na.rm = TRUE), 0)
  expect_equal(max(out, na.rm = TRUE), 1)
})

test_that("zero.one returns NaN for a constant vector", {
  out <- zero.one(rep(2, 4))
  expect_true(all(is.nan(out)))
})
