test_that("populate_unified_labels substitutes every unified label", {
  s <- paste(
    "p2 ~ ar_x*p1 + cl_yx*q1",
    "q2 ~ ar_y*q1 + cl_xy*p1",
    "p2 ~~ d_var_x*p2",
    "q2 ~~ d_var_y*q2",
    "p2 ~~ d_cov_xy*q2",
    sep = "\n"
  )
  out <- populate_unified_labels(s,
    ar_x = 0.31, ar_y = 0.32, cl_xy = 0.11, cl_yx = 0.12,
    d_var_x = 0.41, d_var_y = 0.42, d_cov_xy = 0.05
  )
  for (val in c("0.31\\*p1", "0.32\\*q1", "0.11\\*p1", "0.12\\*q1",
                "0.41\\*p2", "0.42\\*q2", "0.05\\*q2")) {
    expect_match(out, val)
  }
  for (lab in c("ar_x\\*", "ar_y\\*", "cl_xy\\*", "cl_yx\\*",
                "d_var_x\\*", "d_var_y\\*", "d_cov_xy\\*")) {
    expect_no_match(out, lab)
  }
})

test_that("populate_unified_labels uses word boundaries (no partial matches)", {
  s <- "p2 ~ ar_xx*p1"
  out <- populate_unified_labels(s, ar_x = 0.5)
  expect_equal(out, s)
})

test_that("populate_unified_labels handles negative and zero values", {
  s <- "p2 ~ ar_x*p1 + cl_xy*q1\n p2 ~~ d_cov_xy*q2"
  out <- populate_unified_labels(s, ar_x = -0.25, cl_xy = 0, d_cov_xy = -0.1)
  expect_match(out, "-0.25\\*p1")
  expect_match(out, "0\\*q1")
  expect_match(out, "-0.1\\*q2")
})

test_that("populate_unified_labels is a no-op for syntax without unified labels", {
  s <- "y1 ~ x1 + z\n y1 ~~ 1*y1\n"
  expect_identical(populate_unified_labels(s), s)
})

test_that("populated syntax simulates with lavaan", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  syntax <- estimateCLPM(waves = 3)
  populated <- populate_unified_labels(syntax,
    ar_x = 0.3, ar_y = 0.3, cl_xy = 0.1, cl_yx = 0.1,
    d_var_x = 0.5, d_var_y = 0.5, d_cov_xy = 0
  )
  set.seed(3)
  dat <- lavaan::simulateData(model = populated, sample.nobs = 200)
  expect_s3_class(dat, "data.frame")
  expect_equal(nrow(dat), 200)
})
