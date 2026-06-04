test_that("monteCarloLavaan runs ALT under a CLPM DGP and returns unified labels", {
  skip_on_cran()
  skip_if_not_installed("lavaan")

  grid <- data.frame(
    stability_p = 0.3, stability_q = 0.3,
    cross_p = 0.1, cross_q = 0.1,
    variance_p = 0.5, variance_q = 0.5,
    variance_between_x = 0.5, variance_between_y = 0.5,
    cov_pq = 0
  )

  set.seed(42)
  out <- monteCarloLavaan(
    estimator   = "ALT",
    dgp         = "clpm",
    trials      = 2,
    waves       = 5,
    sample_size = 300,
    param_grid  = grid,
    verbose     = FALSE
  )

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 2)
  for (col in c("ar_x", "ar_y", "cl_xy", "cl_yx",
                "se_ar_x", "se_ar_y", "se_cl_xy", "se_cl_yx",
                "std_ar_x", "std_ar_y", "std_cl_xy", "std_cl_yx",
                "n_obs", "n_miss", "converged",
                "chisq", "df", "cfi", "tli", "rmsea", "srmr",
                "true_ar_x", "true_ar_y", "true_cl_xy", "true_cl_yx",
                "estimator", "dgp", "trial")) {
    expect_true(col %in% names(out), info = paste("missing column:", col))
  }
  expect_true(all(out$estimator == "ALT"))
  expect_true(all(out$dgp == "CLPM"))
})

test_that("monteCarloLavaan runs CLPM under a CLPM DGP and recovers truth approximately", {
  skip_on_cran()
  skip_if_not_installed("lavaan")

  grid <- data.frame(
    stability_p = 0.4, stability_q = 0.4,
    cross_p = 0.15, cross_q = 0.15,
    variance_p = 0.5, variance_q = 0.5,
    variance_between_x = 0.5, variance_between_y = 0.5,
    cov_pq = 0
  )

  set.seed(7)
  out <- monteCarloLavaan(
    estimator   = "CLPM",
    dgp         = "clpm",
    trials      = 5,
    waves       = 5,
    sample_size = 1000,
    param_grid  = grid
  )

  expect_equal(nrow(out), 5)
  conv <- out[out$converged %in% TRUE, , drop = FALSE]
  if (nrow(conv) > 0) {
    expect_lt(abs(mean(conv$ar_x, na.rm = TRUE) - 0.4), 0.15)
    expect_lt(abs(mean(conv$cl_xy, na.rm = TRUE) - 0.15), 0.15)
  }
})

test_that("monteCarloLavaan errors on unknown estimator / dgp", {
  expect_error(monteCarloLavaan("BOGUS"), "estimator must be one of")
  expect_error(monteCarloLavaan("CLPM", dgp = "bogus"), "dgp must be one of")
})

test_that("run_mc_sims dispatches ALT through the generic wrapper", {
  skip_on_cran()
  skip_if_not_installed("lavaan")

  grid <- data.frame(
    stability_p = 0.3, stability_q = 0.3,
    cross_p = 0.1, cross_q = 0.1,
    variance_p = 0.5, variance_q = 0.5,
    variance_between_x = 0.5, variance_between_y = 0.5,
    cov_pq = 0
  )

  set.seed(11)
  out <- run_mc_sims(
    estimator       = "ALT",
    param_grid      = grid,
    trials          = 1,
    waves           = 4,
    sample_size     = 300,
    verbose         = FALSE,
    data_generation = "clpm"
  )

  expect_s3_class(out, "data.frame")
  expect_true("ar_x" %in% names(out))
  expect_true(any(out$estimator == "ALT"))
})
