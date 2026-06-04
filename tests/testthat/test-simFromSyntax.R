test_that("simFromSyntax returns bivariate data for every unified estimator", {
  skip_on_cran()
  skip_if_not_installed("lavaan")

  for (est in c("CLPM", "RICLPM", "ALT", "LGM", "LCMSR",
                "LCHANGE", "BB", "TSO")) {
    set.seed(13)
    dat <- simFromSyntax(estimator = est, waves = 4, sample_size = 200,
                         ar_x = 0.3, ar_y = 0.3, cl_xy = 0.2, cl_yx = 0.1)
    expect_s3_class(dat, "data.frame")
    expect_equal(nrow(dat), 200)
    needed <- c(paste0("x", 1:4), paste0("y", 1:4))
    expect_true(all(needed %in% names(dat)),
                label = paste(est, "missing:",
                              paste(setdiff(needed, names(dat)), collapse = ",")))
  }
})

test_that("populate_unified_labels substitutes ar/cl values", {
  s <- "p2 ~ ar_x*p1 + cl_yx*q1\np2 ~~ d_var_x*p2\n"
  pop <- populate_unified_labels(s, ar_x = 0.42, cl_yx = -0.13, d_var_x = 0.9)
  expect_true(grepl("0.42\\*p1", pop))
  expect_true(grepl("-0.13\\*q1", pop))
  expect_true(grepl("0.9\\*p2", pop))
  expect_false(grepl("ar_x\\*", pop))
})

test_that("run_mc_sims accepts estimator-as-DGP combinations", {
  skip_on_cran()
  skip_if_not_installed("lavaan")

  grid <- data.frame(
    stability_p = 0.3, stability_q = 0.3,
    cross_p = 0.1, cross_q = 0.2,
    variance_p = 0.5, variance_q = 0.5,
    variance_between_x = 0.5, variance_between_y = 0.5,
    cov_pq = 0
  )

  set.seed(99)
  res <- run_mc_sims(
    estimator       = "CLPM",
    data_generation = "ALT",   ## estimator-as-DGP route
    param_grid      = grid,
    trials          = 2,
    waves           = 4,
    sample_size     = 300,
    verbose         = FALSE
  )

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  for (col in c("ar_x", "ar_y", "cl_xy", "cl_yx", "cfi", "converged",
                "n_obs", "true_ar_x")) {
    expect_true(col %in% names(res), info = paste("missing col:", col))
  }
  expect_true(all(res$dgp == "ALT"))
  expect_true(all(res$estimator == "CLPM"))
})
