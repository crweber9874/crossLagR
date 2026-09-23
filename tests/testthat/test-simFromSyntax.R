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

test_that("simFromSyntax fixes the random-intercept variance when asked", {
  skip_on_cran()
  skip_if_not_installed("lavaan")

  icc_x <- function(d) {
    w <- paste0("x", 1:4); cv <- cov(d[, w])
    mean(cv[lower.tri(cv)]) / mean(diag(cv))
  }
  set.seed(101)
  d_lo <- simFromSyntax("RICLPM", waves = 4, sample_size = 3000,
                        ar_x = .3, ar_y = .2, cl_xy = .15, cl_yx = .1,
                        var_between_x = 0.2, var_between_y = 0.2, cov_between_xy = 0.06)
  set.seed(101)
  d_hi <- simFromSyntax("RICLPM", waves = 4, sample_size = 3000,
                        ar_x = .3, ar_y = .2, cl_xy = .15, cl_yx = .1,
                        var_between_x = 3.0, var_between_y = 3.0, cov_between_xy = 0.9)
  ## More between-person variance -> higher ICC.
  expect_gt(icc_x(d_hi), icc_x(d_lo))
  expect_gt(icc_x(d_hi), 0.6)
})

test_that("Bollen-Brand DGP induces dynamic-panel bias in the RI-CLPM", {
  skip_on_cran()
  skip_if_not_installed("lavaan")

  grid <- data.frame(
    stability_p = 0.30, stability_q = 0.20,
    cross_q = 0.15, cross_p = 0.10,
    variance_p = 0.5, variance_q = 0.5,
    variance_between_x = 2, variance_between_y = 2, cov_pq = 0
  )
  set.seed(202)
  ri <- run_mc_sims("RICLPM", data_generation = "BB", param_grid = grid,
                    trials = 6, waves = 4, sample_size = 1200, verbose = FALSE)
  bb <- run_mc_sims("BB", data_generation = "BB", param_grid = grid,
                    trials = 6, waves = 4, sample_size = 1200, verbose = FALSE)

  ## BB recovers ar_x (unified label present, near truth); RI-CLPM is biased up.
  expect_false(all(is.na(bb$ar_x)))
  expect_lt(abs(mean(bb$ar_x, na.rm = TRUE) - 0.30), 0.10)
  expect_gt(mean(ri$ar_x, na.rm = TRUE), 0.50)
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
