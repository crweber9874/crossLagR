## Unit tests for diagnoseFit / lavaanDiagnose / heywoodHelp.

clean_fit <- function() {
  set.seed(401)
  dat <- simCLPM(waves = 4, sample_size = 500, seed = 401)$data
  lavaan::lavaan(estimateCLPM(waves = 4), data = dat, meanstructure = TRUE)
}

test_that("diagnoseFit returns 0 issues for a clean fit", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  fit <- clean_fit()
  diag <- diagnoseFit(fit, quiet = TRUE)
  expect_s3_class(diag, "crossLagR_diagnosis")
  expect_equal(diag$n_issues, 0L)
  expect_match(diag$summary, "No Heywood")
})

test_that("diagnoseFit errors on non-lavaan input", {
  expect_error(diagnoseFit(list()), "lavaan object")
  expect_error(diagnoseFit(1), "lavaan object")
})

test_that("diagnoseFit picks up NPD-covariance warning text", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  fit <- clean_fit()
  w <- c("lavaan WARNING: the covariance matrix is not positive definite")
  diag <- diagnoseFit(fit, warnings_seen = w, quiet = TRUE)
  expect_true("npd_covariance" %in% diag$codes)
})

test_that("diagnoseFit picks up vcov-not-PD warning text", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  fit <- clean_fit()
  w <- c("the variance-covariance matrix of the estimated parameters (vcov) does not appear positive definite")
  diag <- diagnoseFit(fit, warnings_seen = w, quiet = TRUE)
  expect_true("vcov_not_pd" %in% diag$codes)
})

test_that("diagnoseFit picks up iter-limit warning text", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  fit <- clean_fit()
  w <- c("optimizer reached the iteration limit; the model did not converge")
  diag <- diagnoseFit(fit, warnings_seen = w, quiet = TRUE)
  expect_true("iter_limit" %in% diag$codes)
})

test_that("diagnoseFit emits a single warning by default when issues exist", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  fit <- clean_fit()
  expect_warning(
    diagnoseFit(fit, warnings_seen = "covariance matrix is not positive definite"),
    "Heywood-class"
  )
})

test_that("diagnoseFit is silent with quiet = TRUE", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  fit <- clean_fit()
  expect_no_warning(
    diagnoseFit(fit, warnings_seen = "not positive definite", quiet = TRUE)
  )
})

test_that("diagnoseFit avoids duplicate npd entry when negative_variance fires", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  ## Fit RICLPM to data with zero between-person variance -> often Heywood.
  set.seed(909)
  dat <- simCLPM(waves = 4, sample_size = 200,
                 beta_x = 0.3, beta_y = 0.3,
                 omega_xy = 0.1, omega_yx = 0.1,
                 seed = 909)$data
  fit <- suppressWarnings(
    lavaan::lavaan(estimateRICLPM(waves = 4), data = dat, meanstructure = TRUE)
  )
  diag <- diagnoseFit(
    fit,
    warnings_seen = c("not positive definite", "lavaan warning"),
    quiet = TRUE
  )
  ## If we caught a negative variance, npd_covariance should be suppressed.
  if ("negative_variance" %in% diag$codes) {
    expect_false("npd_covariance" %in% diag$codes)
  } else {
    succeed()
  }
})

test_that("print.crossLagR_diagnosis runs and prints expected sections", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  fit <- clean_fit()
  diag <- diagnoseFit(fit, warnings_seen = "not positive definite",
                      quiet = TRUE)
  out <- capture.output(print(diag))
  expect_true(any(grepl("issue", out, ignore.case = TRUE)))
  expect_true(any(grepl("Remedies", out)))
  expect_true(any(grepl("blavaan", out)))
})

test_that("summary.crossLagR_diagnosis prints one-line summary", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  fit <- clean_fit()
  diag <- diagnoseFit(fit, quiet = TRUE)
  out <- capture.output(summary(diag))
  expect_match(out, "No Heywood|issue", all = FALSE)
})

test_that("heywoodHelp prints all catalog entries", {
  out <- capture.output(heywoodHelp())
  expected <- c("non_convergence", "negative_variance", "std_loading_gt_one",
                "std_correlation_gt_one", "npd_covariance", "vcov_not_pd",
                "iter_limit")
  for (code in expected) {
    expect_true(any(grepl(code, out)),
                info = paste("missing code in heywoodHelp() output:", code))
  }
})

test_that("heywoodHelp can filter to a single code", {
  out <- capture.output(heywoodHelp("negative_variance"))
  expect_true(any(grepl("negative_variance", out)))
  expect_false(any(grepl("vcov_not_pd", out)))
})

test_that("heywoodHelp rejects unknown codes", {
  expect_error(heywoodHelp("bogus_code"), "Unknown code")
})

test_that("lavaanDiagnose returns a crossLagR_fit bundle", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  set.seed(11)
  dat <- simCLPM(waves = 4, sample_size = 300, seed = 11)$data
  out <- suppressWarnings(
    lavaanDiagnose(estimateCLPM(waves = 4), dat, meanstructure = TRUE)
  )
  expect_s3_class(out, "crossLagR_fit")
  expect_named(out, c("fit", "summary", "diagnosis", "icc", "warnings"))
  expect_s4_class(out$fit, "lavaan")
  expect_s3_class(out$summary, "crossLagR_summary")
  expect_s3_class(out$diagnosis, "crossLagR_diagnosis")
  expect_null(out$icc)
})

test_that("lavaanDiagnose with icc_vars attaches an ICC report", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  set.seed(12)
  dat <- simCLPM(waves = 4, sample_size = 300, seed = 12)$data
  out <- suppressWarnings(
    lavaanDiagnose(estimateCLPM(waves = 4), dat,
                   meanstructure = TRUE, icc_vars = c("x", "y"))
  )
  expect_s3_class(out$icc, "crossLagR_icc")
  expect_equal(nrow(out$icc), 2L)
})

test_that("print.crossLagR_fit runs without error", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  set.seed(13)
  dat <- simCLPM(waves = 4, sample_size = 200, seed = 13)$data
  out <- suppressWarnings(
    lavaanDiagnose(estimateCLPM(waves = 4), dat,
                   meanstructure = TRUE, icc_vars = c("x", "y"))
  )
  printed <- capture.output(print(out))
  expect_true(any(grepl("crossLagR model summary", printed)))
  expect_true(any(grepl("Fit statistics", printed)))
})

test_that("monteCarloLavaan adds heywood_codes column and condensed warning", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  grid <- data.frame(
    stability_p = 0.3, stability_q = 0.3,
    cross_p = 0.1, cross_q = 0.1,
    variance_p = 0.5, variance_q = 0.5,
    variance_between_x = 0.5, variance_between_y = 0.5,
    cov_pq = 0
  )
  set.seed(55)
  out <- suppressWarnings(monteCarloLavaan(
    estimator = "CLPM", dgp = "CLPM",
    trials = 2, waves = 4, sample_size = 250,
    param_grid = grid, verbose = FALSE
  ))
  expect_true("heywood_n_issues" %in% names(out))
  expect_true("heywood_codes" %in% names(out))
  expect_type(out$heywood_n_issues, "integer")
})
