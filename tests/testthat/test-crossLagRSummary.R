## Unit tests for crossLagRSummary + unifiedLabels.

fit_clean_clpm <- function() {
  set.seed(801)
  dat <- simCLPM(waves = 4, sample_size = 500, seed = 801)$data
  lavaan::lavaan(estimateCLPM(waves = 4), data = dat, meanstructure = TRUE)
}

fit_clean_riclpm <- function() {
  set.seed(802)
  dat <- simRICLPM(waves = 4, sample.nobs = 500)$data
  suppressWarnings(
    lavaan::lavaan(estimateRICLPM(waves = 4), data = dat, meanstructure = TRUE)
  )
}

test_that("crossLagRSummary returns a structured object for a CLPM fit", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  fit <- fit_clean_clpm()
  s <- crossLagRSummary(fit)
  expect_s3_class(s, "crossLagR_summary")
  expect_true(all(c("parameters", "other", "fit", "n_obs", "converged")
                  %in% names(s)))
  expect_s3_class(s$parameters, "data.frame")
  expect_true(all(c("label", "description", "est", "se", "z", "p",
                    "ci.lower", "ci.upper", "std.all")
                  %in% names(s$parameters)))
})

test_that("crossLagRSummary picks up unified labels (ar_x, cl_xy, etc.)", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  fit <- fit_clean_clpm()
  s <- crossLagRSummary(fit)
  expect_true("ar_x"  %in% s$parameters$label)
  expect_true("ar_y"  %in% s$parameters$label)
  expect_true("cl_xy" %in% s$parameters$label)
  expect_true("cl_yx" %in% s$parameters$label)
})

test_that("crossLagRSummary descriptions match the unified-label catalog", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  fit <- fit_clean_clpm()
  s <- crossLagRSummary(fit)
  ar_x_row <- s$parameters[s$parameters$label == "ar_x", ]
  expect_match(ar_x_row$description, "Autoregressive")
  cl_xy_row <- s$parameters[s$parameters$label == "cl_xy", ]
  expect_match(cl_xy_row$description, "Cross-lagged")
})

test_that("crossLagRSummary picks up unified labels on RICLPM (ar_x, cl_xy)", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  fit <- fit_clean_riclpm()
  s <- crossLagRSummary(fit)
  expect_true("ar_x"  %in% s$parameters$label)
  expect_true("cl_xy" %in% s$parameters$label)
})

test_that("crossLagRSummary returns fit indices as named numeric", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  fit <- fit_clean_clpm()
  s <- crossLagRSummary(fit)
  expect_type(s$fit, "double")
  expect_true(all(c("chisq", "df", "cfi", "tli", "rmsea", "srmr", "aic", "bic")
                  %in% names(s$fit)))
})

test_that("crossLagRSummary reports n_obs / n_miss (CLAUDE.md preference)", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  fit <- fit_clean_clpm()
  s <- crossLagRSummary(fit)
  expect_equal(s$n_obs, 500L)
  expect_true(is.integer(s$n_miss) || is.na(s$n_miss))
})

test_that("crossLagRSummary rejects non-lavaan input", {
  expect_error(crossLagRSummary(list()), "lavaan")
  expect_error(crossLagRSummary(1), "lavaan")
})

test_that("crossLagRSummary errors on bad ci", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  fit <- fit_clean_clpm()
  expect_error(crossLagRSummary(fit, ci = 0))
  expect_error(crossLagRSummary(fit, ci = 1))
  expect_error(crossLagRSummary(fit, ci = 1.5))
})

test_that("print.crossLagR_summary runs and shows expected sections", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  fit <- fit_clean_clpm()
  s <- crossLagRSummary(fit)
  out <- capture.output(print(s))
  expect_true(any(grepl("crossLagR model summary", out)))
  expect_true(any(grepl("Parameters \\(unified labels\\)", out)))
  expect_true(any(grepl("Fit statistics", out)))
  expect_true(any(grepl("chi-square", out)))
})

test_that("unifiedLabels returns a tidy catalog and prints by group", {
  cat_df <- unifiedLabels(print = FALSE)
  expect_s3_class(cat_df, "data.frame")
  expect_named(cat_df, c("label", "description", "group"))
  for (lab in c("ar_x", "cl_xy", "I_x", "T_var_x", "S_x", "d_var_x")) {
    expect_true(lab %in% cat_df$label, info = paste("missing:", lab))
  }
  out <- capture.output(unifiedLabels())
  expect_true(any(grepl("dynamics", out)))
  expect_true(any(grepl("random_intercepts", out)))
})
