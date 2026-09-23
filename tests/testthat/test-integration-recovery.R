## End-to-end parameter-recovery integration tests: simFromSyntax() -> lavaan
## fit -> check the AR/CL estimates land near the true values. Existing
## test-estimateXXX.R files only check that the generated syntax string
## parses; these tests instead confirm the whole simulate-then-fit pipeline
## actually recovers known truth for each unified estimator family. Single
## large-N fit per estimator (not a full Monte Carlo) -- enough to catch a
## broken sim/estimate pairing without being a slow suite.
##
## Recovery checks use an absolute error bound (expect_lt(abs(est - truth), ...))
## rather than expect_equal(tolerance = ), because expect_equal's tolerance is
## RELATIVE (base R all.equal semantics): for true values in the 0.1-0.3 range,
## a "tight-looking" tolerance = 0.05 actually demands accuracy far beyond what
## a single N = 2000 draw can deliver.

recovers <- function(est, truth, abs_tol) {
  expect_lt(abs(est - truth), abs_tol)
}

test_that("CLPM recovers known AR/CL truth", {
  skip_on_cran()
  dat <- simFromSyntax("CLPM", waves = 4, sample_size = 2000,
                       ar_x = 0.3, ar_y = 0.3, cl_xy = 0.15, cl_yx = 0.1,
                       seed = 1)
  fit <- suppressWarnings(lavaan::lavaan(estimateCLPM(waves = 4), data = dat,
                                        meanstructure = TRUE))
  expect_true(lavaan::lavInspect(fit, "converged"))
  pe <- lavaan::parameterEstimates(fit)
  est <- function(l) pe$est[pe$label == l][1]
  recovers(est("ar_x"),  0.3,  0.05)
  recovers(est("ar_y"),  0.3,  0.05)
  recovers(est("cl_xy"), 0.15, 0.05)
  recovers(est("cl_yx"), 0.1,  0.05)
})

test_that("RICLPM recovers known AR/CL truth with trait variance present", {
  skip_on_cran()
  dat <- simFromSyntax("RICLPM", waves = 4, sample_size = 2000,
                       ar_x = 0.3, ar_y = 0.3, cl_xy = 0.15, cl_yx = 0.1,
                       var_between_x = 1, var_between_y = 1, seed = 2)
  fit <- suppressWarnings(lavaan::lavaan(estimateRICLPM(waves = 4), data = dat,
                                        meanstructure = TRUE))
  expect_true(lavaan::lavInspect(fit, "converged"))
  pe <- lavaan::parameterEstimates(fit)
  est <- function(l) pe$est[pe$label == l][1]
  recovers(est("ar_x"),  0.3,  0.05)
  recovers(est("ar_y"),  0.3,  0.05)
  recovers(est("cl_xy"), 0.15, 0.05)
  recovers(est("cl_yx"), 0.1,  0.05)
})

test_that("ALT recovers known AR/CL truth", {
  skip_on_cran()
  dat <- simFromSyntax("ALT", waves = 4, sample_size = 2000,
                       ar_x = 0.3, ar_y = 0.3, cl_xy = 0.15, cl_yx = 0.1,
                       seed = 3)
  fit <- suppressWarnings(lavaan::lavaan(estimateALT(waves = 4), data = dat,
                                        meanstructure = TRUE))
  expect_true(lavaan::lavInspect(fit, "converged"))
  pe <- lavaan::parameterEstimates(fit)
  est <- function(l) pe$est[pe$label == l][1]
  recovers(est("ar_x"),  0.3,  0.08)
  recovers(est("ar_y"),  0.3,  0.08)
  recovers(est("cl_xy"), 0.15, 0.08)
  recovers(est("cl_yx"), 0.1,  0.08)
})

test_that("LCMSR recovers known AR/CL truth", {
  skip_on_cran()
  dat <- simFromSyntax("LCMSR", waves = 4, sample_size = 2000,
                       ar_x = 0.3, ar_y = 0.3, cl_xy = 0.15, cl_yx = 0.1,
                       seed = 4)
  fit <- suppressWarnings(lavaan::lavaan(estimateLCMSR(waves = 4), data = dat,
                                        meanstructure = TRUE))
  expect_true(lavaan::lavInspect(fit, "converged"))
  pe <- lavaan::parameterEstimates(fit)
  est <- function(l) pe$est[pe$label == l][1]
  recovers(est("ar_x"),  0.3,  0.08)
  recovers(est("ar_y"),  0.3,  0.08)
  recovers(est("cl_xy"), 0.15, 0.08)
  recovers(est("cl_yx"), 0.1,  0.08)
})

test_that("Bollen & Brand recovers known AR/CL truth when coefficients are constrained", {
  skip_on_cran()
  args <- list(constrain_coefficients = TRUE)
  dat <- simFromSyntax("BB", waves = 4, sample_size = 2000,
                       ar_x = 0.3, ar_y = 0.3, cl_xy = 0.15, cl_yx = 0.1,
                       estimator_args = args, seed = 5)
  fit <- suppressWarnings(lavaan::lavaan(
    estimateBollen_and_Brand(waves = 4, constrain_coefficients = TRUE),
    data = dat, meanstructure = TRUE
  ))
  expect_true(lavaan::lavInspect(fit, "converged"))
  pe <- lavaan::parameterEstimates(fit)
  est <- function(l) pe$est[pe$label == l][1]
  recovers(est("ar_x"),  0.3,  0.05)
  recovers(est("ar_y"),  0.3,  0.05)
  recovers(est("cl_xy"), 0.15, 0.05)
  recovers(est("cl_yx"), 0.1,  0.05)
})

test_that("LGM and TSO fit their simulated data (no ar_x/cl_xy -- these are dynamics-free models)", {
  skip_on_cran()
  dat_lgm <- simFromSyntax("LGM", waves = 4, sample_size = 1000,
                           estimator_args = list(variable_type = "bivariate"),
                           seed = 6)
  fit_lgm <- suppressWarnings(lavaan::lavaan(
    estimateLGM(waves = 4, variable_type = "bivariate"),
    data = dat_lgm, meanstructure = TRUE
  ))
  expect_true(lavaan::lavInspect(fit_lgm, "converged"))

  dat_tso <- simFromSyntax("TSO", waves = 4, sample_size = 1000, seed = 7)
  fit_tso <- suppressWarnings(lavaan::lavaan(estimateTSO(waves = 4),
                                            data = dat_tso, meanstructure = TRUE))
  expect_true(lavaan::lavInspect(fit_tso, "converged"))
})

test_that("LCHANGE fits its simulated data and converges", {
  skip_on_cran()
  ## LCHANGE's ar_x/cl_xy are proportional-change and coupling parameters on
  ## change scores, not level-CLPM AR/CL -- signs and magnitudes are not
  ## directly comparable to the level-model tests above, so this checks
  ## convergence only, not point recovery.
  dat <- simFromSyntax("LCHANGE", waves = 5, sample_size = 1000,
                       ar_x = -0.2, ar_y = -0.2, cl_xy = 0.1, cl_yx = 0.1,
                       estimator_args = list(variable_type = "bivariate"),
                       seed = 8)
  fit <- suppressWarnings(lavaan::lavaan(
    estimateLChange(waves = 5, variable_type = "bivariate"),
    data = dat, meanstructure = TRUE
  ))
  expect_true(lavaan::lavInspect(fit, "converged"))
})
