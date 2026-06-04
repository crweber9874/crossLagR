## Unit tests for iccReport + integrations.

make_wide_panel <- function(n = 300, T = 4, between_sd = 1.0,
                            within_sd = 0.5, seed = 1) {
  set.seed(seed)
  alpha_x <- stats::rnorm(n, 0, between_sd)
  alpha_y <- stats::rnorm(n, 0, between_sd)
  x <- alpha_x + matrix(stats::rnorm(n * T, 0, within_sd), n, T)
  y <- alpha_y + matrix(stats::rnorm(n * T, 0, within_sd), n, T)
  d <- as.data.frame(cbind(x, y))
  names(d) <- c(paste0("x", 1:T), paste0("y", 1:T))
  d
}

test_that("iccReport wide format returns expected structure", {
  d <- make_wide_panel(seed = 1)
  out <- iccReport(d, vars = c("x", "y"), quiet = TRUE)
  expect_s3_class(out, "crossLagR_icc")
  expect_named(out, c("variable", "between_var", "within_var",
                      "total_var", "icc", "n_obs", "n_waves",
                      "exceeds_threshold"))
  expect_equal(nrow(out), 2)
  expect_true(all(out$icc > 0 & out$icc < 1))
  expect_equal(attr(out, "threshold"), 0.65)
})

test_that("iccReport recovers approximate ICC for known DGP", {
  ## True ICC = 1^2 / (1^2 + 0.5^2) = 0.8
  d <- make_wide_panel(n = 1000, T = 6,
                       between_sd = 1.0, within_sd = 0.5, seed = 7)
  out <- iccReport(d, vars = c("x", "y"), quiet = TRUE)
  expect_true(all(abs(out$icc - 0.8) < 0.08))
})

test_that("iccReport warns when ICC exceeds threshold and cites Shiny", {
  d <- make_wide_panel(n = 500, T = 5,
                       between_sd = 1.5, within_sd = 0.4, seed = 21)
  expect_warning(
    iccReport(d, vars = c("x", "y"), threshold = 0.65),
    "launch_shiny_sim"
  )
})

test_that("iccReport quiet = TRUE suppresses warning", {
  d <- make_wide_panel(n = 500, T = 5,
                       between_sd = 1.5, within_sd = 0.4, seed = 22)
  expect_no_warning(
    iccReport(d, vars = c("x", "y"), threshold = 0.65, quiet = TRUE)
  )
})

test_that("iccReport does not warn when ICC is below threshold", {
  d <- make_wide_panel(n = 500, T = 5,
                       between_sd = 0.3, within_sd = 1.5, seed = 23)
  expect_no_warning(
    iccReport(d, vars = c("x", "y"), threshold = 0.65)
  )
})

test_that("iccReport flags exceeds_threshold correctly", {
  d <- make_wide_panel(n = 500, T = 5,
                       between_sd = 1.5, within_sd = 0.4, seed = 24)
  out <- iccReport(d, vars = c("x", "y"), threshold = 0.65, quiet = TRUE)
  expect_true(all(out$exceeds_threshold))
})

test_that("iccReport long format routes through withinBetween", {
  d_wide <- make_wide_panel(n = 200, T = 4, seed = 31)
  d_wide$id <- seq_len(nrow(d_wide))
  long <- stats::reshape(
    d_wide,
    direction = "long",
    varying   = list(x = paste0("x", 1:4), y = paste0("y", 1:4)),
    v.names   = c("x", "y"),
    timevar   = "wave",
    idvar     = "id"
  )
  out <- iccReport(long, vars = c("x", "y"), format = "long",
                   id = "id", time = "wave", quiet = TRUE)
  expect_equal(nrow(out), 2)
  expect_true(all(out$icc > 0 & out$icc < 1))
})

test_that("iccReport rejects invalid inputs", {
  d <- make_wide_panel(seed = 1)
  expect_error(iccReport(list(), vars = "x"))
  expect_error(iccReport(d, vars = character(0)))
  expect_error(iccReport(d, vars = "x", threshold = NA))
  expect_error(iccReport(d, vars = "x", format = "long"))   # missing id/time
})

test_that("iccReport warns when stem matches nothing", {
  d <- make_wide_panel(seed = 1)
  expect_warning(
    iccReport(d, vars = "z", quiet = TRUE),
    "no columns matching"
  )
})

test_that("print.crossLagR_icc shows threshold and flagged variables", {
  d <- make_wide_panel(n = 500, T = 4,
                       between_sd = 1.5, within_sd = 0.4, seed = 25)
  out <- iccReport(d, vars = c("x", "y"), threshold = 0.65, quiet = TRUE)
  printed <- capture.output(print(out))
  expect_true(any(grepl("threshold", printed)))
  expect_true(any(grepl("launch_shiny_sim", printed)))
})

test_that("monteCarloLavaan output gains icc_x and icc_y columns", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  grid <- data.frame(
    stability_p = 0.3, stability_q = 0.3,
    cross_p = 0.1, cross_q = 0.1,
    variance_p = 0.5, variance_q = 0.5,
    variance_between_x = 0.5, variance_between_y = 0.5,
    cov_pq = 0
  )
  set.seed(909)
  out <- suppressWarnings(monteCarloLavaan(
    estimator = "CLPM", dgp = "CLPM",
    trials = 2, waves = 4, sample_size = 250,
    param_grid = grid, verbose = FALSE
  ))
  expect_true(all(c("icc_x", "icc_y", "icc_exceeds") %in% names(out)))
  expect_true(all(is.finite(out$icc_x)))
  expect_true(all(is.finite(out$icc_y)))
})
