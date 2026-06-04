make_long_panel <- function(n_units = 200, n_waves = 4, between_sd = 1.0,
                            within_sd = 0.5, seed = 1) {
  set.seed(seed)
  alpha_x <- stats::rnorm(n_units, 0, between_sd)
  alpha_y <- stats::rnorm(n_units, 0, between_sd)
  do.call(rbind, lapply(seq_len(n_units), function(i) {
    data.frame(
      id   = i,
      wave = seq_len(n_waves),
      x    = alpha_x[i] + stats::rnorm(n_waves, 0, within_sd),
      y    = alpha_y[i] + stats::rnorm(n_waves, 0, within_sd)
    )
  }))
}

test_that("withinBetween dplyr method returns expected structure", {
  d <- make_long_panel()
  out <- withinBetween(d, id = "id", time = "wave",
                       vars = c("x", "y"), method = "dplyr")
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 2)
  expect_named(out, c("variable", "between_var", "within_var",
                      "total_var", "icc", "method"))
  expect_true(all(out$method == "dplyr"))
  expect_true(all(out$between_var > 0))
  expect_true(all(out$within_var > 0))
  expect_true(all(out$icc >= 0 & out$icc <= 1))
})

test_that("withinBetween recovers approximate ICC for high-ICC DGP", {
  d <- make_long_panel(between_sd = 1.5, within_sd = 0.3, seed = 22)
  out <- withinBetween(d, id = "id", time = "wave",
                       vars = c("x", "y"), method = "dplyr")
  expect_true(all(out$icc > 0.7))
})

test_that("withinBetween tso method requires exactly two variables", {
  d <- make_long_panel()
  expect_error(
    withinBetween(d, id = "id", time = "wave", vars = "x", method = "tso"),
    "exactly two variables"
  )
})

test_that("withinBetween tso method runs and returns lavaan-based estimates", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  d <- make_long_panel(n_units = 300, n_waves = 4,
                       between_sd = 1.2, within_sd = 0.5, seed = 33)
  out <- suppressWarnings(
    withinBetween(d, id = "id", time = "wave",
                  vars = c("x", "y"), method = "tso")
  )
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 2)
  expect_true(all(out$method == "tso"))
  expect_true(all(out$between_var > 0))
})

test_that("simWithinBetween + withinBetween round-trip yields sane bootstrap", {
  d <- make_long_panel(seed = 99)
  obs <- withinBetween(d, id = "id", time = "wave",
                       vars = c("x", "y"), method = "dplyr")
  sims <- simWithinBetween(obs, n_units = 100, n_waves = 4,
                           n_sims = 25, seed = 4)
  expect_s3_class(sims, "data.frame")
  expect_equal(nrow(sims), 25 * 2)
  expect_named(sims, c("variable", "sim", "between_var", "within_var"))
  expect_true(all(sims$between_var >= 0))
  expect_true(all(sims$within_var  >= 0))
})

test_that("simWithinBetween requires n_units and n_waves", {
  obs <- data.frame(variable = "x", between_var = 1, within_var = 0.5,
                    total_var = 1.5, icc = 2/3, method = "dplyr",
                    stringsAsFactors = FALSE)
  expect_error(simWithinBetween(obs), "n_units")
})
