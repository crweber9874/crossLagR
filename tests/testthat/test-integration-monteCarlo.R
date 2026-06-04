## End-to-end integration tests for the unified-label Monte Carlo wrappers.
## Kept small (few trials, modest N) so the suite remains CI-friendly.

grid_default <- function() {
  data.frame(
    stability_p = 0.3, stability_q = 0.3,
    cross_p = 0.1, cross_q = 0.1,
    variance_p = 0.5, variance_q = 0.5,
    variance_between_x = 0.5, variance_between_y = 0.5,
    cov_pq = 0
  )
}

test_that("monteCarloLavaan returns one row per trial across estimator/DGP combos", {
  skip_on_cran()
  skip_if_not_installed("lavaan")

  combos <- list(
    c("CLPM",   "CLPM"),
    c("RICLPM", "RICLPM"),
    c("ALT",    "CLPM"),
    c("CLPM",   "RICLPM")
  )
  for (combo in combos) {
    set.seed(31)
    out <- suppressWarnings(monteCarloLavaan(
      estimator   = combo[1],
      dgp         = combo[2],
      trials      = 2,
      waves       = 4,
      sample_size = 300,
      param_grid  = grid_default(),
      verbose     = FALSE
    ))
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 2,
                 info = paste(combo, collapse = "/"))
    expect_true("ar_x"      %in% names(out))
    expect_true("converged" %in% names(out))
    expect_true("n_obs"     %in% names(out))
    expect_true("true_ar_x" %in% names(out))
    expect_true(all(out$estimator == combo[1]))
    expect_true(all(out$dgp       == combo[2]))
  }
})

test_that("monteCarloLavaan honors legacy lowercase dgp codes", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  set.seed(13)
  out <- suppressWarnings(monteCarloLavaan(
    estimator = "CLPM", dgp = "clpm",
    trials = 1, waves = 4, sample_size = 200,
    param_grid = grid_default(), verbose = FALSE
  ))
  expect_equal(unique(out$dgp), "CLPM")
})

test_that("monteCarloLavaan with clpmu DGP populates true_ar_x and cl_yx", {
  skip_on_cran()
  skip_if_not_installed("lavaan")

  grid <- data.frame(
    stability_p = 0.4, stability_q = 0.4,
    cross_p = 0.15, cross_q = 0.15,
    variance_p = 0.5, variance_q = 0.5,
    variance_between_x = 0.5, variance_between_y = 0.5,
    cov_pq = 0,
    confounder_p = 0.3, confounder_q = 0.3,
    confounder_variance = 1, confounder_stability = 0.4
  )
  set.seed(101)
  out <- suppressWarnings(monteCarloLavaan(
    estimator = "CLPM", dgp = "clpmu",
    trials = 1, waves = 4, sample_size = 300,
    param_grid = grid, verbose = FALSE
  ))
  expect_equal(out$true_ar_x[1], 0.4)
  expect_equal(out$true_cl_yx[1], 0.15)
})

test_that("monteCarloLavaan respects multi-row param grids", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  g <- rbind(
    cbind(grid_default(), id = 1),
    cbind(grid_default(), id = 2)
  )
  g$cross_q <- c(0.1, 0.25)
  set.seed(91)
  out <- suppressWarnings(monteCarloLavaan(
    estimator = "CLPM", dgp = "CLPM",
    trials = 1, waves = 4, sample_size = 250,
    param_grid = g[, setdiff(names(g), "id")], verbose = FALSE
  ))
  expect_equal(nrow(out), 2)
  expect_equal(unique(out$param_combo), c(1L, 2L))
  expect_equal(sort(out$true_cl_xy), c(0.10, 0.25))
})

test_that("monteCarloLavaan reports n_obs and converged columns", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  set.seed(44)
  out <- suppressWarnings(monteCarloLavaan(
    estimator = "CLPM", dgp = "CLPM",
    trials = 2, waves = 4, sample_size = 250,
    param_grid = grid_default(), verbose = FALSE
  ))
  expect_true(all(out$n_obs[out$converged %in% TRUE] == 250))
  expect_true(is.logical(out$converged))
})

test_that("monteCarloLavaan estimator_args flow through to syntax builder", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  set.seed(77)
  out <- suppressWarnings(monteCarloLavaan(
    estimator      = "CLPM",
    dgp            = "CLPM",
    trials         = 1,
    waves          = 4,
    sample_size    = 250,
    param_grid     = grid_default(),
    estimator_args = list(constrain_beta = FALSE, constrain_omega = FALSE),
    verbose        = FALSE
  ))
  expect_s3_class(out, "data.frame")
  expect_true("ar_x" %in% names(out))
})
