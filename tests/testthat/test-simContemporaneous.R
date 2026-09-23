# Tests for simCLPM_t and simRICLPM_t: simulators that add a non-recursive
# contemporaneous causal block to the standard lagged CLPM / RI-CLPM DGPs.

test_that("simCLPM_t returns expected list shape and column names", {
  sim <- simCLPM_t(waves = 4, sample_size = 200, seed = 1)
  expect_named(sim, c("model", "data", "parameters", "reduced_form"))
  expect_equal(ncol(sim$data), 8L)
  expect_equal(nrow(sim$data), 200L)
  expect_equal(colnames(sim$data),
               c("x1", "y1", "x2", "y2", "x3", "y3", "x4", "y4"))
})

test_that("simCLPM_t reduces to lagged CLPM when contemp_xy = contemp_yx = 0", {
  sim <- simCLPM_t(waves = 5,
                   beta_x = 0.3, beta_y = 0.3,
                   omega_xy = 0.05, omega_yx = 0.05,
                   contemp_xy = 0, contemp_yx = 0,
                   var_x = 1, var_y = 1, cov_xy = 0,
                   sample_size = 200, seed = 2)
  expect_equal(sim$reduced_form$M,
               matrix(c(0.3, 0.05, 0.05, 0.3), 2, 2, byrow = TRUE),
               tolerance = 1e-12)
  expect_equal(sim$reduced_form$Sigma_reduced, diag(2), tolerance = 1e-12)
})

test_that("simCLPM_t reduced-form M matches (I - Gamma)^-1 B by construction", {
  sim <- simCLPM_t(waves = 4, beta_x = 0.3, beta_y = 0.4,
                   omega_xy = 0.05, omega_yx = 0.10,
                   contemp_xy = 0.10, contemp_yx = 0.15,
                   sample_size = 100, seed = 3)
  Gamma <- matrix(c(0, 0.15, 0.10, 0), 2, byrow = TRUE)
  B     <- matrix(c(0.3, 0.10, 0.05, 0.4), 2, byrow = TRUE)
  expect_equal(sim$reduced_form$M,
               solve(diag(2) - Gamma) %*% B, tolerance = 1e-12)
})

test_that("simCLPM_t errors when (I - Gamma) is singular", {
  expect_error(
    simCLPM_t(waves = 3, contemp_xy = 1.0, contemp_yx = 1.0,
              sample_size = 100, seed = 4),
    "singular"
  )
})

test_that("simCLPM_t warns on non-stationary reduced form", {
  expect_warning(
    simCLPM_t(waves = 3, beta_x = 0.95, beta_y = 0.95,
              omega_xy = 0.2, omega_yx = 0.2,
              contemp_xy = 0.5, contemp_yx = 0.5,
              sample_size = 50, seed = 5),
    "non-stationary"
  )
})

test_that("simCLPM_t empirical lag-1 cross-products approximate M", {
  skip_on_cran()
  sim <- simCLPM_t(waves = 10,
                   beta_x = 0.30, beta_y = 0.30,
                   omega_xy = 0.05, omega_yx = 0.05,
                   contemp_xy = 0.10, contemp_yx = 0.15,
                   var_x = 1, var_y = 1, cov_xy = 0,
                   sample_size = 5000, seed = 6)
  # Stack consecutive pairs across waves 2..10
  d <- sim$data
  N <- nrow(d)
  idx <- 2:10
  y_cur <- cbind(
    as.vector(as.matrix(d[, paste0("x", idx)])),
    as.vector(as.matrix(d[, paste0("y", idx)]))
  )
  y_lag <- cbind(
    as.vector(as.matrix(d[, paste0("x", idx - 1)])),
    as.vector(as.matrix(d[, paste0("y", idx - 1)]))
  )
  # OLS lag-1 regression: M_hat = (y_lag' y_lag)^-1 y_lag' y_cur, transposed
  M_hat <- solve(crossprod(y_lag), crossprod(y_lag, y_cur))
  expect_equal(t(M_hat), sim$reduced_form$M, tolerance = 0.05)
})

# ---- simRICLPM_t -----------------------------------------------------------

test_that("simRICLPM_t returns expected shape", {
  sim <- simRICLPM_t(waves = 4, sample_size = 200, seed = 7)
  expect_named(sim, c("model", "data", "parameters", "reduced_form"))
  expect_equal(ncol(sim$data), 8L)
  expect_equal(colnames(sim$data),
               c("x1", "y1", "x2", "y2", "x3", "y3", "x4", "y4"))
})

test_that("simRICLPM_t reduces to RI-CLPM when contemp params are 0", {
  sim <- simRICLPM_t(waves = 5,
                     beta_x = 0.3, beta_y = 0.3,
                     omega_xy = 0.05, omega_yx = 0.05,
                     contemp_xy = 0, contemp_yx = 0,
                     var_x = 1, var_y = 1, cov_xy = 0,
                     var_BX = 1, var_BY = 1, cov_BXBY = 0.5,
                     sample_size = 200, seed = 8)
  expect_equal(sim$reduced_form$M,
               matrix(c(0.3, 0.05, 0.05, 0.3), 2, 2, byrow = TRUE),
               tolerance = 1e-12)
  expect_equal(sim$reduced_form$Sigma_reduced, diag(2), tolerance = 1e-12)
})

test_that("simRICLPM_t recovers between-person variance under MC", {
  skip_on_cran()
  sim <- simRICLPM_t(waves = 6, beta_x = 0.3, beta_y = 0.3,
                     omega_xy = 0, omega_yx = 0,
                     contemp_xy = 0, contemp_yx = 0,
                     var_x = 1, var_y = 1, cov_xy = 0,
                     var_BX = 2, var_BY = 2, cov_BXBY = 1,
                     sample_size = 3000, seed = 9)
  # Person-level means across waves should reflect trait variance
  d <- sim$data
  px <- rowMeans(d[, paste0("x", 1:6)])
  py <- rowMeans(d[, paste0("y", 1:6)])
  # Var(person mean) = trait variance + within variance / T (T = 6 waves)
  # Within stationary variance = 1 / (1 - 0.3^2) = 1.0989, so the
  # within contribution to person-mean variance is ~0.18
  expect_gt(var(px), 1.5)   # mostly trait
  expect_lt(var(px), 2.6)
  expect_gt(cov(px, py), 0.5)  # trait covariance survives
})

test_that("simRICLPM_t errors on non-PD inputs", {
  expect_error(
    simRICLPM_t(var_BX = 1, var_BY = 1, cov_BXBY = 2,
                sample_size = 100, seed = 10),
    "positive definite"
  )
})
