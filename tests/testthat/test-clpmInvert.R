test_that("clpmInvert recovers parameters used to generate the estimates", {
  truths <- list(c(0.30, 1.50, 0.5), c(0.45, 0.60, 0.7), c(0.20, 4.00, 0.3))
  for (tr in truths) {
    pop <- clpmCLBias(bw_ratio = tr[2], ar = tr[1], cl = 0, rho_i = tr[3])
    got <- clpmInvert(pop$ar_pop, pop$cl_pop, rho_i = tr[3])
    expect_equal(got$ar_implied, tr[1], tolerance = 1e-3)
    expect_equal(got$bw_ratio_implied, tr[2], tolerance = 1e-3)
  }
})

test_that("clpmInvert is vectorised over rho_i and returns the documented shape", {
  out <- clpmInvert(0.692, 0.092, rho_i = c(0.5, 0.7, 0.9))
  expect_s3_class(out, "data.frame")
  expect_identical(names(out), c("rho_i", "ar_implied", "bw_ratio_implied", "icc_implied"))
  expect_equal(nrow(out), 3L)
  expect_equal(out$icc, out$bw_ratio_implied / (1 + out$bw_ratio_implied))
})

test_that("clpmInvert returns NA when no admissible solution exists", {
  out <- clpmInvert(0.692, 0.092, rho_i = 0.3)
  expect_true(is.na(out$bw_ratio_implied))
})

test_that("clpmInvert validates its inputs", {
  expect_error(clpmInvert(0.6, 0, rho_i = 0.5), "non-zero")
  expect_error(clpmInvert(0.6, 0.09, rho_i = 1), "strictly within")
})

test_that("clpmInvert round-trips against the forward map", {
  out <- clpmInvert(0.692, 0.092, rho_i = 0.6)
  fwd <- clpmCLBias(bw_ratio = out$bw_ratio_implied, ar = out$ar_implied, cl = 0, rho_i = 0.6)
  expect_equal(fwd$ar_pop, 0.692, tolerance = 1e-3)
  expect_equal(fwd$cl_pop, 0.092, tolerance = 1e-3)
})

test_that("clpmInvert is exact, not grid-limited", {
  pop <- clpmCLBias(bw_ratio = 1.5, ar = 0.30, cl = 0, rho_i = 0.5)
  got <- clpmInvert(pop$ar_pop, pop$cl_pop, rho_i = 0.5)
  expect_equal(got$ar_implied, 0.30, tolerance = 1e-10)
  expect_equal(got$bw_ratio_implied, 1.50, tolerance = 1e-10)
})

test_that("clpmInvert rejects algebraic solutions outside the parameter space", {
  # At rho_i = 0.3 the closed form solves at b = -63.8, which is not a
  # stationary autoregression; the admissibility guard must return NA.
  expect_true(is.na(clpmInvert(0.692, 0.092, rho_i = 0.3)$ar_implied))
  # Relaxing the guard exposes the raw algebraic root.
  loose <- clpmInvert(0.692, 0.092, rho_i = 0.3, b_max = 100)
  expect_lt(loose$ar_implied, -1)
  expect_gt(loose$bw_ratio_implied, 0)
})

test_that("clpmInvert recovers asymmetric truths from an asymmetric report", {
  fwd <- function(bx, by, rx, ry, p) {
    Om <- matrix(c(rx, p * sqrt(rx * ry), p * sqrt(rx * ry), ry), 2)
    V  <- Om %*% solve(Om + diag(2))
    B  <- diag(c(bx, by))
    B + (diag(2) - B) %*% V
  }
  B <- fwd(0.30, 0.50, 1.5, 0.8, 0.5)
  out <- clpmInvert(ar_obs = c(B[1, 1], B[2, 2]),
                    cl_obs = c(B[1, 2], B[2, 1]), rho_i = 0.5)
  expect_equal(out$ar_implied_x, 0.30, tolerance = 1e-5)
  expect_equal(out$ar_implied_y, 0.50, tolerance = 1e-5)
  expect_equal(out$bw_ratio_implied_x, 1.5, tolerance = 1e-5)
  expect_equal(out$bw_ratio_implied_y, 0.8, tolerance = 1e-5)
})

test_that("clpmInvert dispatches on symmetry and returns the matching shape", {
  sym <- clpmInvert(0.692, 0.092, rho_i = 0.5)
  expect_identical(names(sym), c("rho_i", "ar_implied", "bw_ratio_implied", "icc_implied"))

  # length-2 input that happens to be symmetric takes the closed-form path
  sym2 <- clpmInvert(c(0.692, 0.692), c(0.092, 0.092), rho_i = 0.5)
  expect_identical(names(sym2), c("rho_i", "ar_implied", "bw_ratio_implied", "icc_implied"))
  expect_equal(sym2$bw_ratio_implied, sym$bw_ratio_implied, tolerance = 1e-10)

  asym <- clpmInvert(c(0.70, 0.62), c(0.09, 0.04), rho_i = 0.5)
  expect_true(all(c("ar_implied_x", "bw_ratio_implied_y", "icc_implied_x") %in% names(asym)))
})

test_that("asymmetric solutions agree with the symmetric closed form when equal", {
  B <- local({
    p <- 0.5; r <- 1.5; b <- 0.3
    Om <- matrix(c(r, p * r, p * r, r), 2)
    V  <- Om %*% solve(Om + diag(2))
    b * diag(2) + (1 - b) * V
  })
  # nudge one entry so the asymmetric solver is used, then compare
  asym <- crossLagR:::.clpm_invert_asym(c(B[1, 1], B[2, 2]),
                                        c(B[1, 2], B[2, 1]), 0.5, 1)
  expect_equal(asym$bw_ratio_implied_x, 1.5, tolerance = 1e-5)
  expect_equal(asym$ar_implied_y, 0.3, tolerance = 1e-5)
})
