test_that("clpm_bias_matrix returns no bias when there is no trait variance", {
  psi <- diag(2)
  b   <- diag(0.3, 2)
  bias <- clpm_bias_matrix(matrix(0, 2, 2), psi, b)
  expect_equal(bias, matrix(0, 2, 2))
})

test_that("clpm_bias_matrix equals (I - B) * Omega * (Omega + Psi)^-1 directly", {
  omega <- matrix(c(0.8, 0.4, 0.4, 0.8), 2)
  psi   <- matrix(c(1.0, 0.2, 0.2, 1.0), 2)
  b     <- matrix(c(0.4, 0.1, 0.0, 0.5), 2, byrow = TRUE)
  expect_equal(
    clpm_bias_matrix(omega, psi, b),
    (diag(2) - b) %*% omega %*% solve(omega + psi)
  )
})

test_that("clpm_artifact_grid matches the general matrix form for the symmetric design", {
  rho_i <- 0.5
  ar    <- 0.3
  for (r in c(0.5, 1, 2, 4, 9)) {
    sig <- r # within var = 1, so trait var = r
    omega <- matrix(c(sig, rho_i * sig, rho_i * sig, sig), 2)
    psi   <- diag(2)
    b     <- diag(ar, 2)
    m     <- clpm_bias_matrix(omega, psi, b)
    g     <- clpm_artifact_grid(r, ar = ar, rho_i = rho_i)
    expect_equal(g$ar_bias, m[1, 1], tolerance = 1e-10)
    expect_equal(g$cl_bias, m[1, 2], tolerance = 1e-10)
  }
})

test_that("clpm_artifact_grid reproduces the chapter's closed-form numbers", {
  g_peak <- clpm_artifact_grid(cl_bias_peak_ratio(0.5), ar = 0.3, rho_i = 0.5)
  expect_equal(round(g_peak$cl_bias, 3), 0.094)

  g_high <- clpm_artifact_grid(9, ar = 0.3, rho_i = 0.5)
  expect_equal(round(g_high$ar_bias, 2), 0.61)
})

test_that("cl_bias is zero when traits are uncorrelated and grows with rho_i", {
  expect_equal(clpm_artifact_grid(2, ar = 0.3, rho_i = 0)$cl_bias, 0)
  b_low  <- clpm_artifact_grid(2, ar = 0.3, rho_i = 0.3)$cl_bias
  b_high <- clpm_artifact_grid(2, ar = 0.3, rho_i = 0.7)$cl_bias
  expect_gt(b_high, b_low)
})

test_that("cl_bias is maximised at r_star = 1/sqrt(1 - rho^2)", {
  rho_i  <- 0.5
  r_star <- cl_bias_peak_ratio(rho_i)
  rs     <- seq(0.05, 12, by = 0.01)
  cl     <- clpm_artifact_grid(rs, ar = 0.3, rho_i = rho_i)$cl_bias
  expect_equal(rs[which.max(cl)], r_star, tolerance = 0.02)
})

test_that("ar_bias is monotone increasing and saturates at (1 - ar)", {
  g <- clpm_artifact_grid(seq(0, 5000, length.out = 50), ar = 0.3, rho_i = 0.5)
  expect_true(all(diff(g$ar_bias) >= 0))
  expect_equal(tail(g$ar_bias, 1), 0.7, tolerance = 1e-2)
})

test_that("plot helpers return ggplot objects", {
  expect_s3_class(plot_cl_bias_peak(ar = 0.3, rho_i = 0.5), "ggplot")
  expect_s3_class(
    plot_cl_bias_surface(rho_i = 0.5, n_r = 20, n_ar = 20, contours = FALSE),
    "ggplot"
  )
})

test_that("clpmBias returns the expected shape and columns", {
  out <- clpmBias(bw_ratio = c(1, 2, 4), ar = 0.3, rho_i = 0.5)
  expect_named(out, c("bw_ratio", "ar", "rho_i", "icc",
                      "ar_bias", "ar_bias_pct",
                      "cl_bias", "cl_bias_pct"))
  expect_equal(nrow(out), 3L)
  expect_equal(out$icc, c(1, 2, 4) / (1 + c(1, 2, 4)))
})

test_that("clpmBias percentage columns follow the convention", {
  out <- clpmBias(bw_ratio = c(1, 2), ar = 0.3, rho_i = 0.5)
  expect_equal(out$ar_bias_pct, 100 * out$ar_bias / out$ar, tolerance = 1e-12)
  # true CL is 0 here, so cl_bias_pct is NA (division by zero)
  expect_true(all(is.na(out$cl_bias_pct)))
})

test_that("clpmBias matches the internal closed form", {
  out <- clpmBias(bw_ratio = c(1, 2, 4), ar = 0.3, rho_i = 0.5)
  g   <- clpm_artifact_grid(c(1, 2, 4), ar = 0.3, rho_i = 0.5)
  expect_equal(out$cl_bias, g$cl_bias)
  expect_equal(out$ar_bias, g$ar_bias)
})

test_that("clpmBias expand evaluates the full crossing", {
  out <- clpmBias(bw_ratio = c(1, 2), ar = c(0, 0.5), rho_i = 0.5, expand = TRUE)
  expect_equal(nrow(out), 4L)
})

test_that("clpmBias validates its inputs", {
  expect_error(clpmBias(bw_ratio = -1))
  expect_error(clpmBias(rho_i = 1.5))
})

test_that("clpmBias gives zero cross-lag bias when traits are uncorrelated", {
  expect_equal(clpmBias(bw_ratio = 2, ar = 0.3, rho_i = 0)$cl_bias, 0)
})

# ---- As-if artifact (population CLPM under RI-CLPM truth) ------------------

test_that("artifact grid matches matrix-form bias for non-zero true cl", {
  ar <- 0.3; cl <- 0.05; rho_i <- 0.5; r <- 2
  psi   <- .psi_stationary_symmetric(ar, cl)
  s2w   <- psi$sigma2_w
  Psi_w <- diag(s2w, 2)
  B     <- matrix(c(ar, cl, cl, ar), 2, byrow = TRUE)
  # build stationary Psi explicitly: Psi = Psi_w + B Psi B' iterated, or direct
  Psi_full <- matrix(c(1, psi$rho_w, psi$rho_w, 1), 2)
  expect_equal(Psi_full - B %*% Psi_full %*% t(B), Psi_w, tolerance = 1e-10)
  Omega <- r * matrix(c(1, rho_i, rho_i, 1), 2)
  slope <- (Omega + B %*% Psi_full) %*% solve(Omega + Psi_full)
  g     <- clpm_artifact_grid(r, ar = ar, cl = cl, rho_i = rho_i)
  expect_equal(g$ar_pop, slope[1, 1], tolerance = 1e-10)
  expect_equal(g$cl_pop, slope[1, 2], tolerance = 1e-10)
})

test_that("CLPM cross-lag inflates on the rising side of the trait-share curve", {
  # cl_pop is non-monotone in r (rises, peaks, decays); check rising side only.
  g <- clpm_artifact_grid(c(0, 0.3, 0.6, 1.0), ar = 0.3, cl = 0.02, rho_i = 0.5)
  expect_true(all(diff(g$cl_pop) > 0))
  expect_equal(g$cl_pop[1], 0.02, tolerance = 1e-10)
  # cl_pop > cl_true for r > 0 (inflation)
  expect_gt(g$cl_pop[4], 0.02)
})

test_that("tipping point round-trips through the closed form", {
  observed <- 0.08
  r_star   <- clpm_artifact_tipping(target = observed, ar = 0.3, cl = 0,
                                    rho_i = 0.5)
  expect_true(is.finite(r_star))
  back     <- clpm_artifact_grid(r_star, ar = 0.3, cl = 0, rho_i = 0.5)$cl_pop
  expect_equal(back, observed, tolerance = 1e-6)
})

test_that("tipping point returns NA when target is unreachable", {
  # very large target unreachable: max CL inflation is bounded
  r_star <- clpm_artifact_tipping(target = 5, ar = 0.3, cl = 0, rho_i = 0.5,
                                  interval = c(0, 50))
  expect_true(is.na(r_star))
})

test_that("closed-form tipping point finds falling-side roots the old rising-crossing rule missed", {
  # cl_pop(0) = cl > target here, and cl_pop(r) falls monotonically through
  # target as r grows -- there is no "rising" crossing, only a falling one.
  ar <- 0.4187; cl <- 0.0365; rho_i <- 0.0468; rel <- 0.9925
  target <- 0.00459
  r_star <- clpm_artifact_tipping(target = target, ar = ar, cl = cl,
                                  rho_i = rho_i, reliability = rel,
                                  interval = c(0, 50))
  expect_true(is.finite(r_star))
  back <- clpm_artifact_grid(r_star, ar = ar, cl = cl, rho_i = rho_i,
                             reliability = rel)$cl_pop
  expect_equal(back, target, tolerance = 1e-6)
})

test_that("closed-form tipping point picks the smaller of two valid roots", {
  # cl = 0 case: cl_pop(r) rises to a peak then decays, so a reachable target
  # below the peak has two non-negative roots; the smaller is reported.
  observed <- 0.08
  r_star <- clpm_artifact_tipping(target = observed, ar = 0.3, cl = 0, rho_i = 0.5)
  r_star_high <- clpm_artifact_grid(
    seq(r_star + 0.01, 20, length.out = 2000),
    ar = 0.3, cl = 0, rho_i = 0.5
  )
  other_crossing <- r_star_high$bw_ratio[which.min(abs(r_star_high$cl_pop - observed))]
  expect_lt(r_star, other_crossing)
})

test_that("clpmCLBias returns expected columns and shape", {
  out <- clpmCLBias(bw_ratio = c(1, 2), ar = 0.3, cl = c(0, 0.05),
                      rho_i = 0.5, expand = TRUE)
  expect_named(out, c("bw_ratio", "ar", "cl", "rho_i", "reliability",
                      "icc", "rho_w", "ar_pop", "cl_pop",
                      "ar_bias", "ar_bias_pct",
                      "cl_bias", "cl_bias_pct"))
  expect_equal(nrow(out), 4L)
})

test_that("percentage bias columns follow the documented convention", {
  out <- clpmCLBias(bw_ratio = 2, ar = 0.3, cl = c(0, 0.05), rho_i = 0.5,
                      expand = TRUE)
  # ar_bias_pct = 100 * ar_bias / ar
  expect_equal(out$ar_bias_pct, 100 * out$ar_bias / out$ar, tolerance = 1e-12)
  # cl_bias_pct = 100 * cl_bias / cl; NA (division by zero) when cl == 0
  expect_true(is.na(out$cl_bias_pct[1]))
  expect_equal(out$cl_bias_pct[2], 100 * out$cl_bias[2] / out$cl[2],
               tolerance = 1e-12)
})

test_that("clpmCLBias verbose=TRUE emits a syntax message", {
  expect_message(
    clpmCLBias(bw_ratio = 0.5, ar = 0.3, cl = 0.05, rho_i = 0.5,
                 verbose = TRUE),
    "between-to-within ratio of 0.5"
  )
})

test_that("clpmCLBias errors on non-stationary B", {
  expect_error(clpmCLBias(bw_ratio = 1, ar = 0.7, cl = 0.5),
               "Non-stationary")
})

test_that("clpmArtifactTipping returns implied trait shares", {
  tp <- clpmArtifactTipping(observed = 0.08, ar = 0.3,
                            cl = c(0, 0.02, 0.05), rho_i = 0.5)
  expect_named(tp, c("cl", "bw_ratio", "icc"))
  # larger true cl needs less trait variance to reach the target
  expect_true(tp$bw_ratio[1] > tp$bw_ratio[3])
})

test_that("artifact MC simulator matches closed form for non-zero true cl", {
  skip_on_cran()
  # use the existing simulator extended to non-zero cl: drive with B = [[ar, cl], [cl, ar]]
  set.seed(11)
  ar_t <- 0.30; cl_t <- 0.05; rho_i <- 0.5; r <- 2; N <- 1000; T_waves <- 20
  psi  <- .psi_stationary_symmetric(ar_t, cl_t)
  sigW <- sqrt(psi$sigma2_w)
  L_t  <- chol(r * matrix(c(1, rho_i, rho_i, 1), 2))
  trait <- matrix(stats::rnorm(2 * N), N, 2) %*% L_t
  # initialise within at stationary covariance
  L_w  <- chol(matrix(c(1, psi$rho_w, psi$rho_w, 1), 2))
  xi   <- array(NA_real_, c(N, T_waves, 2))
  xi[, 1, ] <- matrix(stats::rnorm(2 * N), N, 2) %*% L_w
  B   <- matrix(c(ar_t, cl_t, cl_t, ar_t), 2, byrow = TRUE)
  for (t in 2:T_waves) {
    xi[, t, ] <- xi[, t - 1, ] %*% t(B) +
      sigW * matrix(stats::rnorm(2 * N), N, 2)
  }
  x <- xi
  for (t in seq_len(T_waves)) x[, t, ] <- x[, t, ] + trait

  idx   <- 2:T_waves
  y1    <- as.vector(x[, idx, 1]);    y2 <- as.vector(x[, idx, 2])
  x1lag <- as.vector(x[, idx - 1, 1]); x2lag <- as.vector(x[, idx - 1, 2])
  X     <- cbind(1, x1lag, x2lag)
  b1    <- solve(crossprod(X), crossprod(X, y1))
  b2    <- solve(crossprod(X), crossprod(X, y2))
  ar_hat <- (b1[2] + b2[3]) / 2
  cl_hat <- (b1[3] + b2[2]) / 2
  pop    <- clpm_artifact_grid(r, ar = ar_t, cl = cl_t, rho_i = rho_i)
  expect_equal(ar_hat, pop$ar_pop, tolerance = 0.02)
  expect_equal(cl_hat, pop$cl_pop, tolerance = 0.02)
})

# ---- Bias landscape with contemporaneous effects (lavaan-fitted) ----------

test_that("clpm_bias_landscape_lavaan returns expected columns and shape", {
  skip_on_cran()
  df <- clpm_bias_landscape_lavaan(
    ar_seq = c(0, 0.3), bw_seq = c(0, 1), rho_seq = 0.5,
    cl_values = 0, contemp_xy = 0, contemp_yx = 0,
    waves = 3, sample_size = 300, seed = 1
  )
  expect_named(df, c("ar", "bw_ratio", "rho_i", "cl", "model", "cl_bias"))
  expect_equal(levels(df$model), c("CLPM", "RI-CLPM"))
})

test_that("clpm_bias_landscape_naive_lavaan returns expected columns and shape", {
  skip_on_cran()
  df <- clpm_bias_landscape_naive_lavaan(
    ar_seq = 0.3, bw_seq = c(0, 1), rho_seq = 0.5,
    cl_values = 0, contemp_xy = 0.1, contemp_yx = 0.1,
    waves = 3, sample_size = 300, seed = 1
  )
  expect_named(df, c("ar", "bw_ratio", "rho_i", "cl", "model", "cl_bias"))
})

test_that("plot_clpm_bias_landscape_combined returns a ggplot", {
  skip_on_cran()
  left <- clpm_bias_landscape_lavaan(
    ar_seq = 0.3, bw_seq = c(0, 1), rho_seq = 0.5, cl_values = 0,
    contemp_xy = 0, contemp_yx = 0, waves = 3, sample_size = 300, seed = 1
  )
  right <- clpm_bias_landscape_naive_lavaan(
    ar_seq = 0.3, bw_seq = c(0, 1), rho_seq = 0.5, cl_values = 0,
    contemp_xy = 0.1, contemp_yx = 0.1, waves = 3, sample_size = 300, seed = 1
  )
  expect_s3_class(plot_clpm_bias_landscape_combined(left, right), "ggplot")
})
