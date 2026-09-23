# Bias-landscape helpers used by 08_liss.qmd (lavaan-fitted CLPM/RI-CLPM
# bias surfaces under contemporaneous effects). Called directly from the book,
# so exported.

utils::globalVariables(c(
  "bw_ratio", "ar", "rho_i", "cl_bias", "model", "regime",
  "cl_f", "bw_ratio_f", "rho_i_f"
))

#' Bias landscape with contemporaneous effects, by lavaan SEM fitting
#'
#' Fits the analyst's model in lavaan, as Bakker et al. (2021) do, for a grid
#' of true-parameter combinations. For each cell:
#'
#' \enumerate{
#'   \item Generate data with \code{simRICLPM_t()} (structural DGP with
#'     non-recursive contemporaneous effects). Innovation variance is
#'     normalised so the within stationary diagonal equals 1, keeping the
#'     \code{bw_ratio} interpretation consistent across cells.
#'   \item Fit the analyst's CLPM via \code{estimateCLPM()} + \code{lavaan}.
#'   \item Fit the analyst's RI-CLPM via \code{estimateRICLPM()} + \code{lavaan}.
#'   \item Extract the \code{cl_xy} estimate from each fit; bias is that
#'     estimate minus the structural lagged cross-lag \code{cl}.
#' }
#'
#' Single (large-sample) rep per cell as a population proxy. With
#' \code{sample_size = 1000} and \code{waves = 5} the OLS-equivalent slope
#' precision is well below the figure's visual resolution; converting to
#' multi-rep Monte Carlo just adds visual jitter without changing the picture.
#'
#' @param ar_seq Numeric vector of true autoregressions for the x axis.
#' @param bw_seq Numeric vector of between-to-within variance ratios for
#'   facet columns.
#' @param rho_seq Numeric vector of trait correlations for facet rows.
#' @param cl_values Numeric vector of true within-person cross-lags;
#'   becomes an outer facet column so multiple panels appear side-by-side.
#' @param contemp_xy,contemp_yx Contemporaneous coefficients.
#' @param waves Number of waves per simulated panel.
#' @param sample_size Persons per cell.
#' @param seed Optional integer seed.
#' @param verbose If \code{TRUE}, prints a progress dot per cell.
#' @return Data frame with columns \code{ar}, \code{bw_ratio}, \code{rho_i},
#'   \code{cl}, \code{model} (\code{"CLPM"} or \code{"RI-CLPM"}), and
#'   \code{cl_bias}.
#' @export
clpm_bias_landscape_lavaan <- function(
  ar_seq = seq(0, 0.6, by = 0.1),
  bw_seq = c(0, 0.3, 0.6, 0.9),
  rho_seq = c(0, 0.3, 0.6, 0.9),
  cl_values = c(0, 0.2),
  contemp_xy = 0.10,
  contemp_yx = 0.30,
  waves = 5,
  sample_size = 1000,
  seed = 1,
  verbose = FALSE
) {
  if (!is.null(seed)) set.seed(seed)
  grid <- expand.grid(
    ar       = ar_seq,
    bw_ratio = bw_seq,
    rho_i    = rho_seq,
    cl       = cl_values
  )
  clpm_syn <- estimateCLPM(waves = waves)
  riclpm_syn <- estimateRICLPM(waves = waves)

  fit_cl_xy <- function(syn, dat) {
    fit <- tryCatch(
      suppressWarnings(lavaan::lavaan(syn,
        data = dat,
        meanstructure = FALSE
      )),
      error = function(e) NULL
    )
    if (is.null(fit) || !lavaan::lavInspect(fit, "converged")) {
      return(NA_real_)
    }
    pt <- lavaan::parameterEstimates(fit)
    v <- pt$est[pt$label == "cl_xy"]
    if (length(v)) v[1] else NA_real_
  }

  rows <- vector("list", 2 * nrow(grid))
  k <- 1L
  for (i in seq_len(nrow(grid))) {
    r <- grid[i, ]

    # Pre-compute stationary within Psi[1,1] from structural matrices so we
    # can scale Sigma_eps such that the stationary within DIAGONAL = 1, which
    # makes `bw_ratio` interpretable as the trait-to-within ratio directly.
    Gamma <- matrix(c(0, contemp_yx, contemp_xy, 0), 2, byrow = TRUE)
    Bmat <- matrix(c(r$ar, r$cl, r$cl, r$ar), 2, byrow = TRUE)
    Minv <- solve(diag(2) - Gamma)
    Mhat <- Minv %*% Bmat
    # Stationarity guard: reduced-form M must have all |eigenvalues| < 1.
    # Cells outside the stationary region yield NA bias.
    eig_max <- max(abs(eigen(Mhat, only.values = TRUE)$values))
    if (eig_max >= 0.999) {
      rows[[k]] <- data.frame(
        ar = r$ar, bw_ratio = r$bw_ratio, rho_i = r$rho_i, cl = r$cl,
        model = "CLPM", cl_bias = NA_real_, stringsAsFactors = FALSE
      )
      k <- k + 1L
      rows[[k]] <- data.frame(
        ar = r$ar, bw_ratio = r$bw_ratio, rho_i = r$rho_i, cl = r$cl,
        model = "RI-CLPM", cl_bias = NA_real_, stringsAsFactors = FALSE
      )
      k <- k + 1L
      next
    }
    Sred <- Minv %*% diag(2) %*% t(Minv)
    Mkron <- Mhat %x% Mhat
    Psi_norm_to_one <- matrix(
      solve(diag(4) - Mkron, as.vector(Sred)), 2, 2
    )
    psi11 <- Psi_norm_to_one[1, 1]
    inn_var <- 1 / psi11

    bw_floor <- max(r$bw_ratio, 1e-8)
    sim <- simRICLPM_t(
      waves = waves,
      beta_x = r$ar, beta_y = r$ar,
      omega_xy = r$cl, omega_yx = r$cl,
      contemp_xy = contemp_xy, contemp_yx = contemp_yx,
      var_x = inn_var, var_y = inn_var, cov_xy = 0,
      var_BX = bw_floor, var_BY = bw_floor,
      cov_BXBY = r$rho_i * bw_floor,
      sample_size = sample_size,
      seed = NULL
    )

    cl_clpm <- fit_cl_xy(clpm_syn, sim$data)
    cl_riclpm <- fit_cl_xy(riclpm_syn, sim$data)

    rows[[k]] <- data.frame(
      ar = r$ar, bw_ratio = r$bw_ratio, rho_i = r$rho_i, cl = r$cl,
      model = "CLPM", cl_bias = cl_clpm - r$cl,
      stringsAsFactors = FALSE
    )
    k <- k + 1L
    rows[[k]] <- data.frame(
      ar = r$ar, bw_ratio = r$bw_ratio, rho_i = r$rho_i, cl = r$cl,
      model = "RI-CLPM", cl_bias = cl_riclpm - r$cl,
      stringsAsFactors = FALSE
    )
    k <- k + 1L

    if (isTRUE(verbose) && i %% 10 == 0) {
      message(sprintf("  cell %d / %d", i, nrow(grid)))
    }
  }
  out <- do.call(rbind, rows)
  out$model <- factor(out$model, levels = c("CLPM", "RI-CLPM"))
  out
}

# ---- "Naive" sequential simulator for contemporaneous effects --------------
#
# What happens if you try to extend standard simCLPM/simRICLPM with
# contemporaneous effects by just adding a `gamma * x_t` term to the y
# equation INSIDE the existing sequential loop? You get only ONE direction:
# the model can include x_t -> y_t (because x_t is already drawn when y_t is
# generated), but NOT y_t -> x_t (because y_t doesn't exist yet when x_t is
# drawn). The bidirectional non-recursive system simply cannot be generated
# sequentially -- there is no valid ordering.
#
# This function exists to make that failure mode concrete and quantitatively
# comparable to the proper reduced-form simulator (simRICLPM_t). It is
# deliberately deficient.

.sim_naive_contemp <- function(N, waves, ar, cl, contemp_xy,
                               var_x = 1, var_y = 1, cov_xy = 0,
                               var_BX = 0, var_BY = 0, cov_BXBY = 0) {
  has_trait <- var_BX > 1e-10
  if (has_trait) {
    LB <- chol(matrix(c(var_BX, cov_BXBY, cov_BXBY, var_BY), 2))
    trait <- matrix(stats::rnorm(2 * N), N, 2) %*% LB
  } else {
    trait <- matrix(0, N, 2)
  }

  Sigma_eps <- matrix(c(var_x, cov_xy, cov_xy, var_y), 2)
  L_eps <- chol(Sigma_eps)

  xi <- array(NA_real_, c(N, waves, 2))
  xi[, 1, ] <- matrix(stats::rnorm(2 * N), N, 2) %*% L_eps

  for (t in 2:waves) {
    eps <- matrix(stats::rnorm(2 * N), N, 2) %*% L_eps
    # Step 1: x_t depends only on lagged variables (no contemp y_t -> x_t)
    x_curr <- ar * xi[, t - 1, 1] + cl * xi[, t - 1, 2] + eps[, 1]
    # Step 2: y_t depends on lagged AND on the just-drawn x_t (contemp x_t -> y_t)
    y_curr <- ar * xi[, t - 1, 2] + cl * xi[, t - 1, 1] +
      contemp_xy * x_curr + eps[, 2]
    xi[, t, 1] <- x_curr
    xi[, t, 2] <- y_curr
  }

  z <- xi
  for (t in seq_len(waves)) z[, t, ] <- z[, t, ] + trait

  out <- data.frame(matrix(NA_real_, nrow = N, ncol = 2 * waves))
  col_names <- character(2 * waves)
  for (t in seq_len(waves)) {
    out[, 2 * t - 1] <- z[, t, 1]
    out[, 2 * t] <- z[, t, 2]
    col_names[2 * t - 1] <- paste0("x", t)
    col_names[2 * t] <- paste0("y", t)
  }
  colnames(out) <- col_names
  out
}

#' Bias landscape under a naive sequential contemp simulator
#'
#' Same interface as \code{clpm_bias_landscape_lavaan()}, but uses a
#' deliberately deficient internal simulator (\code{.sim_naive_contemp()})
#' that can only carry one direction of contemporaneous effect because
#' standard sequential simulation has no ordering for the bidirectional
#' loop. The resulting landscape illustrates what an analyst who tried to
#' bolt contemporaneous structure onto \code{simCLPM()} / \code{simRICLPM()}
#' without solving the simultaneous system would actually generate.
#'
#' @inheritParams clpm_bias_landscape_lavaan
#' @return Data frame with columns \code{ar}, \code{bw_ratio}, \code{rho_i},
#'   \code{cl}, \code{model} (\code{"CLPM"} or \code{"RI-CLPM"}), and
#'   \code{cl_bias}.
#' @export
clpm_bias_landscape_naive_lavaan <- function(
  ar_seq = seq(0, 0.5, by = 0.1),
  bw_seq = c(0, 0.3, 0.6, 0.9),
  rho_seq = c(0, 0.3, 0.6, 0.9),
  cl_values = c(0, 0.2),
  contemp_xy = 0.20,
  contemp_yx = 0.20, # accepted but SILENTLY DROPPED -- see details
  waves = 6,
  sample_size = 3000,
  seed = 1
) {
  if (!is.null(seed)) set.seed(seed)
  grid <- expand.grid(
    ar       = ar_seq,
    bw_ratio = bw_seq,
    rho_i    = rho_seq,
    cl       = cl_values
  )
  clpm_syn <- estimateCLPM(waves = waves)
  riclpm_syn <- estimateRICLPM(waves = waves)

  fit_cl_xy <- function(syn, dat) {
    fit <- tryCatch(
      suppressWarnings(lavaan::lavaan(syn,
        data = dat,
        meanstructure = FALSE
      )),
      error = function(e) NULL
    )
    if (is.null(fit) || !lavaan::lavInspect(fit, "converged")) {
      return(NA_real_)
    }
    pt <- lavaan::parameterEstimates(fit)
    v <- pt$est[pt$label == "cl_xy"]
    if (length(v)) v[1] else NA_real_
  }

  rows <- vector("list", 2 * nrow(grid))
  k <- 1L
  for (i in seq_len(nrow(grid))) {
    r <- grid[i, ]
    dat <- .sim_naive_contemp(
      N = sample_size,
      waves = waves,
      ar = r$ar,
      cl = r$cl,
      contemp_xy = contemp_xy, # only direction captured
      var_x = 1, var_y = 1, cov_xy = 0,
      var_BX = r$bw_ratio,
      var_BY = r$bw_ratio,
      cov_BXBY = r$rho_i * r$bw_ratio
    )
    cl_clpm <- fit_cl_xy(clpm_syn, dat)
    cl_riclpm <- fit_cl_xy(riclpm_syn, dat)
    rows[[k]] <- data.frame(
      ar = r$ar, bw_ratio = r$bw_ratio, rho_i = r$rho_i, cl = r$cl,
      model = "CLPM", cl_bias = cl_clpm - r$cl,
      stringsAsFactors = FALSE
    )
    k <- k + 1L
    rows[[k]] <- data.frame(
      ar = r$ar, bw_ratio = r$bw_ratio, rho_i = r$rho_i, cl = r$cl,
      model = "RI-CLPM", cl_bias = cl_riclpm - r$cl,
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
  out <- do.call(rbind, rows)
  out$model <- factor(out$model, levels = c("CLPM", "RI-CLPM"))
  out
}

#' Plot combined Bakker-style bias landscape (left + right halves)
#'
#' Stacks two bias-landscape data frames (e.g. one regime with no
#' contemporaneous effects, one with) into a single faceted figure with an
#' outer column distinguishing the two regimes.
#'
#' @param df_left,df_right Data frames with columns \code{ar}, \code{bw_ratio},
#'   \code{rho_i}, \code{cl}, \code{model}, \code{cl_bias} -- e.g. from
#'   \code{clpm_bias_landscape_lavaan()} or \code{clpm_bias_landscape_naive_lavaan()}.
#' @param left_label,right_label Strings labelling the regime strip.
#' @return A \code{ggplot} object.
#' @export
plot_clpm_bias_landscape_combined <- function(
  df_left, df_right,
  left_label = "beta[Sim] == 0",
  right_label = "beta[Sim] == 0.20"
) {
  df_left$regime <- factor(left_label, levels = c(left_label, right_label))
  df_right$regime <- factor(right_label, levels = c(left_label, right_label))
  df <- rbind(df_left, df_right)

  df$rho_i_f <- factor(df$rho_i,
    levels = sort(unique(df$rho_i)),
    labels = paste0("rho[I] == ", sort(unique(df$rho_i)))
  )
  df$bw_ratio_f <- factor(df$bw_ratio,
    levels = sort(unique(df$bw_ratio)),
    labels = paste0("r == ", sort(unique(df$bw_ratio)))
  )
  df$cl_f <- factor(df$cl,
    levels = sort(unique(df$cl)),
    labels = paste0("beta[CL] == ", sort(unique(df$cl)))
  )

  ggplot2::ggplot(df, ggplot2::aes(ar, cl_bias,
    color = model, shape = model,
    linetype = model
  )) +
    ggplot2::geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
    ggplot2::geom_line(linewidth = 0.45) +
    ggplot2::geom_point(size = 1.0) +
    ggplot2::facet_grid(
      cl_f + rho_i_f ~ regime + bw_ratio_f,
      labeller = ggplot2::label_parsed
    ) +
    ggplot2::scale_color_manual(values = c(
      "CLPM" = "#2C3E50",
      "RI-CLPM" = "#C0392B"
    )) +
    ggplot2::scale_shape_manual(values = c("CLPM" = 15, "RI-CLPM" = 17)) +
    ggplot2::scale_linetype_manual(values = c(
      "CLPM" = "solid",
      "RI-CLPM" = "solid"
    )) +
    ggplot2::labs(
      x = expression(beta[Autoregressive]),
      y = expression(E ~ "[" * hat(beta)[CL] - beta[CL] * "]"),
      color = "Analysis Model", shape = "Analysis Model",
      linetype = "Analysis Model"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      strip.text.x = ggplot2::element_text(size = 9, face = "bold"),
      strip.text.y = ggplot2::element_text(size = 9, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey92", color = NA),
      panel.spacing = ggplot2::unit(0.2, "lines"),
      legend.position = "right"
    )
}
