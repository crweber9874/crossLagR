#' @title estimateHLM
#' @description Generates model for HLM model in brms
#'
#' @param data  Long data frame, following structure from datLong() names() -> x, y, xlag, ylag
#' @param chains Specify the number of sampling chains, default 3
#' @param residual_correlations Logical, should residual correlations be estimated?
#' @param iterations Number of iterations for the model, default 3000
#' @param cores Number of cores to use for the model, default 10
#'
#'
#' @return brms object
#'
#' @details
#' This function generates a model for a Hierarchical Linear Model (HLM) using the brms package.
#' The data must be in long form, and include the lagged variables xlag and ylag.
#' The model includes random intercepts and slopes for x and y.
#'
#'@examples
#'# Not Run

#' @export
#'

estimateHLM = function(
    data = dat_long, # Must follow x, y, xlag, ylag; from datLong
    residual_correlation = TRUE,
    chains = 3,
    cores = 10,
    iterations = 3000,
    variational = FALSE,
    control = list(adapt_delta = 0.95),
    ...
   ) {

  yeqn = brms::bf(y ~ xlag + ylag + (1 + ylag + xlag|p|id))
  xeqn = brms::bf(x ~ xlag + ylag + (1 + ylag + xlag|p|id))

    if (variational) {
    algorithm_val <- "meanfield"
    tol_rel_obj_val <- 1e-9
  }
    else {
    algorithm_val <- "sampling" # or "optimizing"
    tol_rel_obj_val <- NULL # or another relevant value
  }
  fit = brms::brm(yeqn + xeqn + brms::set_rescor(residual_correlation),
             data = data,
             family = gaussian(),
             cores = cores,
             chains = chains,
             iter = iterations,
             control = control,
             algorithm = algorithm_val,
             ...)
  return(fit)
}


