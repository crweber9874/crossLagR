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

estimateSTARTS= function(
    data = dat_long, # Must follow x, y, xlag, ylag; from datLong
    residual_correlation = TRUE,
    ...
   ) {

  fit1 = lm(y ~ xlag + ylag + as.factor(id), data = data)
  fit2 = lm(x ~ xlag + ylag + as.factor(id), data = data)

  return(list(yeqn = fit1, xeqn = fit2))
}


