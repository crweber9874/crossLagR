#' @title compCTSEM
#' @description Create a continuous time model specification using the ctsem package.
#'
#' Wraps \code{ctsem::ctModel} to define a continuous-time structural equation model.
#' Returns the model object for subsequent fitting.
#'
#' @param data A long-format data frame.
#' @param id Character. Name of the ID column. Default is \code{"id"}.
#' @param manifest_var_names Character vector of manifest variable names. Default is \code{c("x", "y")}.
#' @param latent_var_names Character vector of latent variable names. Default is \code{c("x", "y")}.
#' @param LAMBDA Factor loading matrix. Default is \code{diag(2)}.
#' @param type Character. Type of ctsem model. Default is \code{"stanct"}.
#' @param time_points Integer. Number of time points. Default is 10.
#' @param cores Integer. Number of CPU cores for parallel computation.
#' @param chains Integer. Number of MCMC chains.
#' @param ... Additional arguments passed to \code{ctsem::ctModel}.
#' @return A list; returns model and drift matrix
#' #'
#' @examples
#' \dontrun{
#' model <- compCTSEM(data = my_long_data, time_points = 5)
#' }
#'
#' @export
compCTSEM <- function(data = dat_long,
                      id = "id",
                      manifest_var_names = c("x", "y"),
                      latent_var_names = c("x", "y"),
                      LAMBDA = diag(2),
                      type = "stanct",
                      time_points = 10,
                      cores = cores,
                      chains = chains,
                      ...) {
  my_model <- ctsem::ctModel(
    Tpoints = time_points,
    id = id,
    manifestNames = manifest_var_names,
    latentNames = latent_var_names,
    LAMBDA = LAMBDA,
    type = type,
    ...
  )
  return(my_model)
}
