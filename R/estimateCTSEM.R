#' @title estimateCTSEM
#' @description Create a continuous time model and generate a one unit discrete time
#' prediction for auto and cross effects.
#' The function returns the auto and cross effects translated to a discrete time model,
#' and a change of t, dt = 1. These are the AR AND CL parameters in the CLPM.
#'
#' #'
#' @param data Must be a long format data frame.
#' @param id id, which is used to estimate subject specific effects.
#' @importFrom dplyr %>% select mutate row_number group_by ungroup
#' @importFrom tidyr pivot_longer
#' @importFrom readr parse_number
#'
#' @return A list; returns model and drift matrix
#' #'
#' @examples

#'
#' @export
#'
#'

estimateCTSEM <- function(data = dat_long,
                          id = "id",
                          manifest.var.names = c("x", "y"),
                          latent.var.names = c("x", "y"),
                          LAMBDA = diag(2),
                          type = "stanct",
                          time_points = 10,
                          cores = cores,
                          chains = chains,
                          model = NA,
                          ...) {
  my_model <- ctsem::ctModel(
    Tpoints = time_points,
    id = id,
    manifestNames = manifest.var.names,
    latentNames = latent.var.names,
    LAMBDA = LAMBDA,
    type = type,
    ...
  )
  fit <- ctsem::ctStanFit(
    datalong = data,
    iter = 3000,
    ctstanmodel = my_model,
    ...
  )
  summary(fit, timeinterval = 1)$parmatrices %>%
    filter(matrix == "dtDRIFT") %>%
    select(Mean) %>%
    t() %>%
    as.data.frame() -> drift
  names(drift) <- c("xx", "xy", "yx", "yy")
  return(list(parms = fit, drift = drift, model = my_model))
}
