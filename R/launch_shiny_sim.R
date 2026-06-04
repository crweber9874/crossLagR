#' @title launch_shiny_sim
#' @description Launch the bundled \code{crossLagR} Shiny simulator.
#'
#' @details
#' The app lives at \code{inst/shinySim/} inside the package source tree, so it
#' is available under \code{system.file("shinySim", package = "crossLagR")}
#' after \code{install.packages()} and during interactive development under
#' \code{devtools::load_all()}.
#'
#' The app's runtime depends on a number of packages that are listed in
#' \code{Suggests} (not \code{Imports}) so that users who only need the
#' modeling functions are not forced to install the full Shiny + bslib stack.
#' If any of them are missing this function stops with an actionable message
#' telling the user exactly what to install.
#'
#' @param launch.browser Logical. Open the app in a browser (default
#'   \code{TRUE}). Pass \code{FALSE} when running headless / on a server.
#' @param port Integer or \code{NULL}. Port to bind to. \code{NULL} (default)
#'   lets Shiny pick an open port.
#' @param host Character. Host to bind to. Default \code{"127.0.0.1"} (local
#'   only). Use \code{"0.0.0.0"} to expose to the network.
#' @param quiet Logical. Suppress the startup banner printed by Shiny. Default
#'   \code{FALSE}.
#' @param ... Additional arguments forwarded to \code{shiny::runApp()}.
#'
#' @return Called for its side effect of starting the Shiny app. Returns
#'   invisibly.
#'
#' @examples
#' \dontrun{
#' # Standard local launch:
#' launch_shiny_sim()
#'
#' # Headless on a fixed port (e.g. for a remote SSH-forwarded session):
#' launch_shiny_sim(launch.browser = FALSE, port = 4321)
#' }
#'
#' @export
launch_shiny_sim <- function(launch.browser = TRUE,
                             port = NULL,
                             host = "127.0.0.1",
                             quiet = FALSE,
                             ...) {

  ## --- Locate the app directory ----------------------------------------------
  app_dir <- system.file("shinySim", package = "crossLagR")
  if (!nzchar(app_dir) || !dir.exists(app_dir)) {
    stop("Shiny app files not found at system.file('shinySim', package = 'crossLagR'). ",
         "Reinstall crossLagR, or call devtools::load_all() from the package root.",
         call. = FALSE)
  }

  ## --- Verify runtime dependencies are available -----------------------------
  required_pkgs <- c(
    shiny         = "shiny",
    bslib         = "bslib",
    shinyWidgets  = "shinyWidgets",
    DT            = "DT",
    dplyr         = "dplyr",
    ggplot2       = "ggplot2",
    plotly        = "plotly",
    DiagrammeR    = "DiagrammeR",
    digest        = "digest",
    callr         = "callr",
    lavaan        = "lavaan"
  )
  installed <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
  missing   <- required_pkgs[!installed]
  if (length(missing)) {
    stop(
      "Cannot launch the Shiny simulator. The following package",
      if (length(missing) > 1) "s are" else " is",
      " missing:\n  ",
      paste(missing, collapse = ", "),
      "\n\nInstall with:\n  install.packages(c(",
      paste(shQuote(missing), collapse = ", "), "))",
      call. = FALSE
    )
  }

  ## --- Launch -----------------------------------------------------------------
  if (!quiet) {
    message("Launching CrossLagR: Simulations from: ", app_dir)
  }
  invisible(shiny::runApp(
    appDir         = app_dir,
    launch.browser = launch.browser,
    port           = port,
    host           = host,
    quiet          = quiet,
    ...
  ))
}
