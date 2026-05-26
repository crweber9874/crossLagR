## Sourced by Shiny once at app startup, before ui.R / server.R.
## Use this file to initialize the package and any objects the app needs.

## --- Initialize the crossLagR package ---
if (requireNamespace("crossLagR", quietly = TRUE)) {
  library(crossLagR)
} else if (requireNamespace("devtools", quietly = TRUE)) {
  ## Dev mode: app launched from the package root with the package not installed
  pkg_root <- normalizePath(file.path(getwd(), ".."), mustWork = FALSE)
  if (file.exists(file.path(pkg_root, "DESCRIPTION"))) {
    devtools::load_all(pkg_root, quiet = TRUE)
  } else {
    stop("crossLagR is not installed and DESCRIPTION not found at ", pkg_root,
         ". Install the package or launch the app from the package root.")
  }
} else {
  stop("crossLagR is not installed and devtools is unavailable; install one of them.")
}

## --- Sanity check that the simulation entry point is in scope ---
if (!exists("run_mc_sims", mode = "function")) {
  stop("run_mc_sims() not found after loading crossLagR. ",
       "Check the package install / load_all path.")
}
