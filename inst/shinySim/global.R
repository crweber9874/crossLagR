## Sourced by Shiny once at app startup, before ui.R / server.R.

## ---- Load crossLagR ---------------------------------------------------------
## Prefer the in-tree development version when DESCRIPTION is present, since
## the Shiny app uses functions (monteCarloLavaan, simFromSyntax,
## populate_unified_labels) added recently. If the installed crossLagR is
## stale, devtools::load_all() the package root instead so the app always
## sees the latest exports.

required_fns <- c("run_mc_sims", "monteCarloLavaan", "simFromSyntax",
                  "populate_unified_labels",
                  "estimateCLPM", "estimateRICLPM", "estimateALT",
                  "estimateLGM", "estimateLCMSR", "estimateLChange",
                  "estimateBollen_and_Brand", "estimateTSO",
                  "simCLPM", "simRICLPM", "simCLPMu", "simLChange")

has_all_fns <- function() {
  all(vapply(required_fns, exists, logical(1), mode = "function"))
}

## Find the package root: shinySim/ may sit at <pkg_root>/shinySim or
## <pkg_root>/shiny/shinySim. Walk upward until a DESCRIPTION is found.
find_pkg_root <- function(start = getwd(), max_up = 4) {
  d <- normalizePath(start, mustWork = FALSE)
  for (i in 0:max_up) {
    if (file.exists(file.path(d, "DESCRIPTION"))) return(d)
    d <- normalizePath(file.path(d, ".."), mustWork = FALSE)
  }
  NULL
}

pkg_root <- find_pkg_root()

if (!is.null(pkg_root) && requireNamespace("devtools", quietly = TRUE)) {
  ## Dev mode: always load_all the in-tree package so updates are picked up.
  suppressMessages(devtools::load_all(pkg_root, quiet = TRUE))
} else if (requireNamespace("crossLagR", quietly = TRUE)) {
  suppressMessages(library(crossLagR))
} else {
  stop("crossLagR is not installed and no in-tree DESCRIPTION was found near ",
       getwd(), ". Install the package or launch from the package root.")
}

if (!has_all_fns()) {
  ## Installed version is loaded but missing recent functions; force a
  ## load_all() from the in-tree package as a last resort.
  if (!is.null(pkg_root) && requireNamespace("devtools", quietly = TRUE)) {
    suppressMessages(devtools::load_all(pkg_root, quiet = TRUE))
  }
}

missing_fns <- required_fns[!vapply(required_fns,
                                    exists, logical(1), mode = "function")]
if (length(missing_fns)) {
  stop("crossLagR is loaded but the following functions are still missing: ",
       paste(missing_fns, collapse = ", "),
       ".\nIf you have an installed version of crossLagR shadowing the dev tree, run ",
       "remove.packages('crossLagR') and relaunch, or call devtools::load_all() ",
       "from '", pkg_root %||% "<unknown>", "' before launching the app.")
}

## ---- Required packages ------------------------------------------------------
required_pkgs <- c("digest", "DiagrammeR", "lavaan", "plotly", "ggplot2", "DT",
                   "bslib", "shinyWidgets", "dplyr", "callr", "tidyr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs,
                                      requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Missing required packages for the Shiny app: ",
       paste(missing_pkgs, collapse = ", "),
       ".\nInstall with: install.packages(c(",
       paste(shQuote(missing_pkgs), collapse = ", "), "))")
}

## ---- App helpers ------------------------------------------------------------
## DAG generators (mirror the quarto-book diagrams).
dags_path <- file.path(getwd(), "dags.R")
if (!file.exists(dags_path)) {
  ## App launched from package root rather than shinySim/ — locate dags.R.
  candidate <- file.path("shinySim", "dags.R")
  if (file.exists(candidate)) dags_path <- candidate
}
if (!file.exists(dags_path)) {
  stop("dags.R helper not found relative to working directory ", getwd(), ".")
}
source(dags_path, local = FALSE)

`%||%` <- function(a, b) if (!is.null(a)) a else b

## ---- Shared UI/server constants --------------------------------------------
## Defined here so both ui.R and server.R can reference them.

## Estimators available in the book + crossLagR package.
## (Allison-Chamberlain "FI" intentionally omitted.)
estimator_choices <- c(
  "CLPM"     = "CLPM",
  "RI-CLPM"  = "RICLPM",
  "ALT (Autoregressive Latent Trajectory)" = "ALT",
  "LGM (Latent Growth Model)"              = "LGM",
  "LCM-SR / GCLM"                          = "LCMSR",
  "LCS (Latent Change Score)"              = "LCHANGE",
  "Bollen-Brand (Dynamic Panel FE)"        = "BB",
  "TSO (Trait-State-Occasion)"             = "TSO",
  "CTSEM (Continuous-Time SEM)"            = "CTSEM",
  "RI (lmer random intercepts)"            = "RI",
  "OLS (wide-format regression)"           = "OLS"
)

## DGPs: one entry per data-generating model. Each value is the model name
## the backend dispatch understands. The dispatch (in server.R + mcSim.R)
## picks the most robust simulation path automatically:
##
##   - For models with both a hand-written sim*() function AND a syntax-route
##     path (CLPM, RI-CLPM, LCS), prefer the syntax route — it's more tolerant
##     at extreme ICCs and uses the estimator's own labels (no parameter-name
##     translation needed).
##   - For ALT, LGM, LCM-SR, Bollen-Brand, TSO, the syntax route is the only
##     option — those models have no bespoke sim* function.
##   - "CLPM + unmeasured confounder" is the only entry that REQUIRES the
##     parameterized sim (simCLPMu), because no estimator defines a confounder.
dgp_choices <- c(
  "CLPM"                          = "CLPM",
  "RI-CLPM"                       = "RICLPM",
  "ALT"                           = "ALT",
  "LGM"                           = "LGM",
  "LCM-SR"                        = "LCMSR",
  "LCS (Latent Change Score)"     = "LCHANGE",
  "Bollen-Brand"                  = "BB",
  "TSO"                           = "TSO",
  "CLPM + unmeasured confounder"  = "clpmu"
)
