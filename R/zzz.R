.onAttach <- function(libname, pkgname) {
  ver <- utils::packageVersion(pkgname)
  packageStartupMessage(
    "Loading: ", pkgname, " (", ver, ")\n",
    "This package is BETA and you may encounter bugs.\n",
    "Please report any issues at:\n",
    "  https://github.com/crweber9874/crossLagR/issues"
  )
}

# NAMESPACE import order: R resolves conflicting generic names (coef, vcov,
# anova, nobs, filter, lag, ...) to whichever import() comes LAST, and roxygen2
# writes import() directives alphabetically. stats is therefore imported with
# targeted importFrom() calls rather than import(stats), so lavaan/dplyr/lme4
# keep their own methods for those generics. Add stats functions to the
# @importFrom in R/estimateRI.R or R/monteCarloRI.R as needed -- do not
# reintroduce a blanket @import stats.
