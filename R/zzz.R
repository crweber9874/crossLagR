.onAttach <- function(libname, pkgname) {
  ver <- utils::packageVersion(pkgname)
  packageStartupMessage(
    "Loading: ", pkgname, " (", ver, ")\n",
    "This package is BETA and you may encounter bugs.\n",
    "Please report any issues at:\n",
    "  https://github.com/crweber9874/crossLagR/issues"
  )
}
