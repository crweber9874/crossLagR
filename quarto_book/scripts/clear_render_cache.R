#!/usr/bin/env Rscript

# Clear knitr cache artifacts before each render so figures are always rebuilt.
# Do not remove .quarto here: Quarto preview tracks project state in .quarto and
# treats mutations during pre-render as an invalid project change.
cache_dirs <- c(
  "_freeze",
  ".quarto/_freeze",
  ".quarto/project-cache",
  grep(
    "_cache$",
    list.dirs(".", recursive = TRUE, full.names = TRUE),
    value = TRUE
  )
)

cache_dirs <- unique(cache_dirs[file.exists(cache_dirs)])

if (length(cache_dirs) > 0) {
  unlink(cache_dirs, recursive = TRUE, force = TRUE)
  message("Cleared render cache directories:")
  message(paste0(" - ", cache_dirs, collapse = "\n"))
} else {
  message("No render cache directories found.")
}
