#!/usr/bin/env Rscript

# Clear knitr cache artifacts before each render so figures are always rebuilt.
# Do not remove .quarto here: Quarto preview tracks project state in .quarto and
# treats mutations during pre-render as an invalid project change.
# Clear only Quarto's freeze snapshots, not knitr's chunk-level *_cache
# directories. The chunk cache is what makes slow MC chunks (e.g. the
# lavaan-driven bias landscape) tolerable on re-render -- nuking it would
# turn every render into a 10-minute fit-fest.
cache_dirs <- c(
  "_freeze",
  ".quarto/_freeze",
  ".quarto/project-cache"
)

cache_dirs <- unique(cache_dirs[file.exists(cache_dirs)])

if (length(cache_dirs) > 0) {
  unlink(cache_dirs, recursive = TRUE, force = TRUE)
  message("Cleared render cache directories:")
  message(paste0(" - ", cache_dirs, collapse = "\n"))
} else {
  message("No render cache directories found.")
}
