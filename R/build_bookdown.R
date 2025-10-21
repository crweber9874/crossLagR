#' Build the package bookdown documentation
#'
#' @param output_format Output format (default: "bookdown::gitbook")
#' @export
build_book <- function(output_format = "bookdown::gitbook") {
  old_dir <- getwd()
  on.exit(setwd(old_dir), add = TRUE)

  # Move to book directory
  setwd(file.path(old_dir, "book"))

  # Build book
  bookdown::render_book("index.Rmd", output_format = output_format)

  message("Book built successfully!")
  message("View at: ", file.path(old_dir, "book/_book/index.html"))

  invisible(file.path(old_dir, "book/_book"))
}

#' Preview the package book in browser
#'
#' @export
preview_book <- function() {
  book_path <- file.path(getwd(), "book/_book/index.html")

  if (!file.exists(book_path)) {
    message("Building book first...")
    build_book()
  }

  utils::browseURL(book_path)
}
