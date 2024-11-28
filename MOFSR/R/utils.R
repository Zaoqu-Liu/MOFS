#' @noRd
zquiet <- function(..., messages = FALSE, cat_output = FALSE) {
  if (!cat_output) {
    sink(tempfile())
    on.exit(sink())
  }
  output <- if (messages) {
    eval(...)
  } else {
    suppressMessages(eval(...))
  }
  output
}

#' @importFrom dplyr %>%
#' @export
dplyr::`%>%`
