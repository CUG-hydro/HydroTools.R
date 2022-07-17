#' @keywords internal
#' @import magrittr
#' @importFrom graphics abline axis legend lines par points
#' @importFrom methods as
#' @importFrom stats convolve cor.test median qf rnorm runif sd nls predict
#' @importFrom graphics text
#' @importFrom stats aggregate setNames
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL

.onLoad <- function(libname, pkgname) {
  
  if (getRversion() >= "2.15.1") {
    utils::globalVariables(
      c(
        ".", "Eeq", "Evp", "slope"
      )
    )
  }
}
