#' @section Package options:
#'
#' \emph{blastr} uses the following \code{\link{options}} to configure behaviour:
#'
#' \itemize{
#'   \item \code{blastr.blastdb.path}: Path to a local installation of the NCBI
#'   BLAST databases. By default this option is not set. I recommend setting it to
#'   a path of your choice in your .Rprofile file. Run \code{\link{update_blastdb}}
#'   to install these databases locally.
#' }
#' 
#' @docType package
#' @name blastr
NULL

.onLoad <- function(libname, pkgname) {
  ## set global options
  options(verbose = FALSE)
  op <- options()
  op.blastr <- list(
    blastr.blastdb.path = NULL
  )
  toset <- !(names(op.blastr) %in% names(op))
  if (any(toset)) {
    options(op.blastr[toset])
  }
  invisible()
}

.onUnload <- function(libpath) {
  options(verbose = FALSE)
}
