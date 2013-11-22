#' @section Package options:
#'
#' \emph{blastr} uses the following \code{\link{options}} to configure behaviour:
#'
#' \itemize{
#'   \item \code{blastr.blastdb.path}: Path to a local installation of the NCBI
#'   BLAST databases. Defaults to "$HOME/local/db/blast". You can override the
#'   default by setting this option in your .Rprofile file. Run\code{\link{update_blastdb}}
#'   to install these databases.
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
    blastr.blastdb.path = normalizePath("~/local/db/blast", mustWork=FALSE)
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
