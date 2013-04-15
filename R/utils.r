#' @importFrom ape read.dna
#' @importFrom ape write.dna
#' @importFrom Biostrings writeXStringSet
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings readDNAMultipleAlignment
#' @import methods
setOldClass("DNAbin")


setAs("XStringSet", "DNAbin", function (from) {
  tmp <- tempfile()
  on.exit(unlink(tmp))
  writeXStringSet(from, filepath=tmp, format="fasta")
  to <- read.dna(tmp, format="fasta")
  to
})


setAs("DNAbin", "DNAStringSet", function (from) {
  tmp <- tempfile()
  on.exit(unlink(tmp))
  write.dna(from, file=tmp, format="fasta")
  to <- readDNAStringSet(tmp, "fasta")
  to
})


setAs("DNAbin", "DNAMultipleAlignment", function (from) {
  tmp <- tempfile()
  on.exit(unlink(tmp))
  write.dna(from, file=tmp, format="fasta")
  to <- readDNAMultipleAlignment(tmp, "fasta")
  to
})


setAs("XStringSet", "DNAMultipleAlignment", function (from) {
  tmp <- tempfile()
  on.exit(unlink(tmp))
  writeXStringSet(from, filepath=tmp, format="fasta")
  to <- readDNAMultipleAlignment(tmp, "fasta")
  to
})

