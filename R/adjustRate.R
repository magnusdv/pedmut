#' Adjust the overall mutation rate of a model
#'
#' Adjusts the overall mutation rate of a model by scaling the off-diagonal
#' matrix entries.
#'
#' The adjusted matrix is calculated as `a * M + (1-a) * I`, where `M` is the
#' original matrix, `a = newrate/rate`, and `I` is the identity matrix.
#'
#' The maximum allowed value of `newrate` (to avoid negative values in the
#' adjusted matrix) is `rate/(1 - m))`, where `m` is the smallest diagonal
#' element in the original matrix.
#'
#' @param mutmat A mutation matrix with nonzero mutation overall rate.
#' @param newrate The new overall mutation rate.
#' @param afreq The allele frequencies. Extracted from the mutation matrix if
#'   not provided.
#' @param rate The current overall mutation rate. Calculated from the input if
#'   not provided.
#'
#' @seealso [mutRate()]
#'
#' @returns A new mutation matrix with the adjusted rate.
#'
#' @examples
#' m = mutationMatrix("equal", afreq = c(a=0.2, b=0.3, c=0.5), rate = 0.2)
#' m
#' adjustRate(m, 0.4)
#'
#' @export
adjustRate = function(mutmat, newrate, afreq = NULL, rate = NULL) {
  afreq = afreq %||% attr(mutmat, "afreq") %||% stop2("`afreq` not found")
  rate = rate %||% mutRate(mutmat, afreq)
  if(rate < sqrt(.Machine$double.eps))
    stop2("Cannot adjust rate for a model with rate 0")

  M = as.matrix(mutmat)
  a = newrate / rate

  # Check if within theoretical limit
  dg = diag(M)
  if(a > 1/(1 - min(dg)))
     stop2("Maximum adjusted rate is ", rate/(1 - min(dg)))

  newM = a * M + (1-a) * diag(length(afreq))
  newMutationMatrix(newM, model = "custom", afreq = afreq)
}
