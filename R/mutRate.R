#' Overall mutation rate
#'
#' Calculate the overall mutation rate at a locus, given a mutation model an a
#' set of allele frequencies.
#'
#' The mutation rate is found by the formula `1 - sum(diag(mut) * afreq)`.
#'
#' If `mut` is a [mutationModel()], the rate is calculated separately for the
#' male and female matrices.
#'
#' @param mut A [mutationModel()] or [mutationMatrix()].
#' @param afreq A vector of allele frequencies.
#'
#' @return A single number, or (if `mut` is a [mutationModel()] and the female
#'   and male rates differ) a list of two numbers, named "female" and "male".
#'
#' @examples
#' m = mutationMatrix("stepwise", alleles = 1:4, afreq = c(.1,.2,.3,.4),
#'                    rate = 0.01, rate2 = 1e-6, range = 0.1)
#' r = mutRate(m)
#'
#' stopifnot(all.equal(r, 0.01))
#' @export
mutRate = function(mut, afreq = NULL) {
  if(isMutationModel(mut)) {
    r = lapply(mut, function(m) mutRate(m, afreq))
    return(if(all.equal(r$male, r$female)) r$male else unlist(r))
  }

  afreq = afreq %||% attr(mut, "afreq") %||%
      stop2("Argument `afreq` is missing and not present as model attribute")

  1 - sum(diag(mut) * afreq)
}
