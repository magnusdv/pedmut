#' Transformations to reversibility
#'
#' This function implements three methods for transforming a mutation matrix to
#' a reversible model. All methods are based on Metropolis-Hastings (MH)
#' proposal functions.
#'
#' @param mutmat A mutation matrix.
#' @param method A character indicating which transformation to use. Either "BA"
#'   (Barker), "MH" (Metropolis-Hastings) or "PR" (preserved rate).
#' @param afreq A vector of allele frequencies. Extracted from `mutmat` if not
#'   provided.
#' @param adjust Logical. If TRUE, the mutation rate is adjusted to preserve the
#'   original rate; see [adjustRate()].
#' @returns A reversible mutation matrix with the same allele frequencies.
#'
#' @examples
#' m = mutationMatrix("equal", afreq = c(a=0.2, b=0.3, c=0.5), rate = 0.2)
#' makeReversible(m)
#' makeReversible(m, adjust = TRUE)
#' @export
makeReversible = function(mutmat, method = c("BA", "MH", "PR"), adjust = FALSE,
                          afreq = NULL) {
  afreq = afreq %||% attr(mutmat, "afreq") %||% stop2("`afreq` must be provided")

  M = as.matrix(mutmat)
  pm = M * afreq
  pmT = t.default(pm)

  switch(match.arg(method),
    MH = {
      R = M * pmin(1, pmT/pm)
      R[M == 0] = 0
    },
    BA = {
      R = M * pmT/(pm + pmT)
      R[M == 0] = 0
    },
    PR = {
      R = (pm + pmT)/(2*afreq)
    }
  )

  diag(R) = 0
  diag(R) = 1 - rowSums(R)

  if(adjust) {
    origRate = mutRate(M, afreq)
    R = adjustRate(R, newrate = origRate, afreq = afreq)
  }

  newMutationMatrix(R, "custom", afreq = afreq, rate = mutRate(R, afreq))
}
