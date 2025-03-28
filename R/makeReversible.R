#' Transformations to reversibility
#'
#' This function implements three methods for transforming a mutation matrix to
#' a reversible model. All methods are based on Metropolis-Hastings proposal
#' functions.
#'
#' @param mutmat A mutation matrix.
#' @param method A character indicating which transformation to use. Either "BA"
#'   (Barker), "MH" (Metropolis-Hastings) or "PR" (preserved rate).
#' @param afreq A vector of allele frequencies. Extracted from `mutmat` if not
#'   provided.
#' @param adjust Logical. If TRUE, the overall mutation rate is adjusted to
#'   preserve the original rate; see [adjustRate()]. Not relevant for the "PR"
#'   method, which by construction always preserves the overall rate.
#' @returns A reversible mutation matrix with the same allele frequencies.
#'
#' @examples
#' m = mutationMatrix("equal", afreq = c(a=0.2, b=0.3, c=0.5), rate = 0.2)
#' makeReversible(m)
#' makeReversible(m, adjust = TRUE)
#'
#' makeReversible(m, "MH", adjust = TRUE)
#' # makeReversible(m, "PR") # not well-defined!
#'
#' @export
makeReversible = function(mutmat, method = c("BA", "MH", "PR"), adjust = FALSE,
                          afreq = NULL) {

  method = match.arg(method)
  afreq = afreq %||% attr(mutmat, "afreq") %||% stop2("`afreq` must be provided")

  M = as.matrix(mutmat)
  pm = M * afreq
  pmT = t.default(pm)

  # Transform the matrix
  switch(method,
    MH = {
      R = M * pmin(1, pmT/pm)
      R[M == 0] = 0
    },
    BA = {
      R = M * pmT/(pm + pmT)
      R[M == 0] = 0
    },
    PR = {
      # Check if the model is well-defined (formula from paper)
      if(any(rowSums(pm) < afreq * (1 + diag(M))))
        stop2("The PR transformation is not well-defined for this model")
      R = (pm + pmT)/(2*afreq)
    }
  )

  # Modify diagonal so that row sums are 1
  diag(R) = 0
  diag(R) = 1 - rowSums(R)

  # Adjust rate if needed
  if((method != "PR") && adjust) {
    origRate = mutRate(M, afreq)
    R = adjustRate(R, newrate = origRate, afreq = afreq)
  }

  # Return as mutation matrix
  newMutationMatrix(R, "custom", afreq = afreq, rate = mutRate(R, afreq))
}
