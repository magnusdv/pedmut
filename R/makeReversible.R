#' Transformations to reversibility
#'
#' This function implements three methods for transforming a mutation model
#' `(M,p)` into a reversible one, `(R,p)`. All methods are based on Metropolis-
#' Hastings proposal functions.
#'
#' These transformations may also be applied through the `transform` argument of
#' [mutationMatrix()] and [mutationModel()].
#'
#' @param mutmat A [mutationMatrix()] or [mutationModel()].
#' @param method A character indicating which transformation to use. Either "BA"
#'   (Barker), "MH" (Metropolis-Hastings) or "PR" (preserved rate).
#' @param afreq A vector of allele frequencies. Extracted from `mutmat` if not
#'   provided.
#' @param adjust Logical. If TRUE (default), the overall mutation rate is
#'   adjusted to preserve the original rate; see [adjustRate()]. Not relevant
#'   for method "PR", which by construction always preserves the overall rate.
#'
#' @returns A reversible mutation matrix with the same allele frequencies.
#'
#' @examples
#' m = mutationMatrix("equal", afreq = c(a=0.2, b=0.3, c=0.5), rate = 0.2)
#' makeReversible(m, "BA")
#' makeReversible(m, "MH")
#' makeReversible(m, "PR")
#'
#' makeReversible(m, "BA", adjust = FALSE)  # rate differs!
#'
#' # Apply to full model with different female/male rates
#' mod = mutationModel("equal", afreq = c(a=0.2, b=0.3, c=0.5),
#'                     rate = list(female = 0.1, male = 0.2))
#' modR = makeReversible(mod)
#'
#' @export
makeReversible = function(mutmat, method = c("BA", "MH", "PR"), adjust = TRUE,
                          afreq = NULL) {

  if(isMutationModel(mutmat))
    return(mapFullModel(mutmat, makeReversible, method = method, adjust = adjust, afreq = afreq))

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
      if(any(colSums(pm) > afreq * (1 + 2*diag(M))))
        stop2("The PR transformation is not well-defined for this model")
      R = (pm + pmT)/(2*afreq)
    }
  )

  # Modify diagonal so that row sums are 1
  diag(R) = 0
  diag(R) = 1 - rowSums(R)

  # Original overall rate
  rate = mutRate(M, afreq)

  # Adjust rate of R if needed
  if(adjust && method != "PR") {
    R = adjustRate(R, newrate = rate, afreq = afreq)
    newrate = rate
  }
  else {
    newrate = mutRate(R, afreq)
  }

  # Return as mutation matrix
  newMutationMatrix(R, "custom", afreq = afreq, rate = newrate)
}

