#' Transformation (stabilisation) to stationarity
#'
#' For a given mutation model `(M,p)`, transform `M` into another mutation
#' matrix `S` such that `S` is stationary with respect to `p`. Several methods
#' for doing this are described by Simonsson and Mostad (2016); only the "PM"
#' method is included here.
#'
#' These transformations may also be applied by setting `transform = "PM"` in
#' [mutationMatrix()] or [mutationModel()].
#'
#' For details about the transformation, see Simonsson and Mostad (2016).
#'
#' This function is a slightly optimised version of the `stabilize()` method in
#' the Familias R package.
#'
#' @param mutmat A square mutation matrix; typically a [mutationMatrix()] or
#'   [mutationModel()].
#' @param method A character string indicating the method to use. Currently only
#'   "PM" is implemented.
#' @param afreq A vector of allele frequencies. Extracted from `mutmat` if not
#'   provided.
#'
#' @return An object of the same class the input `mutmat`; either a matrix, a
#'   `mutationMatrix` or a `mutationModel`.
#'
#' @references Simonsson & Mostad (2016). *Stationary mutation models*. Forensic
#'   Sci. Int. Genet. 23:217–225. \doi{10.1016/j.fsigen.2016.04.005}
#'
#' @examples
#'
#' afreq = c(`1` = .2, `2` = .3, `3` = .5)
#' m = mutationMatrix("step", afreq = afreq, rate=0.1, rate2=0.01, range=0.1)
#' m
#' makeStationary(m, afreq = c(.3,.3,.4))
#'
#'
#' ### Example with full model (i.e., male and female)
#'
#' M = mutationModel("equal", afreq = afreq, rate = list(male=0.1, female=0.2))
#' M
#' makeStationary(M)
#'
#' @export
makeStationary = function(mutmat, afreq = NULL, method = "PM") {

  if(isMutationModel(mutmat))
    return(mapFullModel(mutmat, makeStationary, afreq = afreq, method = method))

  if(!is.matrix(mutmat))
    stop2("Argument `mutmat` must be either a `matrix`, a `mutationMatrix()` or a `mutationModel()`: ", class(mutmat)[1])

  method = match.arg(method)
  afreq = afreq %||% attr(mutmat, "afreq") %||% stop2("`afreq` must be provided")

  S = .stabilizePM(mutmat, afreq)
  newrate = mutRate(mutmat, afreq)

  newMutationMatrix(S, "custom", afreq = afreq, rate = newrate)
}


# Based on Familias::stabilize
.stabilizePM = function(M, p) {

  n = dim(M)[1]
  dM = diag(M)
  I = diag(n)

  R = 1 - sum(dM * p)

  X = t.default(M) - I
  X[n, ] = 1

  v = solve.default(X, c(rep(0, n-1), 1))
  d = R * v / ((1 - sum(dM * v)) * p)

  if(any(d * (1 - dM) > 1))
    stop2("PM stabilization doesn't exist.")

  (M - I) * d + I
}
