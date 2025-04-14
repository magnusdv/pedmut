#' Stabilization of mutation matrix
#'
#' NB: REPLACED BY [makeStationary]. Produces a mutation matrix close to the
#' input `mutmat`, for which the given frequency vector is the stationary
#' distribution. Several methods for doing this are described by Simonsson and
#' Mostad (2016); only the "PM" method is included here.
#'
#' This function is based on, and reuses code from, the `stabilize()` method of
#' the Familias R package.
#'
#' @param mutmat A mutation matrix.
#' @param afreq A vector of allele frequencies
#' @param method Either "DP", "RM" or "PM". Currently only "PM" is implemented.
#' @param details A logical. If TRUE, the complete Familias output is included.
#'
#' @return An object of the same class the input `mutmat`; either a matrix, a
#'   `mutationMatrix` or a `mutationModel`.
#'
#' @author Petter Mostad, Thore Egeland, Ivar Simonsson, Magnus D. Vigeland
#' @references Simonsson, Mostad: Stationary Mutation models. (FSI: Genetics,
#'   2016).
#'
#' @examples
#'
#' afreq = c(.2, .3, .5)
#' m = mutationMatrix("stepwise", alleles = 1:3, afreq = afreq,
#'                    rate = 0.1, rate2 = 0.01, range = 0.1)
#' m
#' stabilize(m, afreq = c(.3,.3,.4))
#'
#'
#' ### Example with full model (i.e., male and female)
#'
#' M = mutationModel("stepwise", alleles = 1:3, afreq = afreq,
#'                    rate = list(male = 0.1, female = 0.2),
#'                    rate2 = 0.01, range = 0.1)
#' M
#' stabilize(M)
#'
#' @export
stabilize = function(mutmat, afreq = NULL, method = "PM", details = FALSE) {

  if(isMutationModel(mutmat)) {
    Pm = stabilize(mutmat$male, afreq = afreq, method = method)
    Pf = stabilize(mutmat$female, afreq = afreq, method = method)
    return(mutationModel(list(male = Pm, female = Pf)))
  }
  else if(isMutationMatrix(mutmat)) {
    if(is.null(afreq))
      afreq = attr(mutmat, "afreq")
    P = stabilize(as.matrix(mutmat), afreq, method)
    res = mutationMatrix("custom", matrix = P, afreq = afreq)
    return(res)
  }
  else if(!is.matrix(mutmat))
    stop2("Argument `mutmat` must be either a `matrix`, a `mutationMatrix()` or a `mutationModel()`: ", class(mutmat)[1])

  ### Variables used in Familias code below
  M = as.matrix(mutmat)
  pe = afreq %||% stop2("Argument `afreq` is missing")
  t = 1
  ###

  R = 1 - sum(diag(M) * pe)
  n = dim(M)[1]

  X = t.default(M) - diag(n)
  X[n, ] = rep(1, n)
  v = solve.default(X, c(rep(0, n - 1), 1))
  d = R * v/((1 - sum(diag(M) * v)) * pe)

  if (any(d * (1 - diag(M)) > t))
    stop2("PM stabilization doesn't exist.")

  P = diag(d) %*% (M - diag(n)) + diag(n)
  rownames(P) = colnames(P)

  if(details) {
    fratio = max(max(P/M), max(M/P))
    minS = min(diag(P))
    return(list(stabilized = P, fratio = fratio, mindiag = minS, error = ""))
  }

  P
}
