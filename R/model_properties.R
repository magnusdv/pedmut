#' Mutation model properties
#'
#' Functions that check various properties of a mutation model, including
#' stationarity, reversibility and lumpability.
#'
#' The function `isBounded()` checks that a mutation model is *bounded* by the
#' allele frequencies, i.e., that `mutmat[i,j] <= afreq[j]` whenever `i` is not
#' equal to `j`.
#'
#' Lumpability is a property of a mutation model that allows aggregating alleles
#' into groups, or *lumps*, without changing the overall mutation process. The
#' functions `isLumpable()` and `alwaysLumpable()` checks lumpability using the
#' row-sum criterion given by Kemeny & Snell (1976). Note that lumping may be
#' possible even if the model is not generally lumpable; see [lumpMutSpecial()]
#' for details.
#'
#' For each of these functions, if `mutmat` is a `mutationModel` object, i.e.,
#' with male and female components, the output is TRUE if and only if both
#' components satisfy the property in question.
#'
#' @param mutmat A [mutationMatrix()] or a [mutationModel()].
#' @param afreq A vector with allele frequencies, of the same length as the size
#'   of `mutmat`.
#' @param lump A character vector containing a nonempty set of allele labels.
#'
#' @return Each of these functions returns TRUE of FALSE.
#'
#' @references
#' Kemeny & Snell (1976). \emph{Finite Markov Chains}. Springer.
#'
#' @examples
#'
#' # "proportional" models are stationary and reversible
#' afr = c(0.2, 0.3, 0.5)
#' m_prop = mutationMatrix(model = "prop", alleles = 1:3, afreq = afr, rate = 0.1)
#' stopifnot(isStationary(m_prop, afr), isReversible(m_prop, afr))
#'
#' # "equal" model is stationary and reversible only when freqs are equal
#' m_eq = mutationMatrix(model = "eq", alleles = 1:3, rate = 0.1)
#' stopifnot(isStationary(m_eq, rep(1/3, 3)), isReversible(m_eq, rep(1/3, 3)))
#' stopifnot(!isStationary(m_eq, afr), !isReversible(m_eq, afr))
#'
#' # "equal" and "proportional" models allow allele lumping
#' stopifnot(isLumpable(m_eq, lump = 1:2))
#' stopifnot(isLumpable(m_prop, lump = 1:2))
#'
#' # In fact lumpable for any allele subset
#' stopifnot(alwaysLumpable(m_eq), alwaysLumpable(m_prop))
#'
#' @name model_properties
NULL

#' @rdname model_properties
#' @export
isStationary = function(mutmat, afreq = NULL) {
  if(isMutationModel(mutmat)) {
    isStatF = isStationary(mutmat$female, afreq = afreq)
    isStatM = sexEqual(mutmat) || isStationary(mutmat$male, afreq = afreq)
    return(isStatF && isStatM)
  }

  if(is.null(afreq))
    afreq = attr(mutmat, "afreq") %||% stop2("Argument `afreq` is missing")

  prod = as.numeric(afreq %*% mutmat)
  tol = sqrt(.Machine$double.eps)
  all(abs(as.numeric(afreq) - prod) < tol)
}

#' @rdname model_properties
#' @export
isReversible = function(mutmat, afreq = NULL) {
  if(isMutationModel(mutmat)) {
    isRevF = isReversible(mutmat$female, afreq = afreq)
    isRevM = sexEqual(mutmat) || isReversible(mutmat$male, afreq = afreq)
    return(isRevF && isRevM)
  }

  if(is.null(afreq))
    afreq = attr(mutmat, "afreq") %||% stop2("Argument `afreq` is missing")

  pm = afreq * mutmat
  pmT = t.default(pm)
  tol = sqrt(.Machine$double.eps)
  all(abs(as.numeric(pm) - as.numeric(pmT)) < tol)
}


#' @rdname model_properties
#' @export
isBounded = function(mutmat, afreq = NULL) {
  if(isMutationModel(mutmat)) {
    isboundF = isBounded(mutmat$female, afreq = afreq)
    isboundM = sexEqual(mutmat) || isBounded(mutmat$male, afreq = afreq)
    return(isboundF & isboundM)
  }

  if(is.null(afreq))
    afreq = attr(mutmat, "afreq") %||% stop2("Argument `afreq` is missing")

  n = length(afreq)
  M = as.matrix(mutmat)
  lines = rep(FALSE, n)
  for (i in 1:n)
    lines[i] = all(M[-i, i] <= afreq[i])
  all(lines)
}


#' @rdname model_properties
#' @export
isLumpable = function(mutmat, lump) {

  if(isMutationModel(mutmat)) {

    if(isTRUE(attr(mutmat, 'alwaysLumpable')))
      return(TRUE)

    test = isLumpable(mutmat$female, lump)
    if(!sexEqual(mutmat))
      test = test || isLumpable(mutmat$male, lump)

    return(test)
  }

  if(is.null(mutmat))
    return(TRUE)

  als = colnames(mutmat)
  lump = prepLump(lump, als)
  N = length(lump)

  if(N == 0)
    return(TRUE)

  tol = sqrt(.Machine$double.eps)

  if(N == 1) {
    lump = lump[[1]]

    # Exhaustive lump (trivial)
    if(length(lump) == length(als))
      return(TRUE)

    # Row sum criterion (reduced to equalities of entries)
    y = mutmat[lump, .mysetdiff(als, lump), drop = FALSE]
    checks = abs(as.numeric(y) - rep(y[1, ], each = length(lump))) < tol
    return(all(checks))
  }

  ### Multiple lumps (of size > 1) ###

  # Alleles not involved in lumps
  singles = .mysetdiff(als, unlist(lump, use.names = FALSE))

  # Loop through all lumps
  for(i in seq_len(N)) {
    a = lump[[i]]

    # Check against singles
    y = mutmat[a, singles, drop = FALSE]
    checks = abs(as.numeric(y) - rep(y[1, ], each = length(a))) < tol
    if(!all(checks))
      return(FALSE)

    # Check row sums against other lumps
    for(b in lump[-i]) {
      z = mutmat[a, b, drop = FALSE]
      rs = .rowSums(z, length(a), length(b))
      if(!.equalNums(rs, tol))
        return(FALSE)
    }
  }

  return(TRUE)
}

#' @rdname model_properties
#' @export
alwaysLumpable = function(mutmat) {

  if(isMutationModel(mutmat)) {
    alwaysF = alwaysLumpable(mutmat$female)
    return(alwaysF && (sexEqual(mutmat) || alwaysLumpable(mutmat$male)))
  }

  N = dim(mutmat)[1L] %||% 0L  # allows NULL input

  # 2*2 and smaller: trivially lumpable
  if(N <= 2)
    return(TRUE)

  offdiag = mutmat[col(mutmat) != row(mutmat)]
  dim(offdiag) = c(N - 1, N)

  # Kemeny & Snell criterion
  tol = sqrt(.Machine$double.eps)
  all(abs(as.numeric(offdiag) - rep(offdiag[1,], each = N-1)) < tol)
}


#' Find the stationary frequency distribution
#'
#' Finds the stationary distribution of allele frequencies, if it exists, w.r.t. a given mutation matrix.
#' @param mutmat A mutation matrix.
#'
#' @return A vector of length `ncol(mutmat)`, or NULL.
#'
#' @examples
#'
#' m1 = mutationMatrix("equal", alleles = 1:4, rate = 0.1)
#' findStationary(m1)
#'
#' m2 = mutationMatrix("random", alleles = 1:3, seed = 123)
#' a = findStationary(m2)
#'
#' a %*% m2 - a  # check
#'
#' @export
findStationary = function(mutmat) {
  validateMutationMatrix(mutmat)
  n = dim(mutmat)[1L]
  P = diag(n) - mutmat
  A = rbind(t.default(P), rep(1, n))
  b = c(rep(0, n), 1)

  # The following is effectively `qr.solve(A, b)`, but catches singular A gracefully
  QR = qr.default(A)

  if(QR$rank != n) {
    message("No stationary distribution")
    return()
  }

  qr.coef(QR, b)
}
