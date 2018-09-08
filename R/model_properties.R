#' Mutation model properties
#'
#' Functions for checking various properties of a mutation model, including
#' stationarity, reversibility and lumpability.
#'
#' @param mutmat A mutation matrix
#' @param afreq A vector with frequency vector, of the same length as the size
#'   of `mutmat`
#' @param lump A nonempty subset of the colnames of `mutmat` (i.e. the allele
#'   labels)
#' @return Each of these functions returns TRUE of FALSE
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
#' # "equal" and "proportional" models are always lumpable
#' stopifnot(isLumpable(m_eq, lump = 1:2))
#' stopifnot(isLumpable(m_prop, lump = 1:2))
#'
#' @name model_properties
NULL

#' @rdname model_properties
#' @export
isStationary = function(mutmat, afreq) {
  prod = as.numeric(afreq %*% mutmat)
  test = all.equal.numeric(prod, afreq, check.attributes = F)
  isTRUE(test)
}

#' @rdname model_properties
#' @export
isReversible = function(mutmat, afreq) {
  pm = afreq * mutmat
  symmetric = all.equal.numeric(pm, t.default(pm))
  isTRUE(symmetric)
}

#' @rdname model_properties
#' @export
isLumpable = function(mutmat, lump) {
  alleles = colnames(mutmat)
  stopifnot(!anyDuplicated(lump), all(lump %in% alleles))
  if(length(lump) == length(alleles))
    return(TRUE)
  y = mutmat[lump, setdiff(alleles, lump), drop = F]
  test = all.equal(as.numeric(y),
                   rep(y[1, ], each = length(lump)),
                   check.attributes = F)
  isTRUE(test)
}

alwaysLumpable = function(mutmat) {
  N = dim(mutmat)[1L]
  if(N == 1) return(FALSE)
  offdiag = mutmat[col(mutmat) != row(mutmat)]
  dim(offdiag) = c(N - 1, N)
  test = all.equal(as.numeric(offdiag),
                   rep(offdiag[1,], each = N-1),
                   check.attributes = F)
  isTRUE(test)
}
