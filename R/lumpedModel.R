#' Combine alleles in a mutation matrix
#'
#' Reduce a mutation matrix by combining a set of alleles into one "lump", if
#' this can be done without distorting the mutation process of the remaining
#' alleles. Such "allele lumping" can give dramatic efficiency improvements in
#' likelihood computations with multi-allelic markers, in cases where only some
#' of the alleles are observed in the pedigree.
#'
#' @param mutmod A `mutationModel` object.
#' @param mutmat A `mutationMatrix` object.
#' @param afreq A vector with frequency vector, of the same length as the size
#'   of `mutmat`. If not given, the `afreq` attribute of the matrix is used.
#' @param lump A nonempty subset of the column names of `mutmat` (i.e., the allele
#'   labels).
#' @param check A logical indicating if lumpability should be checked before
#'   lumping. Default: TRUE.
#'
#' @return A reduced mutation model. If the original matrix has dimensions
#'   \eqn{n\times n}{n*n}, the result will be \eqn{k\times k}{k*k}, where \eqn{k
#'   = n - length(lump) + 1}.
#'
#' @examples
#' m = mutationMatrix("eq", alleles = 1:10, rate = 0.1)
#' afreq = rep(1/100, 100)
#'
#' # Suppose only alleles 1 and 2 are observed.
#' # The following lumped matrix is then equivalent to `m`:
#' mLump = lumpedMatrix(m, afreq = afreq, lump = 3:10)
#' mLump
#'
#' # Full model
#' mod = mutationModel("prop", alleles = 1:4, rate = 0.1, afreq = c(.1,.2,.3,.4))
#' mod2 = lumpedModel(mod, lump = 3:4)
#'
#' @export
lumpedMatrix = function(mutmat, lump, afreq = attr(mutmat, 'afreq'), check = TRUE) {
  als = colnames(mutmat)

  # If all alleles are lumped, return trivial mutationModel
  if(setequal(lump, als)) {
    newM = matrix(1, ncol=1, nrow=1, dimnames=list("lump", "lump"))
    mod = newMutationMatrix(newM, model = "trivial")
    return(mod)
  }

  if(check && !isLumpable(mutmat, lump))
    stop2("The model is not lumpable for this set of alleles")

  keep = .mysetdiff(als, lump)
  N = length(keep) + 1

  lump_idx = match(lump, als)
  keep_idx = match(keep, als)

  newM = mutmat[c(keep, lump[1]), c(keep, lump[1])]
  newM[, N] = 1 - .rowSums(newM[, 1:(N-1), drop = F], N, N-1)
  colnames(newM)[N] = rownames(newM)[N] = "lump"

  if(!is.null(afreq)) {
    lumpedFreq = c(afreq[keep_idx], sum(afreq[lump_idx]))
    names(lumpedFreq) = colnames(newM)
  }
  else {
    lumpedFreq = NULL
  }

  newMutationMatrix(newM, afreq = lumpedFreq, lumpedAlleles = lump,
                   model = attr(mutmat, 'model'), rate = attr(mutmat, 'rate'),
                   seed = attr(mutmat, 'seed'))
}


#' @rdname lumpedMatrix
#' @export
lumpedModel = function(mutmod, lump, afreq = NULL, check = TRUE) {

  if(is.null(afreq))
    afreq = attr(mutmod$female, "afreq")

  sexeq = sexEqual(mutmod)

  lumpedF = lumpedMatrix(mutmod$female, lump = lump, afreq = afreq, check = check)

  if(sexeq)
    lumpedM = lumpedF
  else
    lumpedM = lumpedMatrix(mutmod$male, lump = lump, afreq = afreq, check = check)

  # Still lumpable?
  lmp = alwaysLumpable(lumpedF) && (sexeq || alwaysLumpable(lumpedM))

  # Create model object
  structure(list(female = lumpedF, male = lumpedM), sexEqual = sexeq,
            alwaysLumpable = lmp, class = "mutationModel")
}
