#' Combine alleles in a mutation matrix
#'
#' Reduce a mutation matrix by combining a set of alleles into one *lump*, if
#' this can be done without distorting the mutation process of the remaining
#' alleles. Such allele lumping can give dramatic efficiency improvements in
#' likelihood computations with multi-allelic markers, in cases where only some
#' of the alleles are observed in the pedigree.
#'
#' The lumping implemented in this function is based on the Markov chain lumping
#' theory by Kemeny & Snell (1976). For other, specialised lumping, see
#' [lumpMutSpecial()].
#'
#' @param mutmat A `mutationMatrix` object, typically made with
#'   [mutationMatrix()].
#' @param mutmod A `mutationModel` object, typically made with
#'   [mutationModel()].
#' @param afreq A vector with allele frequencies, of the same length as the size
#'   of `mutmat`. Extracted from the model if not given.
#' @param lump A vector containing the alleles to be lumped together, or a list
#'   of several such vectors.
#' @param check A logical indicating if lumpability (i.e., the row-sum criterium
#'   of Kemeny & Snell) should be checked before lumping. Default: TRUE.
#' @param labelSep A character used to name lumps by pasting allele labels. (For
#'   debugging.)
#'
#' @return A reduced mutation model. If the original matrix has dimensions
#'   \eqn{n\times n}{n*n}, the result will be \eqn{k\times k}{k*k}, where \eqn{k
#'   = n - length(lump) + 1}.
#'
#' @seealso [lumpMutSpecial()].
#'
#' @references Kemeny & Snell (1976). \emph{Finite Markov Chains}. Springer.

#' @examples
#'
#' af = c(.1, .2, .3, .4)
#' names(af) = 1:4
#'
#' ### Example 1: Lumping a mutation matrix
#' mat = mutationMatrix("eq", afreq = af, rate = 0.1)
#' mat
#'
#' # Lump
#' lumpedMatrix(mat, lump = 3:4)
#' lumpedMatrix(mat, lump = 2:4)
#'
#' # Example 2: Full model, proportional
#' mutrate = list(male = 0.1, female = 0.2)
#' mod = mutationModel("prop", afreq = af, rate = mutrate)
#' mod
#'
#' # Lump
#' lumpedModel(mod, lump = 2:4)
#'
#' @export
lumpedMatrix = function(mutmat, lump, afreq = NULL, check = TRUE, labelSep = NULL) {

  if(!inherits(mutmat, "mutationMatrix"))
    stop2(sprintf("Expected the input to be a `mutationMatrix`, but got a: `%s`", class(mutmat)[1]))

  als = colnames(mutmat)
  lump = prepLump(lump, alleles = als, labelSep = labelSep)
  lumpName = names(lump)

  # If exhaustive lumped, return trivial mutationModel
  if(setequal(lump[[1]], als)) {
    newM = matrix(1, ncol=1, nrow=1, dimnames=list("lump", "lump"))
    mod = newMutationMatrix(newM, model = "trivial")
    return(mod)
  }

  if(check && !isLumpable(mutmat, lump))
    stop2("The model is not lumpable for this set of alleles")

  afreq = if(is.null(afreq)) attr(mutmat, 'afreq') else checkAfreq(afreq, alleles = als)
  lumpedFreq = NULL

  if(length(lump) == 1) {
    lump = lump[[1]]
    keep = .mysetdiff(als, lump)
    newAls = c(keep, lumpName)
    N = length(keep) + 1

    newM = matrix(0, N, N, dimnames = list(newAls, newAls))
    newM[, -N] = mutmat[c(keep, lump[1]), keep]
    newM[, N] = 1 - .rowSums(newM, N, N)

    if(!is.null(afreq))
      lumpedFreq = c(afreq[keep], lump = sum(afreq[lump]))
  }
  else {
    # Multiple nontrivial lumps!

    singles = .mysetdiff(als, unlist(lump, use.names = FALSE))
    lmp1 = vapply(lump, function(lmp) lmp[1], character(1))
    newAls = c(singles, lmp1)
    N = length(newAls)

    singleCols = mutmat[newAls, singles, drop = FALSE]
    lumpedCols = vapply(lump, function(lmp)
      .rowSums(mutmat[newAls, lmp, drop = FALSE], N, length(lmp)), FUN.VALUE = numeric(N))

    newM = cbind(singleCols, lumpedCols)
    rownames(newM) = colnames(newM)

    if(!is.null(afreq)) {
      lmpFrqs = vapply(lump, function(lmp) sum(afreq[match(lmp, als)]), numeric(1))
      lumpedFreq = c(afreq[match(singles, als)], lmpFrqs)
      names(lumpedFreq) = colnames(newM)
    }
  }

  newMutationMatrix(newM, afreq = lumpedFreq, lumpedAlleles = lump,
                   model = attr(mutmat, 'model'), rate = attr(mutmat, 'rate'),
                   seed = attr(mutmat, 'seed'))
}


#' @rdname lumpedMatrix
#' @export
lumpedModel = function(mutmod, lump, afreq = attr(mutmod, "afreq"), check = TRUE) {
  mapFullModel(mutmod, lumpedMatrix, lump = lump, afreq = afreq, check = check)
}


prepLump = function(lump, alleles = NULL, labelSep = NULL) {
  if(!is.list(lump))
    lump = list(lump)

  # Remove empty and trivial lumps
  lump[lengths(lump) < 2] = NULL

  if(!length(lump))
    return(lump)

  # Convert each lump to character vector
  lump = lapply(lump, as.character)

  n = length(lump)
  unl = unlist(lump, use.names = FALSE)

  # Checks
  if(dup <- anyDuplicated(unl))
    stop2("Duplicated entry in `lump`: ", unl[dup])

  if(!is.null(alleles) && !all(unl %in% alleles))
    stop2("Alleles not found in mutation matrix: ", .mysetdiff(unl, alleles))

  # Name lumps
  if(!is.null(labelSep))
    names(lump) = vapply(lump, function(v) paste(v, collapse = labelSep), character(1))

  if(is.null(names(lump))) {
    nms = if(n == 1) "lump" else paste0("lump", seq_len(n))
    if(any(nms %in% alleles))
      nms = make.unique(c(alleles, nms), sep = "")[-seq_along(alleles)]
    names(lump) = nms
  }

  lump
}

