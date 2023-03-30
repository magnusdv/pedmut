#' Combine alleles in a mutation matrix
#'
#' Reduce a mutation matrix by combining a set of alleles into one "lump", if
#' this can be done without distorting the mutation process of the remaining
#' alleles. Such "allele lumping" can give dramatic efficiency improvements in
#' likelihood computations with multi-allelic markers, in cases where only some
#' of the alleles are observed in the pedigree.
#'
#' @param mutmod A `mutationModel` object, typically made with
#'   [mutationModel()].
#' @param mutmat A `mutationMatrix` object, typically made with
#'   [mutationMatrix()].
#' @param afreq A vector with frequency vector, of the same length as the size
#'   of `mutmat`. If not given, the `afreq` attribute of the matrix is used.
#' @param lump A nonempty subset of the alleles (i.e., the column names of
#'   `mutmat`), or a list of several such subsets.
#' @param check A logical indicating if lumpability should be checked before
#'   lumping. Default: TRUE.
#' @param labelSep ((For debugging) A character used to name lumps by pasting
#'   allele labels.
#'
#' @return A reduced mutation model. If the original matrix has dimensions
#'   \eqn{n\times n}{n*n}, the result will be \eqn{k\times k}{k*k}, where \eqn{k
#'   = n - length(lump) + 1}.
#'
#' @seealso [mutationModel()], [mutationMatrix()]
#'
#' @examples
#'
#'
#' ### Example 1: Lumping a mutation matrix
#' mat = mutationMatrix("eq", alleles = 1:5,
#'                      afreq = rep(0.2, 5), rate = 0.1)
#' mat
#'
#' # Lump alleles 3, 4 and 5
#' mat2 = lumpedMatrix(mat, lump = 3:5)
#' mat2
#'
#' # Example 2: Full model, proportional
#' mutrate = list(male = 0.1, female = 0.2)
#' mod = mutationModel("prop", alleles = 1:4,
#'                     rate = mutrate, afreq = c(.1,.2,.3,.4))
#' mod
#'
#' # Lump alleles 3 and 4
#' mod2 = lumpedModel(mod, lump = 3:4)
#' mod2
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

  lumpedFreq = NULL
  afreq = if(is.null(afreq)) attr(mutmat, 'afreq') else checkAfreq(afreq, alleles = als)

  if(length(lump) == 1) {
    lump = lump[[1]]
    keep = .mysetdiff(als, lump)
    N = length(keep) + 1

    newM = mutmat[c(keep, lump[1]), c(keep, lump[1]), drop = FALSE]
    newM[, N] = 1 - .rowSums(newM[, seq_len(N-1), drop = FALSE], N, N-1)
    colnames(newM)[N] = rownames(newM)[N] = lumpName

    if(!is.null(afreq)) {
      lump_idx = match(lump, als)
      keep_idx = match(keep, als)
      lumpedFreq = c(afreq[keep_idx], sum(afreq[lump_idx]))
      names(lumpedFreq) = colnames(newM)
    }
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
lumpedModel = function(mutmod, lump, afreq = NULL, check = TRUE) {

  if(!inherits(mutmod, "mutationModel"))
    stop2(sprintf("Expected the input to be a `mutationModel`, but got a: `%s`", class(mutmod)[1]))

  if(is.null(afreq))
    afreq = attr(mutmod$female, "afreq")

  lumpedF = lumpedMatrix(mutmod$female, lump = lump, afreq = afreq, check = check)

  sexeq = sexEqual(mutmod)
  if(sexeq)
    lumpedM = lumpedF
  else
    lumpedM = lumpedMatrix(mutmod$male, lump = lump, afreq = afreq, check = check)

  # Still lumpable?
  alwaysLmp = alwaysLumpable(lumpedF) && (sexeq || alwaysLumpable(lumpedM))

  # Create model object
  structure(list(female = lumpedF, male = lumpedM), sexEqual = sexeq,
            alwaysLumpable = alwaysLmp, class = "mutationModel")
}


prepLump = function(lump, alleles = NULL, labelSep = NULL) {
  if(!is.list(lump))
    lump = list(lump)

  # Remove empty and trivial lumps
  lump[lengths(lump) < 2] = NULL

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

