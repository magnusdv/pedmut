#' Mutation matrix
#'
#' Construct mutation matrices for pedigree likelihood computations.
#'
#' Descriptions of the models:
#'
#' * "custom" : Allows any mutation matrix to be provided by the user, in the `matrix`
#' parameter
#'
#' * "equal" :  All mutations equally likely; probability \eqn{1-rate} of no
#' mutation
#'
#' * "proportional" : Mutation probabilities are proportional to the target
#' allele frequencies
#'
#' * "random" : This produces a matrix of random numbers, where each row is
#' normalised so that it sums to 1
#'
#' @param model A string: either "custom", "equal", "proportional" or "random"
#' @param matrix When `model` is "custom", this must be a square matrix with
#'   nonnegative real entries and row sums equal to 1
#' @param alleles A character vector (or coercible to character) with allele
#'   labels. Required in all models, except "custom" if `matrix` has dimnames
#' @param afreq A numeric vector of allele frequencies. Required in model
#'   "proportional"
#' @param rate A number between 0 and 1. Required in models "equal" and
#'   "proportional"
#' @param seed A single number. Optional parameter in the "random" model, passed
#'   on to `set.seed()`
#' @param mutmat An object of class `mutationMatrix`
#'
#' @return A square matrix with entries in `[0, 1]`, with the allele labels as
#'   both colnames and rownames.
#'
#' @examples
#' mutationMatrix(alleles = 1:3, model = "equal", rate = 0.05)
#'
#' @importFrom stats runif
#' @export
mutationMatrix = function(model = c("custom", "equal", "proportional", "random"),
                          matrix = NULL, alleles = NULL, afreq = NULL,
                          rate = NULL, seed = NULL) {
  model = match.arg(model)

  if(model == "custom") {
    if(is.null(matrix))
      stop2('The "custom" model requires the argument `matrix` to be non-NULL')
    if(!is.matrix(matrix))
      stop2("Custom matrix must be a matrix, not a ", class(matrix))
    if(typeof(matrix) != "double")
      stop2("Custom matrix must be numeric, not ", typeof(matrix))
    if(nrow(matrix) != ncol(matrix))
      stop2("Custom matrix must be square, not ", nrow(matrix), "*", ncol(matrix))

    if(!is.null(mnames <- colnames(matrix))) {
      if(!identical(mnames, rownames(matrix)))
        stop2("Custom matrix must have identical colnames and rownames")

      # If alleles are given, make sure they match, and use them to permute
      if(!is.null(alleles)) {
        alleles = as.character(alleles)
        if(!setequal(alleles, mnames))
          stop2("Custom matrix names don't match the `alleles` argument")
        matrix = matrix[alleles, alleles]
      }
    }
    else {
      if(is.null(alleles))
        stop2("When custom matrix lacks names, the argument `alleles` cannot be NULL")
      if(length(alleles) != dim(matrix)[1L])
        stop2("Length of `alleles` must equal the dimension of `matrix`")
      dimnames(matrix) = list(alleles, alleles)
    }

    m = newMutationMatrix(matrix, model = model)
    return(validateMutationMatrix(m))  # must be done for custom models!
  }

  ## Check input for non-custom models
  if(!is.null(matrix))
    stop2("The `matrix` argument must be NULL when this model is specified")

  if(is.null(alleles))
    stop2("`alleles` cannot be NULL with this model")
  if(!is.null(afreq)) {
    if(length(afreq) != length(alleles))
      stop2("`afreq` must have the same length as `alleles`")
    if(round(sum(afreq), 3) != 1)
      stop2("Allele frequencies must sum to 1 after rounding to 3 decimals: ", sum(afreq))
  }
  if(model != "random") {
    if(is.null(rate))
      stop2("`rate` cannot be NULL with this model")
    if(!is_number(rate, minimum=0))
      stop2("`rate` must be a nonnegative number: ", rate)
  }

  ## Compute matrix according to model
  nall = length(alleles)
  if(nall < 2)
    return(NULL)
  mutmat = matrix(ncol = nall, nrow = nall, dimnames = list(alleles, alleles))

  if(model == "equal") {
    mutmat[] = rate/(nall - 1)
    diag(mutmat) = 1 - rate
  }
  else if(model == "proportional") {
    if(is.null(afreq))
      stop2("`afreq` cannot be NULL with this model")
    alpha = rate / sum(afreq * (1 - afreq))
    mutmat[] = (1 - alpha) * diag(nall) + alpha * rep(afreq, each = nall)
    if(any(mutmat > 1))
      stop2("Impossible mutation matrix; try reducing `rate`")
  }
  else if(model == "random") {
    if(!is.null(seed))
      set.seed(seed)
    mutmat[] = runif(nall^2, min = 0, max = 1)
    mutmat = mutmat / rowSums(mutmat)
  }

  newMutationMatrix(mutmat, model=model, afreq=afreq, rate=rate, seed=seed)
}

newMutationMatrix = function(mutmat, model = "custom", afreq = NULL,
                            rate = NULL, lumpedAlleles = NULL, seed = NULL) {
  if(!is.null(afreq)) {
    stationary = isStationary(mutmat, afreq)
    reversible = isReversible(mutmat, afreq)
  }
  else {
    stationary = reversible = NULL
  }
  structure(mutmat,
            model = model,
            afreq = afreq,
            rate = rate,
            stationary = stationary,
            reversible = reversible,
            lumpedAlleles = lumpedAlleles,
            seed = seed,
            class = "mutationMatrix")
}


#' @rdname mutationMatrix
#' @export
validateMutationMatrix = function(mutmat, alleles = NULL) {

  stopifnot(is.matrix(mutmat),
            is.numeric(mutmat),
            nrow(mutmat) == ncol(mutmat),
            setequal(colnames(mutmat), rownames(mutmat)),
            inherits(mutmat, "mutationMatrix"))

  if(!is.null(alleles)) stopifnot(setequal(rownames(mutmat), alleles))

  if(any(mutmat < 0))
    stop2("Negative entries found in mutation matrix: ", mutmat[mutmat < 0])
  if(any(mutmat > 1))
    stop2("Entries exceeding 1 found in mutation matrix: ", mutmat[mutmat > 1])

  rs = rowSums(mutmat)
  if (any(round(rs, 3) != 1)) {
    print(rowSums(mutmat))
    stop2("Rows which do not sum to 1 (after rounding to 3 decimal places): ", which(round(rs, 3) != 1))
  }

  mutmat
}


#' Mutation models
#'
#' A mutation model is a list of two mutation matrices, named "female" and
#' "male".
#'
#' @param female Either an object of class `mutationMatrix`, or a character
#'   string. In the latter case, it is passed on to [mutationMatrix()] together
#'   with any arguments in `...`.
#' @param male An object of class `mutationMatrix`. By default, this will be the
#'   same as the female mutation matrix.
#' @param ... Further arguments to [mutationMatrix()], e.g. `alleles`.
#' @param mutmod A `mutationModel` object
#' @param alleles A character vector with allele labels. (The validation method
#'   uses this to check that the matrices have appropriate dimnames.)
#'
#' @return An object of class `mutationModel`.
#'
#' @examples
#' mutationModel("eq", alleles = 1:2, rate = 0.1)
#'
#' @export
mutationModel = function(female, male = female, ...) {

  if(!identical(class(female), class(male)))
    stop2("Arguments `female` and `male` must be of the same class")

  if(class(female) == "mutationMatrix")
    mod = list(female = female, male = male)
  else if(is.character(female)) {
    mutmat = mutationMatrix(model = female, ...)
    mod = list(female = mutmat, male = mutmat)
  }
  structure(mod, sexEqual = identical(female, male), class = "mutationModel")
}

#' @rdname mutationModel
#' @export
validateMutationModel = function(mutmod, alleles = NULL) {

  stopifnot(is.list(mutmod),
            length(mutmod) == 2,
            setequal(names(mutmod), c("male", "female")),
            inherits(mutmod, "mutationModel"))

  validateMutationMatrix(mutmod$male, alleles = alleles)
  validateMutationMatrix(mutmod$female, alleles = alleles)
}


#' @export
print.mutationMatrix = function(x, includeMatrix = TRUE, includeAttrs = TRUE,
                                includeProperties = TRUE, ...) {

  if(includeMatrix) {
    print(format(x), quote=FALSE, right=TRUE)
  }

  attrs = attributes(x)
  afreq = attrs$afreq

  if(includeAttrs) {
    if(includeMatrix) cat("\n")

    model = attrs$model
    rate = attrs$rate
    seed = attrs$seed

    cat("Model:", model, "\n")
    cat("Rate:", if(!is.null(rate)) rate else NA, "\n")
    if(!is.null(afreq))
      cat("Frequencies:", toString(afreq), "\n")
    if(model == "random")
      cat("Seed: ", if(!is.null(seed)) seed else NA, "\n")
  }

  if(includeProperties) {
    if(includeAttrs) cat("\n")

    if(!is.null(afreq)) {
      cat("Stationary:", if(isStationary(x, afreq)) "Yes" else "No", "\n")
      cat("Reversible:", if(isReversible(x, afreq)) "Yes" else "No", "\n")
    }
    cat("Lumpable:", if(alwaysLumpable(x)) "Always" else "Not always", "\n")
  }
}

#' @export
print.mutationModel = function(x, ...) {
  if(attr(x, "sexEqual")) {
    cat("Unisex mutation matrix:\n")
    print(x$female, ...)
  }
  else {
    cat("Female mutation matrix:\n")
    print(x$female, ...)

    cat("\nMale mutation matrix:\n")
    print(x$male, ...)
  }
  invisible(x)
}
