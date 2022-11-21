#' Mutation models
#'
#' Constructor for the class `mutationModel`. An object of this class is
#' essentially a list of two mutation matrices, named "female" and "male".
#'
#' @param model Either:
#'
#'   * a `mutationModel` object (returned unchanged after validation)
#'
#'   * a single `mutationMatrix` object (will be applied to both genders)
#'
#'   * a list of two `mutationMatrix` objects, named "female" and "male"
#'
#'   * a single model name (see [mutationMatrix()] for valid options)
#'
#'   * a list of two model names, named "female" and "male"
#'
#' @param alleles A character vector with allele labels; passed on to
#'   [mutationMatrix()].
#' @param afreq A numeric vector of allele frequencies; passed on to
#'   [mutationMatrix()].
#' @param matrix A matrix, or a list of two (named "female" and "male")
#' @param rate A numeric mutation rate, or a list of two (named "female" and
#'   "male")
#' @param rate2 A numeric mutation rate, or a list of two (named "female" and
#'   "male"). Required in the "stepwise" model; see [mutationMatrix()] for
#'   details.
#' @param range A positive number, or a list of two (named "female" and "male").
#'   Required in the "stepwise" model; see [mutationMatrix()] for details.
#' @param seed An integer, or a list of two (named "female" and "male").
#' @param validate A logical, by default TRUE.
#' @param mutmod A `mutationModel` object.
#'
#' @return An object of class `mutationModel`. This is a list of two
#'   `mutationMatrix` objects, named "female" and "male", and the following
#'   attributes:
#'
#'   * `sexEqual` : TRUE if both genders have identical models, otherwise FALSE
#'
#'   * `alwaysLumpable` : TRUE if both genders have models that are lumpable for
#'   any allele subset, otherwise FALSE
#'
#' @examples
#' # "Equal" model, same parameters for both genders
#' M1 = mutationModel("eq", alleles = 1:2, rate = 0.1)
#' M1
#'
#' # Different mutation rates
#' M2 = mutationModel("eq", alleles = 1:2, rate = list(male = 0.1, female = 0.01))
#' M2
#'
#' stopifnot(identical(M1$male, M1$female), identical(M2$male, M1$male))
#'
#' # A custom mutation matrix:
#' mat = matrix(c(0,0,1,1), ncol = 2, dimnames = list(1:2, 1:2))
#' M3 = mutationModel(model = "custom", matrix = mat)
#'
#' # Under the hood arguments are passed to `mutationMatrix()`.
#' # Alternatively, this can be done explicitly in the `model` argument
#' M4 = mutationModel(model = mutationMatrix("custom", matrix = mat))
#'
#' stopifnot(identical(M3, M4))
#'
#' # The latter strategy is needed e.g. in pedtools::marker(), which gives the
#' # user access to `model`, but not `matrix`.
#'
#' @export
mutationModel = function(model, alleles = NULL, afreq = NULL, matrix = NULL,
                         rate = NULL, rate2 = NULL, range = NULL, seed = NULL,
                         validate = TRUE) {

  if(isMutationModel(model)) {
    mod = enforceAlleleOrder(model, alleles)
  }
  else if(is.list(model) &&
          length(model) == 2 &&
          setequal(names(model), c("female", "male")) &&
          all(sapply(model, isMutationMatrix))) {
    mod = list(female = enforceAlleleOrder(model$female, alleles),
               male = enforceAlleleOrder(model$male, alleles))
  }
  else if(isMutationMatrix(model)) {
    model = enforceAlleleOrder(model, alleles)
    mod = list(female = model, male = model)
  }
  else if((is.character(model) && length(model) == 1) || is.list(model) && length(model) == 2) {
    # Build models from parameters
    model = validateSingleInput(model, "character")
    matrix = validateMatrixInput(matrix)
    rate = validateSingleInput(rate, "numeric")
    rate2 = validateSingleInput(rate2, "numeric")
    range = validateSingleInput(range, "numeric")
    seed = validateSingleInput(seed, "numeric")

    female = mutationMatrix(model = model$female, alleles = alleles,
                            afreq = afreq, matrix = matrix$female,
                            rate = rate$female, rate2 = rate2$female,
                            range= range$female, seed = seed$female)

    male = mutationMatrix(model = model$male, alleles = alleles,
                          afreq = afreq, matrix = matrix$male,
                          rate = rate$male, rate2 = rate2$male,
                          range = range$male, seed = seed$male)

    mod = list(female = female, male = male)
  }
  else
    stop2("`model` must be either:\n",
          " * an object of class `mutationModel`\n",
          " * an object of class `mutationMatrix`` (applied to both genders)`\n",
          " * a list of two `mutationMatrix` objects, named 'female' and 'male'\n",
          " * a character string (see ?mutationModel for valid options)\n",
          " * a list of two characters, named 'female' and 'male'")

  sexEqual = identical(mod$female, mod$male)
  lumpable = alwaysLumpable(mod$female) && (sexEqual || alwaysLumpable(mod$male))

  mutmod = structure(mod, sexEqual = sexEqual, alwaysLumpable = lumpable,
                     class = "mutationModel")

  if(validate)
    validateMutationModel(mutmod)

  mutmod
}

validateSingleInput = function(x, mode) {
  if(is.null(x))
    return(x)

  if(is.vector(x, mode = mode) && length(x) == 1)
    return(list(female = x, male = x))

  if(is.list(x) && length(x) == 2 && setequal(names(x), c("female", "male")))
    return(x)

  stop2("Argument `", deparse(substitute(x)), '` must be either\n',
        " * a single ", mode, "\n",
        ' * a list of length 2, named "female" and "male"')
}

validateMatrixInput = function(x) {
  if(is.null(x))
    return(x)

  if(is.matrix(x))
    return(list(female = x, male = x))

  if(is.list(x) && length(x) == 2 && setequal(names(x), c("female", "male")))
    return(x)

  stop2("Argument `", deparse(substitute(x)), '` must be either\n',
        " * a single matrix\n",
        ' * a list of length 2 matrices, named "female" and "male"')
}

#' @rdname mutationModel
#' @export
validateMutationModel = function(mutmod, alleles = NULL) {
  if(!is.list(mutmod) || length(mutmod) != 2 || !setequal(names(mutmod), c("male", "female"))) {
    stop2("`mutmod` must be a list with elements 'male' and 'female'")
  }

  male = mutmod$male
  female = mutmod$female

  if(is.null(alleles))
    alleles = colnames(male)

  validateMutationMatrix(male, alleles = alleles)

  sexEq = attr(mutmod, "sexEqual")
  if(is.null(sexEq))
    stop2("Mutation model attribute `sexEqual` is not set")

  if(sexEq && !identical(male, female))
    stop2("Mutation model attribute `sexEqual` is falsely set to TRUE")

  if(!sexEq) {
    validateMutationMatrix(female, alleles = alleles)

    if(!identical(attr(male, 'afreq'), attr(female, 'afreq')))
      stop2("Mutation model attribute `afreq` differs for males and females")
  }

  mutmod
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
    cat("* * * * * * * * * * *\n")
    cat("Male mutation matrix:\n")
    print(x$male, ...)
  }
  invisible(x)
}

#' @export
toString.mutationModel = function(x, ...) {
  if(attr(x, "sexEqual"))
    toString(x$female)
  else
    sprintf("%s (male); %s (female)", toString(x$male), toString(x$female))
}

isMutationModel = function(x) {
  class(x) ==  "mutationModel"
}

# Permute the allele order of a model or matrix
enforceAlleleOrder = function(m, alleles) {
  if(is.null(alleles))
    return(m)

  if(isMutationModel(m)) {
    m$male = enforceAlleleOrder(m$male, alleles)
    m$female = enforceAlleleOrder(m$female, alleles)
    return(m)
  }

  if(!isMutationMatrix(m))
    stop2("Expected a mutation matrix, not a: ", class(m))

  nms = dimnames(m)
  alleles = as.character(alleles)

  # If already correct order - return
  if(identical(alleles, nms[[1]]))
    return(m)

  if(!setequal(alleles, nms[[1]]) || !setequal(alleles, nms[[2]]))
    stop2("Alleles differ from names of mutation matrix")

  if(length(alleles) > length(nms[[1]]))
    stop2("Duplicated alleles indicated: ", alleles[duplicated[alleles]])

  # Permute
  new_m = m[alleles, alleles]

  # Keep all attributes of original, except `dimnames` and `afreq`
  attrs = attributes(m)
  attrs$dimnames = list(alleles, alleles)
  if(!is.null(attrs$afreq))
    attrs$afreq = attrs$afreq[alleles]

  attributes(new_m) = attrs

  new_m
}


#' @rdname mutationModel
#' @export
sexEqual = function(mutmod) {
  isTRUE(attr(mutmod, "sexEqual"))
}
