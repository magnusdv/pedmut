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
#' @param matrix A matrix, or a list of two (named "female" and "male")
#' @param rate A numeric mutation rate, or a list of two (named "female" and
#'   "male")
#' @param seed An integer, or a list of two (named "female" and "male")
#' @param ... Further arguments to [mutationMatrix()], which are reused for both
#'   models. Typically `alleles` and/or `afreq`.
#' @param mutmod A `mutationModel` object
#' @param alleles A character vector with allele labels. (The validation method
#'   uses this to check that the matrices have appropriate dimnames.)
#'
#' @return An object of class `mutationModel`.
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
mutationModel = function(model, matrix = NULL, rate = NULL, seed = NULL, ...) {

  if(isMutationModel(model)) {
    return(validateMutationModel(model))
  }

  if(is.list(model) &&
          length(model) == 2 &&
          setequal(names(model), c("female", "male")) &&
          all(sapply(model, isMutationMatrix))) {
    mod = model
  }
  else if(isMutationMatrix(model)) {
    mod = list(female = model, male = model)
  }
  else if((is.character(model) && length(model) == 1) || is.list(model) && length(model) == 2) {
    # Build models from parameters
    model = validateSingleInput(model, "character")
    matrix = validateMatrixInput(matrix)
    rate = validateSingleInput(rate, "numeric")
    seed = validateSingleInput(seed, "integer")

    female = mutationMatrix(model = model$female, matrix = matrix$female,
                            rate = rate$female, seed = seed$female, ...)

    male = mutationMatrix(model = model$male, matrix = matrix$male,
                            rate = rate$male, seed = seed$male, ...)

    mod = list(female = female, male = male)
  }
  else
    stop2("`model` must be either:\n",
          " * an object of class `mutationModel`\n",
          " * an object of class `mutationMatrix`` (applied to both genders)`\n",
          " * a list of two `mutationMatrix` objects, named 'female' and 'male'\n",
          " * a character string (see ?mutationModel for valid options)\n",
          " * a list of two characters, named 'female' and 'male'")

  mutmod = structure(mod,
                     sexEqual = identical(mod$female, mod$male),
                     class = "mutationModel")
  validateMutationModel(mutmod)
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
  stopifnot(is.list(mutmod),
            length(mutmod) == 2,
            setequal(names(mutmod), c("male", "female")),
            identical(attr(mutmod, "sexEqual"), identical(mutmod$male, mutmod$female)),
            inherits(mutmod, "mutationModel"))

  validateMutationMatrix(mutmod$male, alleles = alleles)
  validateMutationMatrix(mutmod$female, alleles = alleles)

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

    cat("\nMale mutation matrix:\n")
    print(x$male, ...)
  }
  invisible(x)
}

isMutationModel = function(x) {
  class(x) ==  "mutationModel"
}
