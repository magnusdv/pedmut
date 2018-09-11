#' Mutation models
#'
#' Constructor for the class `mutationModel`. An object of this class is
#' essentially a list of two mutation matrices, named "female" and "male".
#'
#' @param model Either:
#'
#'   * a single `mutationMatrix`
#'
#'   * a list of two `mutationMatrix` objects, named "male" and "female"
#'
#'   * a single character string (see [mutationMatrix()] for valid options)
#'
#'   * a list of two strings, named (named "male" and "female")
#'
#' @param matrix Either a single matrix, or a list of two (named "male" and
#'   "female")
#' @param rate Either a single number or a list of two (named "male" and
#'   "female")
#' @param seed Either a single integer or a list of two (named "male" and
#'   "female")
#' @param ... Further arguments to [mutationMatrix()], which are reused for both
#'   models. Typically `alleles` and/or `afreq`.
#' @param mutmod A `mutationModel` object
#' @param alleles A character vector with allele labels. (The validation method
#'   uses this to check that the matrices have appropriate dimnames.)
#'
#' @return An object of class `mutationModel`.
#'
#' @examples
#' # Same model for both genders
#' M1 = mutationModel("eq", alleles = 1:2, rate = 0.1)
#' M1
#' stopifnot(identical(M1$male, M1$female))
#'
#' # Different mutation rates
#' M2 = mutationModel("eq", alleles = 1:2, rate = list(male = 0.1, female = 0.01))
#' M2
#' stopifnot(identical(M2$male, M1$male))
#'
#' @export
mutationModel = function(model, matrix = NULL, rate = NULL, seed = NULL, ...) {
  # If "model" is a mutationMatrix:
  # return mutationModel with this for both males and females
  if(isMutationMatrix(model)) {
    mod = structure(list(female = model, male = model),
                    sexEqual = TRUE, class = "mutationModel")
    return(validateMutationModel(mod))
  }

  # If "model" is list of two mutationMatrix objects:
  # validate and return
  if(is.list(model) && all(sapply(model, inherits, "mutationMatrix"))) {
    stopifnot(length(model) == 2,
              setequal(names(model), c("female", "male")))
    mod = structure(model, sexEqual = identical(model$male, model$female),
                    class = "mutationModel")
    return(validateMutationModel(mod))
  }

  # Otherwise: build models from parameters
  model = validateSingleInput(model, "character")
  matrix = validateMatrixInput(matrix)
  rate = validateSingleInput(rate, "numeric")
  seed = validateSingleInput(seed, "integer")

  female = mutationMatrix(model = model$female, matrix = matrix$female,
                          rate = rate$female, seed = seed$female, ...)

  male = mutationMatrix(model = model$male, matrix = matrix$male,
                          rate = rate$male, seed = seed$male, ...)

  mod = structure(list(female = female, male = male),
                  sexEqual = identical(female, male),
                  class = "mutationModel")
  validateMutationModel(mod)
}

validateSingleInput = function(x, mode) {
  if(is.null(x))
    return(x)

  if(is.vector(x, mode = mode))
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
