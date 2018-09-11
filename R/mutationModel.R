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
mutationModel = function(female, male = NULL, ...) {

  if(class(female) != "mutationMatrix")
    female = mutationMatrix(model = female, ...)

  if(is.null(male))
    male = female

  mod = structure(list(female = female, male = male),
                  sexEqual = identical(female, male),
                  class = "mutationModel")
  validateMutationModel(mod)
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
