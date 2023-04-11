stop2 = function(...) {
  a = lapply(list(...), toString)
  a = append(a, list(call. = FALSE))
  do.call(stop, a)
}

# Test that input is a single positive (or similar) integer.
isCount = function(x, minimum = 1) {
  isTRUE(length(x) == 1 &&
           (is.integer(x) || (is.numeric(x) && x == as.integer(x))) &&
           x >= minimum)
}

# Test that input is a single number, with optional range constraints
isNumber = function(x, minimum = NA, maximum = NA) {
  isTRUE(length(x) == 1 &&
           is.numeric(x) &&
           (is.na(minimum) || x >= minimum) &&
           (is.na(maximum) || x <= maximum))
}


`%||%` = function(x, y) {
  if(is.null(x)) y else x
}

# A safer version of base::sample
safe_sample <- function(x, ...) x[sample.int(length(x), ...)]

# Order numerically if appropriate otherwise lex
smartOrder = function(x) {
  if (!is.numeric(x) && !anyNA(suppressWarnings(as.numeric(x))))
    ord = order(as.numeric(x))
  else
    ord = order(x)
  ord
}

# Fast setdiff
.mysetdiff = function(x, y) unique.default(x[match(x, y, 0L) == 0L])

# Fast intersection. NB: assumes no duplicates!
.myintersect = function(x, y) y[match(x, y, 0L)]

.equalNums = function(v, tol = sqrt(.Machine$double.eps)) {
  all(v == v[1]) || abs(max(v) - min(v)) < tol
}

checkAfreq = function(afreq, alleles = NULL, len = NULL) {
  if(is.null(afreq))
    return(afreq)

  if(!is.numeric(afreq))
    stop2("Expected frequency vector to be a numeric, not ", class(afreq))

  if(!is.null(len) && length(afreq) != len)
    stop2(sprintf("Expected frequency vector to have length %d, not %d", len, length(afreq)))

  if(round(sum(afreq), 3) != 1)
    stop2("Allele frequencies do not sum to 1: ", round(sum(afreq),3))

  if(!is.null(alleles)) {

    if(length(afreq) != length(alleles))
      stop2(sprintf("Frequency vector does not match the number of alleles (%d)", length(alleles)))

    if(!is.null(nms <- names(afreq))) {
      if(!setequal(nms, alleles))
        stop2("Names of frequency vectors do not match alleles: ", setdiff(nms, alleles))

      # Sort
      afreq = afreq[alleles]
    }
    else
      names(afreq) = alleles
  }

  afreq
}

#' @export
as.matrix.mutationMatrix = function(x, ...) {
  attributes(x) = list(dim = dim(x), dimnames = dimnames(x))
  x
}

#' @export
as.matrix.mutationModel = function(x, ...) {
  if(!sexEqual(x))
    stop2("Male and female matrices are not equal")
  m = x$female
  attributes(m) = list(dim = dim(m), dimnames = dimnames(m))
  m
}

#' Test for mutation matrix/model
#'
#' @param x Any object.
#' @returns TRUE or FALSE
#' @examples
#'
#' mat = mutationMatrix("equal", alleles = 1:2, rate = 0.1)
#' isMutationMatrix(mat)
#'
#' isMutationModel(mat) # FALSE (not a complete model)
#'
#' mod = mutationModel(mat)
#' isMutationModel(mod)
#'
#' @export
isMutationModel = function(x) {
  inherits(x, "mutationModel")
}

#' @rdname isMutationModel
#' @export
isMutationMatrix = function(x) {
  inherits(x, "mutationMatrix")
}
