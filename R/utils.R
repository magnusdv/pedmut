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
