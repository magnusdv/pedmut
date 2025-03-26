#' Mutation matrix
#'
#' Construct mutation matrices for pedigree likelihood computations.
#'
#' Descriptions of the models:
#'
#' * `custom` : Allows any mutation matrix to be provided by the user, in the
#' `matrix` parameter.
#'
#' * `dawid` : A reversible model for integer-valued markers, proposed by Dawid
#' et al. (2002).
#'
#' * `equal` :  All mutations equally likely; probability \eqn{1-rate} of no
#' mutation.
#'
#' * `proportional` : Mutation probabilities are proportional to the target
#' allele frequencies.
#'
#' * `random` : This produces a matrix of random numbers, where each row is
#' normalised so that it sums to 1.
#'
#' * `onestep`: A mutation model for microsatellite markers, allowing mutations
#' only to the nearest neighbours in the allelic ladder. For example, '10' may
#' mutate to either '9' or '11', unless '10' is the lowest allele, in which case
#' '11' is the only option. This model is not applicable to loci with
#' non-integral microvariants.
#'
#' * `stepwise`: A common model in forensic genetics, allowing different
#' mutation rates between integer alleles (like '16') and non-integer
#' "microvariants" like '9.3'). Mutations also depend on the size of the
#' mutation if the parameter 'range' differs from 1.
#'
#' * `trivial` : The identity matrix; i.e. no mutations are possible.
#'
#' @param model A string: either "custom", "dawid", "equal", "proportional",
#'   "random", "stepwise" or "onestep".
#' @param matrix When `model` is "custom", this must be a square matrix with
#'   nonnegative real entries and row sums equal to 1.
#' @param alleles A character vector (or coercible to character) with allele
#'   labels. Required in all models, except "custom" if `matrix` has dimnames.
#' @param afreq A numeric vector of allele frequencies. Required in model
#'   "proportional".
#' @param rate A number between 0 and 1. Required in models "equal",
#'   "proportional", "stepwise" and "onestep".
#' @param seed A single number. Optional parameter in the "random" model, passed
#'   on to `set.seed()`.
#' @param rate2 A number between 0 and 1. The mutation rate between integer
#'   alleles and microvariants. Required in the "stepwise" model.
#' @param range A positive number. The relative probability of mutating n+1
#'   steps versus mutating n steps. Required  in the "stepwise" model.
#' @param mutmat An object of class `mutationMatrix`.
#'
#' @return An object of class `mutationMatrix`, essentially a square matrix with
#'   entries in `[0, 1]`, and whose colnames and rownames are the allele labels.
#'
#' @examples
#' mutationMatrix(alleles = 1:3, model = "equal", rate = 0.05)
#'
#' @importFrom stats runif
#' @export
mutationMatrix = function(model = c("custom", "dawid", "equal", "proportional",
                                    "random", "onestep", "stepwise", "trivial"),
                          matrix = NULL, alleles = NULL, afreq = NULL,
                          rate = NULL, seed = NULL, rate2 = NULL, range = NULL) {

  model = match.arg(model)
  alleles = if(!is.null(alleles)) as.character(alleles) else names(afreq)
  afreq = checkAfreq(afreq, alleles)

  # A few ad hocs for custom models
  if(model == "custom")
    alleles = alleles %||% colnames(matrix)
  else if(!is.null(matrix))
    stop2(sprintf("`matrix` cannot be used with the `%s` model", model))

  mutmat = switch(model,
    custom = .custom(matrix, alleles) |> validateMutationMatrix(),
    dawid = .dawid(alleles, afreq, rate, range),
    equal = .equal(alleles, rate),
    proportional = .proportional(alleles, afreq, rate),
    random = .random(alleles, seed),
    stepwise = .stepwise(alleles, rate, rate2, range),
    onestep = .onestep(alleles, rate),
    trivial = .trivial(alleles)
  )

  newMutationMatrix(mutmat, model = model, afreq = afreq, rate = rate,
                    rate2 = rate2, range = range, seed = seed)
}

newMutationMatrix = function(mutmat, model = "custom", afreq = NULL,
                            rate = NULL, rate2 = NULL, range = NULL,
                            lumpedAlleles = NULL, seed = NULL) {
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
            rate2 = rate2,
            range = range,
            stationary = stationary,
            reversible = reversible,
            lumpedAlleles = lumpedAlleles,
            seed = seed,
            class = "mutationMatrix")
}


#' @rdname mutationMatrix
#' @export
validateMutationMatrix = function(mutmat, alleles = NULL) {

  dm = dim(mutmat)
  if(is.null(dm) || dm[1] != dm[2] || !is.numeric(mutmat))
    stop2("Mutation matrix is not a square numeric matrix")

  if(anyNA(mutmat))
    stop2("Mutation matrix contains NA elements")

  nms = dimnames(mutmat)
  if(!identical(nms[[1]], nms[[2]]))
    stop2("Mutation matrix has unequal row and column names")

  als = nms[[1]]
  if(!is.null(alleles) && !identical(als, as.character(alleles)))
    stop2("Dimnames of mutation matrix not consistent with indicated alleles")

  afreq = attr(mutmat, "afreq")
  if(!is.null(afreq)) {
    if(!identical(names(afreq), als)) {
      print(afreq)
      stop2("Attribute `afreq` must be named with the allele labels: ", als)
    }
    if(round(sum(afreq), 3) != 1)
      stop2("Allele frequencies do not sum to 1 (after rounding to 3 decimal places): ",
            round(afreq, 3))
  }

  if(any(mutmat < 0))
    stop2("Negative entries found in mutation matrix: ", mutmat[mutmat < 0])
  if(any(mutmat > 1))
    stop2("Entries exceeding 1 found in mutation matrix: ", mutmat[mutmat > 1])

  rs = round(.rowSums(mutmat, dm[1], dm[1]), 3)
  if (any(rs != 1)) {
    rw = which(rs != 1)[1]
    stop2(sprintf("Row %d doesn't sum to 1 (after rounding to 3 decimals), but to: %s", rw, rs[rw]))
  }

  mutmat
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
    if(is.null(rate) && !is.null(afreq))
      rate = mutRate(x, afreq)
    rate2 = attrs$rate2
    range = attrs$range
    seed = attrs$seed

    cat("Model:", model, "\n")
    cat("Rate:", if(!is.null(rate)) rate else NA, "\n")
    if(!is.null(afreq))
      cat("Frequencies:", toString(afreq), "\n")
    if(model == "random")
      cat("Seed:", if(!is.null(seed)) seed else NA, "\n")
    if(model == "stepwise"){
      cat("Rate2:", rate2, "\n")
      cat("Range:", range,  "\n")
    }
  }

  if(includeProperties) {
    if(includeAttrs) cat("\n")

    if(!is.null(afreq)) {
      cat("Bounded:", if(isBounded(x, afreq)) "Yes" else "No", "\n")
      cat("Stationary:", if(isStationary(x, afreq)) "Yes" else "No", "\n")
      cat("Reversible:", if(isReversible(x, afreq)) "Yes" else "No", "\n")
    }
    cat("Lumpable:", if(alwaysLumpable(x)) "Always" else "Not always", "\n")
  }
}

#' @export
toString.mutationMatrix = function(x, ...) {
  attrs = attributes(x)
  mod = attrs$model
  param = switch(mod,
                 equal =, proportional =, onestep = paste("rate =", attrs$rate),
                 random = paste("seed =", attrs$seed %||% "NULL"),
                 stepwise = sprintf("rate = %g, rate2 = %g, range = %g", attrs$rate, attrs$rate2, attrs$range),
                 NULL)
  if(!is.null(param))
    mod = sprintf("%s, %s", mod, param)

  mod
}



# Specific model definitions ----------------------------------------------

.trivial = function(alleles) {
  mutmat = diag(length(alleles))
  dimnames(mutmat) = list(alleles, alleles)
  mutmat
}

.equal = function(alleles, rate) {
  checkNullArg(alleles, "equal")
  checkRate(rate, "equal")

  nall = length(alleles)
  m = rate/(nall - 1)
  mutmat = matrix(m, ncol = nall, nrow = nall, dimnames = list(alleles, alleles))
  diag(mutmat) = 1 - rate
  mutmat
}

.proportional = function(alleles, afreq, rate) {
  checkNullArg(alleles, "proportional")
  checkRate(rate, "proportional")
  afreq = checkAfreq(afreq, alleles, checkNULL = TRUE, model = "proportional")

  # Check if mutation rate is too large
  mxr = (1 - sum(afreq^2)) / (1 - min(afreq))
  if(rate > mxr)
    stop2("Model undefined: max `rate` for the given input is: ", mxr)

  nall = length(alleles)
  a = rate / sum(afreq * (1 - afreq))
  mutmat = (1 - a) * diag(nall) + a * rep(afreq, each = nall)

  dimnames(mutmat) = list(alleles, alleles)
  mutmat
}

.random = function(alleles, seed = NULL) {
  checkNullArg(alleles, "random")
  if(!is.null(seed))
    set.seed(seed)

  nall = length(alleles)
  m = runif(nall^2, min = 0, max = 1)
  mutmat = matrix(m, ncol = nall, nrow = nall, dimnames = list(alleles, alleles))
  mutmat / rowSums(mutmat)
}

.stepwise = function(alleles, rate, rate2, range) {
  checkNullArg(alleles, "stepwise")
  checkRate(rate, "stepwise")
  checkRate2(rate2, "stepwise")
  checkRange(range, "stepwise")
  if(rate + rate2 > 1)
    stop2("The total mutation rate `rate + rate2` must be in [0,1]: ", rate + rate2)

  nall = length(alleles)
  if(nall == 0)
    return(NULL)

  alsNum = checkNumericAlleles(alleles, "stepwise")

  # Bug fix: round these!
  microgroup = round((alsNum - round(alsNum))*10)

  # Initialise matrix
  mutmat = matrix(0, ncol = nall, nrow = nall, dimnames = list(alleles, alleles))

  for (i in 1:nall) {
    microcompats = (microgroup == microgroup[i])
    for (j in 1:nall) {
      if (i == j) {
        if (all(microcompats)) mutmat[i,j] = 1 - rate
        else if (sum(microcompats) == 1) mutmat[i,j] = 1 - rate2
        else mutmat[i,j] = 1 - rate - rate2
      } else if (microcompats[j])
        mutmat[i,j] = range^abs(alsNum[i] - alsNum[j])
      else
        mutmat[i,j] = rate2/(nall - sum(microcompats))
    }
    microcompats[i] = FALSE
    if (any(microcompats))
      mutmat[i, microcompats] = mutmat[i, microcompats]/sum(mutmat[i, microcompats]) * rate
  }

  mutmat
}

.onestep = function(alleles, rate) {
  alsNum = checkIntegerAlleles(alleles, "onestep")
  checkRate(rate, "onestep")

  # Initialise matrix
  nall = length(alleles)
  mutmat = matrix(0, ncol = nall, nrow = nall, dimnames = list(alleles, alleles))

  # Loop over rows
  for(i in seq_along(alsNum)) {
    neigh = abs(alsNum - alsNum[i]) == 1 # onestep neighbours
    mutmat[i, i] = if(sum(neigh) > 0) 1 - rate else 1
    mutmat[i, neigh] = rate/sum(neigh) # either rate or rate/2
  }

  mutmat
}

.dawid = function(alleles, afreq, rate, range) {
  alsNum = checkIntegerAlleles(alleles, "dawid")
  afreq = checkAfreq(afreq, alleles, checkNULL = TRUE, model = "dawid")

  checkRate(rate, "dawid")
  checkRange(range, "dawid")

  n = length(alleles)
  a = (1 - range^n)/(1 - range)

  stepsize = outer(1:n, 1:n, function(i,j) abs(i-j))
  R = rate/afreq * (1 - range) / (2*range*(n - a)) * range^stepsize
  diag(R) = 0

  # Diagonal values
  dg = 1 - rowSums(R)

  # If undefined (negative diagonal values), return info on max rate
  if(any(dg < 0)) {
    mxr = 1/max(rowSums(R/rate)) # See formulas in Egeland & Vigeland (2025)
    stop2("Model undefined; max `rate` for the given input is: ", mxr)
  }

  diag(R) = dg
  dimnames(R) = list(alleles, alleles)
  R
}

.custom = function(matrix, alleles) {
  checkNullArg(matrix, "custom")

  if(!is.matrix(matrix))
    stop2("Custom matrix must be a matrix, not a ", class(matrix))
  if(typeof(matrix) != "double")
    stop2("Custom matrix must be numeric, not ", typeof(matrix))
  if(nrow(matrix) != ncol(matrix))
    stop2("Custom matrix must be square, not ", nrow(matrix), "*", ncol(matrix))

  nullAls = is.null(alleles)
  if(!nullAls) alleles = as.character(alleles)

  # If both dimnames and `alleles` are given: Check match and use to permute matrix
  if(!is.null(dmn <- dimnames(matrix)) && !nullAls) {
    chck = setequal(alleles, dmn[[1L]]) && setequal(alleles, dmn[[2L]])
    if(!chck)
      stop2("Custom matrix names don't match the `alleles` argument")

    matrix = matrix[alleles, alleles]
  }
  else if(!nullAls) {
    if(length(alleles) != dim(matrix)[1L])
      stop2("Length of `alleles` must equal the dimension of `matrix`")

    dimnames(matrix) = list(alleles, alleles)
  }
  else if(is.null(dmn) && nullAls)
    stop2("When custom matrix lacks names, the argument `alleles` cannot be NULL")

  matrix
}



# Utilities for checking parameters-------------------------------------------

checkNullArg = function(arg, model = NULL) {
  if(is.null(arg)) {
    argname = deparse(substitute(arg))
    if(is.null(model))
      msg = sprintf("`%s` cannot be NULL with this model", argname)
    else
      msg = sprintf("`%s` cannot be NULL with the `%s` model", argname, model)

    stop2(msg)
  }
}

checkIntegerAlleles = function(alleles, model) {
  checkNullArg(alleles, model)
  alsNum = suppressWarnings(as.numeric(alleles))
  if(any(is.na(alsNum)) || any(round(alsNum) != alsNum)) {
    msg = sprintf("The `%s` model requires all alleles to be integers", model)
    stop2(msg)
  }
  alsNum
}

checkNumericAlleles = function(alleles, model) {
  alsNum = suppressWarnings(round(as.numeric(alleles), 2))
  if(any(is.na(alsNum)))
    stop2(sprintf("The `%s` model requires all alleles to be numeric", model))
  oneDec = round(alsNum, 1) == alsNum
  if(!all(oneDec))
    stop2("Microvariants must be named as a decimal number with one decimal: ", alleles[!oneDec])
  alsNum
}

checkAfreq = function(afreq, alleles = NULL, len = NULL, checkNULL = FALSE, model = NULL) {
  if(checkNULL)
    checkNullArg(afreq, model)

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

checkRate = function(rate, model) {
  checkNullArg(rate, model)
  if(!isNumber(rate, minimum = 0, maximum = 1))
    stop2("`rate` must be a number in the interval `[0,1]`: ", rate)
}

checkRate2 = function(rate2, model) {
  checkNullArg(rate2, model)
  if(!isNumber(rate2, minimum = 0, maximum = 1))
    stop2("`rate2` must be a number in the interval `[0,1]`: ", rate2)
}

checkRange = function(range, model) {
  checkNullArg(range, model)
  if(!isNumber(range, minimum = 0) || range == 0)
    stop2("`range` must be a positive number: ", range)
}
