#' Mutation matrix
#'
#' Construct mutation matrices for pedigree likelihood computations.
#'
#' Descriptions of the models:
#'
#' * `custom` : Allows any mutation matrix to be provided by the user, in the
#' `matrix` parameter
#'
#' * `equal` :  All mutations equally likely; probability \eqn{1-rate} of no
#' mutation
#'
#' * `proportional` : Mutation probabilities are proportional to the target
#' allele frequencies
#'
#' * `random` : This produces a matrix of random numbers, where each row is
#' normalised so that it sums to 1
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
#' @param model A string: either "custom", "equal", "proportional", "random",
#'   "stepwise" or "onestep"
#' @param matrix When `model` is "custom", this must be a square matrix with
#'   nonnegative real entries and row sums equal to 1
#' @param alleles A character vector (or coercible to character) with allele
#'   labels. Required in all models, except "custom" if `matrix` has dimnames
#' @param afreq A numeric vector of allele frequencies. Required in model
#'   "proportional"
#' @param rate A number between 0 and 1. Required in models "equal",
#'   "proportional", "stepwise" and "onestep"
#' @param seed A single number. Optional parameter in the "random" model, passed
#'   on to `set.seed()`
#' @param rate2 A number between 0 and 1. The mutation rate between integer
#'   alleles and microvariants. Required in the "stepwise" model
#' @param range A positive number. The relative probability of mutating n+1
#'   steps versus mutating n steps. Required  in the "stepwise" model
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
mutationMatrix = function(model = c("custom", "equal", "proportional",
                                    "random", "onestep", "stepwise", "trivial"),
                          matrix = NULL, alleles = NULL, afreq = NULL,
                          rate = NULL, seed = NULL, rate2 = NULL, range = NULL) {
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

      # If alleles are given, make sure they match - and use to permute matrix
      if(!is.null(alleles)) {
        if(!setequal(alleles, mnames))
          stop2("Custom matrix names don't match the `alleles` argument")

        alleles = as.character(alleles)
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
  if(model %in% c("equal", "proportional", "stepwise", "onestep")) {
    if(is.null(rate))
      stop2("`rate` cannot be NULL with this model")
    if(!isNumber(rate, minimum = 0))
      stop2("`rate` must be a nonnegative number: ", rate)
  }
  if(model %in% c("stepwise")) {
    if(is.null(rate2))
      stop2("`rate2` cannot be NULL with the `stepwise` model")
    if(!isNumber(rate2, minimum = 0))
      stop2("`rate2` must be a nonnegative number: ", rate2)
    if(rate + rate2 > 1)
      stop2("The total mutation rate `rate + rate2` must be in [0,1]: ", rate + rate2)
    if(is.null(range))
      stop2("`range` cannot be NULL with the `stepwise` model")
    if(!(isNumber(range, minimum = 0) && range > 0))
      stop2("`range` must be a positive number: ", range)
  }

  ## Compute matrix according to model
  nall = length(alleles)
  if(nall == 0)
    return(NULL)
  if(nall == 1)
    model = "trivial"

  mutmat = matrix(0, ncol = nall, nrow = nall, dimnames = list(alleles, alleles))

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
  else if(model == "trivial") {
    mutmat[] = diag(nall)
  }
  else if (model == "stepwise") {
    alsNum = suppressWarnings(as.numeric(alleles))
    if (any(is.na(alsNum)))
      stop2("The `stepwise` mutation model requires all alleles to have numerical names")
    if (any(round(alsNum, 1) != alsNum))
      stop2("Microvariants must be named as a decimal number with one decimal")
    microgroup = (alsNum - round(alsNum))*10
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
  }
  else if (model == "onestep") {
    alsNum = suppressWarnings(as.numeric(alleles))
    if (any(is.na(alsNum)) || any(round(alsNum) != alsNum))
      stop2("The `onestep` mutation model requires all alleles to be integers")

    # Loop over rows
    for(i in seq_along(alsNum)) {
      # Columns with neighbour alleles
      nei = abs(alsNum - alsNum[i]) == 1

      mutmat[i, i] = if(sum(nei) > 0) 1 - rate else 1

      # Insert either rate or rate/2
      mutmat[i, nei] = rate/sum(nei)
    }
  }

  if(!is.null(afreq))
    names(afreq) = alleles

  newMutationMatrix(mutmat, model=model, afreq=afreq, rate=rate,
                    rate2 = rate2, range = range, seed=seed)
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

  rs = .rowSums(mutmat, dm[1], dm[1])
  if (any(round(rs, 3) != 1)) {
    print(rowSums(mutmat))
    stop2("Rows which do not sum to 1 (after rounding to 3 decimal places): ",
          which(round(rs, 3) != 1))
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

isMutationMatrix = function(x) {
  inherits(x, "mutationMatrix")
}
