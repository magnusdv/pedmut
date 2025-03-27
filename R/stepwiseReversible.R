#' Dawid's reversible stepwise model
#'
#' A reversible stepwise mutation model is created following the approach of
#' Dawid et al. (2002). NB: This function is now implemented in the
#' `mutationMatrix()` function with the argument `model = "dawid"`.
#'
#' For the stepwise reversible model, the mutation rate \eqn{r_{i,j},\,  i\neq
#' j} is proportional to the overall mutation rate \eqn{\lambda} for given
#' values of the range, the allele frequency \eqn{p_i} and n, the number of
#' alleles. Hence, one can determine bounds UW and UB so that the model is well
#' defined if \eqn{\lambda \leq UW} and bounded, i.e., \eqn{r_{i,j} \leq p_j,\,
#' i\neq j}, if \eqn{\lambda \leq UB}, The bounds UW and UB are computed.
#'
#' @param alleles A vector of integer integers.
#' @param afreq A numeric vector of allele frequencies.
#' @param rate A numeric mutation rate.
#' @param range A positive number.
#' @param maxRateOnly A logical, by default FALSE. See Value.
#'
#' @return A reversible stepwise mutation model with overall mutation rate equal
#'   to `rate`.
#'
#'   If `maxRateOnly` is TRUE, the function returns a vector of two numbers
#'   named `UW` and `UB`. The first of these is the maximum overall mutation
#'   rate for a well-defined stepwise reversible mutation matrix with the given
#'   input. The latter (UB) is the maximum rate under the additional restraint
#'   that the model is bounded by `afreq`.
#'
#' @author Thore Egeland.
#'
#' @export
#'
#' @examples
#' stepwiseReversible(alleles = 1:3,
#'                    afreq = c(0.2, 0.3,  0.5),
#'                    rate = 0.001,
#'                    range = 0.1)
#'
#' stepwiseReversible(alleles = 1:3,
#'                    afreq = c(0.2, 0.3,  0.5),
#'                    range = 0.1,
#'                    maxRateOnly = TRUE)
#'
#' # Model not well defined:
#' \dontrun{
#' stepwiseReversible(alleles = 1:3,
#'                    afreq = c(0.2, 0.3,  0.5),
#'                    rate = 0.7,
#'                    range = 0.1)
#' }
stepwiseReversible = function(alleles, afreq, rate, range, maxRateOnly = FALSE) {
  if(!is.integer(alleles) && !(is.numeric(alleles) && all(alleles == as.integer(alleles))))
    stop2("Non-integer alleles detected")
  if(!is.numeric(range) || (range <= 0 || range >= 1))
    stop2("`range` must be in the interval (0,1): ", range)

  mxr = maxRate(alleles, afreq,  range)

  if(maxRateOnly)
    return(mxr)

  if(mxr[["UW"]] < rate)
    stop2("Model not well defined; max `rate` for the given input is: ", mxr[["UW"]])

  # remaining checking will be taken care of by `mutationModel` below
  n = length(afreq)
  a = (1 - range^n)/(1 - range)

  R = matrix(ncol = n, nrow = n, 0)
    for (i in 1:n){
      for(j in (1:n)[-i]) {

          R[i,j] = rate * (1 - range) * range^{abs(i-j)}/
                   (2*range*(n - a))*(1/afreq[i])
        }
      R[i,i] = 1 - sum(R[i,-i])
    }

  alleles = as.character(alleles)
  dimnames(R) = list(alleles, alleles)
  mutationMatrix(matrix = R, model = "custom", afreq = afreq, alleles = alleles)
}



#' Upper limits for overall mutation rate for the stepwise reversible model.
#'
#' @param alleles A character vector with allele labels.
#' @param afreq A numeric vector of allele frequencies.
#' @param range A positive number.
#'
#' @return A vector of two numbers named `UW` and `UB`. The first of these is
#'   the maximum overall mutation rate for a well-defined stepwise reversible
#'   mutation matrix with the given input. The latter (UB) is the upper limit of
#'   the overall mutation rate under the additional restraint that the model is
#'   bounded by `afreq`.
#'
#' @author Thore Egeland.
#'
maxRate = function(alleles, afreq,  range){
  n = length(afreq)

  R1 = matrix(ncol = n, nrow = n, 0)

  for (i in 1:n){
    for(j in setdiff(1:n, i)){
      a = (1 - range^n)/(1 - range)
      R1[i,j] = 1 * (1 - range) * range^{abs(i-j)}/
        (2*range*(n - a))*(1/afreq[i])
    }
  }

  # Essential that diag(R1) = 0
  linesums = apply(R1, 1, sum)
  boundDefined = 1/max(linesums)
  maks = apply(R1, 2, max)
  c(UW = boundDefined,
    UB = min(afreq/maks))
}

# Faster version of maxRate
maxRate2 = function(afreq, range){
  n = length(afreq)
  stepsize = outer(1:n, 1:n, function(i,j) abs(i-j))

  a = (1 - range^n) / (1 - range)
  R = (1 - range) / (2*range*(n - a)) * range^stepsize / afreq
  diag(R) = 0

  maxDefined = 1/max(rowSums(R))
  maxBounded = min(afreq / apply(R, 2, max))
  c(UW = maxDefined, UB = maxBounded)
}
