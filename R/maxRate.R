#' Upper limits for overall mutation rate for the stepwise reversible model.
#'
#' @param alleles A character vector with allele labels.
#' @param afreq A numeric vector of allele frequencies.
#' @param range A positive number.
#'
#' @return Bounds on overall mutation rate for mutation matrix to be (i) well defined and (ii) regular.
#'
#' @details For the stepwise reversible model, the mutation rate \eqn{r_{i,j},\,  i\neq j}
#' is proportional to the overall mutation rate \eqn{\lambda} for given values of the range, the allele
#' frequency \eqn{p_i} and n, the number of alleles. Hence, we can determine bounds UW
#' and UB so that the model is well defined if \eqn{\lambda \leq UW} and bounded,
#' i.e., \eqn{r_{i,j} \leq p_j,\, i\neq j}, if \eqn{\lambda \leq UB}, The bounds UW and UB are computed.
#'
#' @author Thore Egeland.
#'
#' @export
#'
#' @examples
#' alleles = 1:3; afreq = c(0.2, 0.3,  0.5);  range = 0.1
#' bounds = maxRate(alleles, afreq ,range)
#' R1 = stepwiseReversible(alleles, afreq, rate = bounds[1], range)
#' isBounded(R1, afreq)
#' R2 = stepwiseReversible(alleles, afreq, rate = bounds[2], range)
#' isBounded(R2, afreq)
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
