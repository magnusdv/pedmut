#' Reversible stepwise mutation model
#'
#' A reversible stepwise mutation model is created following the approach of
#' Dawid et al. (2002).
#'
#' @param alleles A vector of integer integers.
#' @param afreq A numeric vector of allele frequencies.
#' @param rate A numeric mutation rate.
#' @param range A positive number.
#'
#' @return A reversible stepwise mutation model with expected mutation rate
#'   equal input rate.
#'
#' @details The function [boundsGamma()] checks if it is possible to construct
#'   the mutation matrix.
#'
#' @author Thore Egeland.
#'
#' @export
#'
#' @examples
#' stepwiseReversible(alleles = 1:3, afreq = c(0.2, 0.3,  0.5), rate = 0.001, range = 0.1)
#'
#' # Model not well defined:
#' \dontrun{
#' stepwiseReversible(alleles = 1:3, afreq = c(0.2, 0.3,  0.5), rate = 0.7, range = 0.1)
#' }
stepwiseReversible = function(alleles, afreq, rate, range){
  if(!is.integer(alleles) && !(is.numeric(alleles) && all(alleles == as.integer(alleles))))
    stop2("Non-integer alleles detected")
  if(!is.numeric(range) || (range <= 0 || range >= 1))
    stop2("`range` must be in the interval (0,1): ", range)

  # remaining checking will be taken care of by `mutationModel` below
  n = length(afreq)

  if(boundsGamma(alleles, afreq,  range)[[1]] < rate)
    stop2("Model not well defined")

  R = matrix(ncol = n, nrow = n, 0)
    for (i in 1:n){
      for(j in (1:n)[-i]) {
          a = (1 - range^n)/(1 - range)
          R[i,j] = rate * (1 - range) * range^{abs(i-j)}/
                   (2*range*(n - a))*(1/afreq[i])
        }
      R[i,i] = 1 - sum(R[i,-i])
    }

  dimnames(R) = list(alleles, alleles)
  mutationModel(matrix = R, model = "custom",  afreq = afreq, alleles = alleles)
}

