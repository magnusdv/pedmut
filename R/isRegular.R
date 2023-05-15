#' Checks if mutation model is regular
#' 
#' Checks regularity, i.e, if m_ij <= p_j.
#' 
#'  
#' @param M a mutation model.
#' @param afreq A numeric vector of allele frequencies.
#' 
#' @return A logical.
#' 
#' 
#' @author Thore Egeland.
#' 
#' @export
#' 
#' @examples
#' 
#' p = c(0.4, 0.6)
#' R = stepwiseReversible(alleles = 1:2, afreq = p, rate = 0.6, range = 0.1)
#' isRegular(R, afreq = p)


isRegular = function(M, afreq = NULL){
  if(is.null(afreq))
    afreq = attr(M, "afreq")
  n = length(afreq)
  M = as.matrix(M)
  lines = rep(FALSE, n)
  for (i in 1:n)
    lines[i] = all(M[-i, i] <= afreq[i])
  all(lines)
}

