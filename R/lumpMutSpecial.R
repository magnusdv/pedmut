#' Special lumping of mutation models
#'
#' This function implements methods for special, or pedigree-aware, allele
#' lumping. This is typically attempted if the model is not generally lumpable
#' as determined by [alwaysLumpable]. Note that the resulting lumped model is
#' tailor-made for a specific likelihood calculation, and may violate the
#' properties of a well-defined mutation model.
#'
#' The lumping procedure depend on the location of untyped individuals in the
#' pedigree, summarised by the so-called U-signature:
#'
#' * F-depth: The length of the longest chain of untyped, starting with a founder
#' * F-width: The maximum number of children of an untyped founder
#' * N-depth: The length of the longest chain of untyped, starting with a nonfounder
#' * N-width: The maximum number of children of an untyped nonfounder
#'
#' @param mut A square mutation matrix; typically a [mutationMatrix()] or
#'   [mutationModel()].
#' @param lump A vector containing the alleles to be lumped together.
#' @param uSign The U-signature of the pedigree for which lumping is attempted.
#'   See Details.
#' @param afreq A vector with allele frequencies, of the same length as the size
#'   of `mut`. Extracted from the model if not given.
#' @param verbose A logical.
#'
#' @returns A reduced mutation model, if lumping was possible, otherwise the
#'   original model is returned unchanged.
#'
#' @seealso [lumpedModel()].
#'
#' @examples
#'
#' af = rep(0.05, 20)
#' names(af) = 1:20
#' m = mutationMatrix("random", afreq = af, rate = 0.1, seed = 1)
#'
#' # Degree 1 lumping
#' mL = lumpMutSpecial(m, lump = 3:20, uSign = c(1,1,0,0))
#' mL
#'
#' # Check
#' afL = attr(mL, "afreq")
#' stopifnot(sum(af * m[, 1]) == sum(afL * mL[, 1]))
#'
#' # Degree 2
#' mL2 = lumpMutSpecial(m, lump = 3:20, uSign = c(1,2,0,0))
#' mL2
#' afL2 = attr(mL2, "afreq")
#'
#' stopifnot(all.equal(af %*% m[, 1]^2, afL2 %*% mL2[, 1]^2),
#'           all.equal(af %*% m[, 2]^2, afL2 %*% mL2[, 2]^2),
#'           all.equal(af %*% ( m[, 1]*m[, 2]),
#'                     afL2 %*% (mL2[, 1]*mL2[, 2])))
#'
#' @importFrom lpSolve lp
#' @export
lumpMutSpecial = function(mut, lump, uSign, afreq = NULL, verbose = TRUE) {

  isMod = isMutationModel(mut)
  sexeq = isMod && sexEqual(mut)

  # Extract frequencies (a bit ad hoc)
  afreq = afreq %||% attr(mut, "afreq") %||% attr(mut$female, "afreq") %||%
    stop2("`afreq` must be provided")

  als = if(isMod) colnames(mut$male) else colnames(mut)
  lump = prepLump(lump, alleles = als)[[1]]

  if(verbose) {
    nms = names(uSign)
    us = if(is.null(nms)) toString(uSign) else paste(nms, uSign, sep = "=", collapse = ", ")
    cat("U-signature:", us, "\n")
  }

  # Founder chains > 1, and nonfounder-chains not implemented yet
  if(uSign[1] > 1 || uSign[3] + uSign[4] > 0) {
    if(verbose) cat("Lumping not implemented for this signature - returning unchanged\n")
    return(mut)
  }

  # By now: uSign = c(1,d,0,0). All untyped are founders, with max d children
  d = uSign[2]

  if(d == 1) {
    if(verbose) cat("Lumping method: Direct\n")

    .directLump = function(mutmat, afreq, lump) {
      keep = .mysetdiff(names(afreq), lump)
      newAls = c(keep, "lump")
      newAfreq = c(afreq[keep], lump = sum(afreq[lump]))

      N = length(newAls)
      wei = afreq[lump]/sum(afreq[lump])

      newM = matrix(0, ncol = N, nrow = N, dimnames = list(newAls, newAls))
      newM[keep, keep] = mutmat[keep, keep]
      newM[N, -N] = as.numeric(wei %*% mutmat[lump, keep, drop = FALSE])
      newM[, N] = 1 - .rowSums(newM, N, N)
      newMutationMatrix(newM, afreq = newAfreq, lumpedAlleles = lump,
                   model = paste0(attr(mutmat, 'model'), ", lumped"))
    }

    if(!isMod) res = .directLump(mut, afreq, lump)
    else res = mapFullModel(mut, .directLump, afreq = afreq, lump = lump)
    return(res)
  }
  else {
    if(verbose) cat("Lumping method: LP\n")

    keep = .mysetdiff(als, lump)
    nk = length(keep)
    p = afreq[lump]

    # Main computation: Select subset of rows to keep
    A = (if(isMod) mut$female else mut)[lump, keep, drop = FALSE]
    B = if(isMod) mut$male[lump, keep, drop = FALSE] else NULL
    red = reduceLP(p, A = A, B = B, d = d, verbose = FALSE)

    # Was the reduction successful?
    if(red$status != 0) {
      if(verbose) cat("Lumping unsuccessful: No solution found\n")
      return(mut)
    }

    np = length(red$p)
    if(np == length(lump)) {
      if(verbose) cat("Lumping unsuccessful: No reduction\n")
      return(mut)
    }

    newAls = c(keep, paste0(".lmp", 1:np))
    newAfreq = c(afreq[keep], red$p)
    names(newAfreq) = newAls
    N = length(newAls)

    # Create reduced matrix; R is "lumped rectangle"
    .createLPmat = function(mutmat, R) {
      newM = matrix(0, ncol = N, nrow = N, dimnames = list(newAls, newAls))
      newM[1:nk, 1:nk] = mutmat[keep, keep]
      newM[nk + 1:np, 1:nk] = R
      newM[ , N] = 1 - rowSums(newM)
      newMutationMatrix(newM, afreq = newAfreq, lumpedAlleles = lump,
                        model = paste0(attr(mutmat, 'model'), ", lumped"))
    }

    # If single matrix, just returned the reduced matrix
    if(!isMod)
      return(.createLPmat(mut, red$A))

    # If full model, return new model
    matF = .createLPmat(mut$female, red$A)
    matM = .createLPmat(mut$male, red$B)
    lumpable = alwaysLumpable(matF) && (sexeq || alwaysLumpable(matM))
    res = structure(list(female = matF, male = matM),
              afreq = newAfreq, sexEqual = sexeq,
              alwaysLumpable = lumpable, class = "mutationModel")
    return(res)
  }
}


# Helper functions for special lumping ------------------------------------

reduceLP = function(p, A, B = NULL, d, verbose = TRUE) {
  n = length(p); scale = sum(p); p1 = p/scale

  # Monomial exponents: All tuples summing to at most d
  comb = as.matrix(expand.grid(rep(list(0:d), ncol(A))))
  comb = comb[rowSums(comb) %in% 1:d, , drop=FALSE]

  # Find monomials (of columns of A) of degree up to d
  monomsFix0 = function(Y) {
    if(all(Y > 0))
      return(exp(comb %*% t(log(Y))))

    # Temporarily replace 0 with 1 in log(Y) to avoid -Inf
    logY = log(replace(Y, Y == 0, 1))
    zeros = (comb > 0) %*% t.default(Y == 0) > 0
    monoms = exp(comb %*% t.default(logY))
    monoms[zeros] = 0
    monoms
  }

  monomsA = monomsFix0(A)
  M = rbind(1, monomsA)
  rhs = c(1, monomsA %*% p1)

  # If B is not NULL, add monomials of B
  if(!is.null(B)) {
    monomsB = monomsFix0(B)
    M = rbind(M, monomsB)
    rhs = c(rhs, monomsB %*% p1)
  }

  # Row-scale each constraint
  scales = apply(abs(M), 1, max)
  scales[scales == 0] = 1
  M1 = M / scales
  rhs1 = rhs / scales

  sol = lpSolve::lp("min", rep(0, n), M1, "=", rhs1, scale = 0)
  keep = which(sol$solution > 1e-10)

  res = list(status = sol$status,
             p = scale * sol$solution[keep],
             A = A[keep, , drop=FALSE])

  if(!is.null(B))
    res$B = B[keep, , drop=FALSE]

  res
}

# Stroud lumping  ---------------------------------------------------------

# Not working well: Problems with rank deficiency and numerical instability)

# A centred simplex in R^n
simplex = function(n) {
  N = n + 1
  # Vertices in R^N: e_i - (1/N)*1, for i = 1,..., N
  V = diag(N) - 1/N
  # Orthonormal basis for the subspace where sum = 0
   Q = qr.Q(qr(matrix(1, N, 1)), complete = TRUE)[, 2:N]  # Q is N x n.
  # Vertices in R^n: n x N
  coords = t(Q) %*% V
  # Scale so that, with equal weights 1/(n+1), the covariance is I_n:
  coords * sqrt(N)
}

# Not used: Non-centered simplex. See Stroud's book page 161.
simplex2 = function(n) {
  r = (n + 2 - sqrt(n+2))/((n+1)*(n+2))
  s = (n + 2 + sqrt(n+2))/((n+1)*(n+2))
  res = matrix(r, n, n+1)
  res[cbind(1:n, 2:(n+1))] = s
  res
}

stroud = function(U, P) {
  sumP = sum(P)
  P = P / sumP
  m = as.numeric(U %*% P)
  n = length(m)
  C = U %*% diag(P) %*% t(U) - m %*% t(m)
  # C = C + 1e-10*diag(n)
  # C = as.matrix(nearPD(C)$mat)
  L = t(chol(C))
  nodes = L %*% simplex(n) + m
  list(nodes = nodes, weights = rep(sumP/(n+1), n+1))
}

# Code used in lumpMutSpecial:
#
# if(FALSE && d == 2) {
#   if(length(als) <= 2*nKeep + 1) {
#     if(verbose) cat("Nothing to gain by lumping\n")
#     return(mutmat)
#   }
#
#   # Main computation: Find nodes
#   U = mutmat[lump, keep, drop = FALSE]
#   reduced = stroud(t(U), afreq[lump])
#
#   # Replace lumped alleles with nKeep + 1 nodes
#   nNodes = nKeep + 1
#   nodeNames = paste0(".lmp", 1:nNodes)
#   newAls = c(keep, nodeNames)
#
#   # Compute the new mutation matrix
#   N = nKeep + nNodes
#   newM = matrix(0, ncol = N, nrow = N, dimnames = list(newAls, newAls))
#   newM[1:nKeep, 1:nKeep] = mutmat[keep, keep]
#   newM[nKeep + (1:nNodes), 1:nKeep] = t(reduced$nodes)
#
#   # New frequency vector
#   lumpFreq = (1 - sum(afreq[keep]))/nNodes
#   newAfreq = as.numeric(c(afreq[keep], rep(lumpFreq, nNodes)))
#   names(newAfreq) = newAls
# }
