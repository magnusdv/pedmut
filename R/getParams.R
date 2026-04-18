#' Get model parameters
#'
#' Extract model parameters of a mutation matrix/model.
#'
#' @param mut A [mutationModel()] or [mutationMatrix()].
#' @param params A vector contain some or all of the words "model", "rate",
#'   "range", "rate2", "seed". If NULL (default), all present parameters are
#'   included.
#' @param format A numeric code indicating the wanted output format. See Value.
#' @param sep A separator character used to paste male and female values.
#'   Ignored unless `format = 3`.
#'
#' @returns When `mut` is a `mutationModel`, the output format is dictated by the `format`
#'   option, with the following possibilities:
#'
#'   1. A data frame with 2 rows labelled 'female' and 'male'.
#'
#'   2. A data frame with 1 row and female/male columns suffixed by .F/.M
#'   respectively.
#'
#'   3. A data frame with 1 row, in which female/male values are pasted together
#'   (separated with `sep`) if different.
#'
#'   4. A list version of format 1, suitable for programmatic use with `mutationModel()`.
#'
#'   If `mut` is a `mutationMatrix` the output is always a data frame with 1 row.
#'
#' @examples
#' M = mutationModel("equal", 1:2, rate = list(female = 0.2, male = 0.1))
#' getParams(M)
#' getParams(M, format = 2)
#' getParams(M, format = 3)
#' getParams(M, format = 3, sep = "|")
#'
#' pars = getParams(M, format = 4)
#' pars$alleles = 1:2  # Not part of the parameters, but needed to reconstruct the model.
#' M2 = do.call(mutationModel, pars)
#' stopifnot(identical(M, M2))
#'
#' @export
getParams = function(mut, params = NULL, format = 1, sep = "/") {
  if(isMutationModel(mut)) {

    parlist = lapply(mut, function(mat) getParams(mat, params = params))
    params = names(parlist$female)

    res = switch(format,
            # 1
            do.call(rbind, parlist),
            # 2
            {
              names(parlist$female) = paste0(names(parlist$female), ".F")
              names(parlist$male) = paste0(names(parlist$male), ".M")
              res = cbind(parlist$female, parlist$male)
              # Interleave
              npar = length(params)
              idx = as.numeric(rbind(1:npar, 1:npar + npar))
              res[idx]
            },
            # 3
            {
              vals = lapply(params, function(p) {
                valF = parlist$female[[p]]
                valM = parlist$male[[p]]
                if(identical(valF, valM)) valF else paste(valF, valM, sep = sep)
              })
              names(vals) = params
              as.data.frame(vals)
            },
            # 4
            {
              vals = lapply(params, \(v) lapply(parlist, `[[`, v))
              names(vals) = params
              vals
            }
          )
    return(res)
  }

  attrs = attributes(mut)

  if(is.null(params))
    params = intersect(names(attrs), c("model", "rate", "range", "rate2", "seed"))
  else {
    unkn = setdiff(params, c("model", "rate", "range", "rate2", "seed"))
    if(length(unkn))
      stop2("Unknown param: ", unkn)
  }

  vals = lapply(params, function(p) attrs[[p]] %||% NA)
  names(vals) = params
  as.data.frame(vals)
}
