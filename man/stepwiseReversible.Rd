% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stepwiseReversible.R
\name{stepwiseReversible}
\alias{stepwiseReversible}
\title{Dawid's reversible stepwise model}
\usage{
stepwiseReversible(alleles, afreq, rate, range, maxRateOnly = FALSE)
}
\arguments{
\item{alleles}{A vector of integer integers.}

\item{afreq}{A numeric vector of allele frequencies.}

\item{rate}{A numeric mutation rate.}

\item{range}{A positive number.}

\item{maxRateOnly}{A logical, by default FALSE. See Value.}
}
\value{
A reversible stepwise mutation model with overall mutation rate equal
to \code{rate}.

If \code{maxRateOnly} is TRUE, the function returns a vector of two numbers
named \code{UW} and \code{UB}. The first of these is the maximum overall mutation
rate for a well-defined stepwise reversible mutation matrix with the given
input. The latter (UB) is the maximum rate under the additional restraint
that the model is bounded by \code{afreq}.
}
\description{
#' A reversible stepwise mutation model is created following the approach of
Dawid et al. (2002).
}
\details{
\strong{NB: This function is deprecated: Use \code{mutationMatrix(model = "dawid", ...)}
instead.}

For the stepwise reversible model, the mutation rate \eqn{r_{i,j},\,  i\neq
j} is proportional to the overall mutation rate \eqn{\lambda} for given
values of the range, the allele frequency \eqn{p_i} and n, the number of
alleles. Hence, one can determine bounds UW and UB so that the model is well
defined if \eqn{\lambda \leq UW} and bounded, i.e., \eqn{r_{i,j} \leq p_j,\,
i\neq j}, if \eqn{\lambda \leq UB}, The bounds UW and UB are computed.
}
\examples{
stepwiseReversible(alleles = 1:3,
                   afreq = c(0.2, 0.3,  0.5),
                   rate = 0.001,
                   range = 0.1)

stepwiseReversible(alleles = 1:3,
                   afreq = c(0.2, 0.3,  0.5),
                   range = 0.1,
                   maxRateOnly = TRUE)

# Model not well defined:
\dontrun{
stepwiseReversible(alleles = 1:3,
                   afreq = c(0.2, 0.3,  0.5),
                   rate = 0.7,
                   range = 0.1)
}
}
\author{
Thore Egeland.
}
