% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeReversible.R
\name{makeReversible}
\alias{makeReversible}
\title{Transformations to reversibility}
\usage{
makeReversible(
  mutmat,
  method = c("BA", "MH", "PR"),
  adjust = FALSE,
  afreq = NULL
)
}
\arguments{
\item{mutmat}{A mutation matrix.}

\item{method}{A character indicating which transformation to use. Either "BA"
(Barker), "MH" (Metropolis-Hastings) or "PR" (preserved rate).}

\item{adjust}{Logical. If TRUE, the overall mutation rate is adjusted to
preserve the original rate; see \code{\link[=adjustRate]{adjustRate()}}. Not relevant for the "PR"
method, which by construction always preserves the overall rate.}

\item{afreq}{A vector of allele frequencies. Extracted from \code{mutmat} if not
provided.}
}
\value{
A reversible mutation matrix with the same allele frequencies.
}
\description{
This function implements three methods for transforming a mutation matrix to
a reversible model. All methods are based on Metropolis-Hastings proposal
functions.
}
\examples{
m = mutationMatrix("equal", afreq = c(a=0.2, b=0.3, c=0.5), rate = 0.2)
makeReversible(m)
makeReversible(m, adjust = TRUE)

makeReversible(m, "MH", adjust = TRUE)
# makeReversible(m, "PR") # not well-defined!

}
