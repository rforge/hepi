\name{is.applicable.hf}
\alias{is.applicable.hf}
\title{Test Applicability of Hedonic Function to Data}
\description{
  This function tests whether a \code{"\linkS4class{hedonic.function}"} object is applicable to a given data set.
}
\usage{
is.applicable.hf(hf, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hf}{ an object of class \code{"hedonic.function"} }
  \item{data}{ a named list or data frame }
}
\details{
  Hedonic functions are usually estimated for a very specific data structure and cannot
  be applied to other situations.  This function checks whether the given hedonic function
  \code{hf} is applicable to a certain data set \code{data}.  This is done by checking whether
  \code{data} contains at least the variables mentioned in the slot \code{characteristics.names}  of the
  object \code{hf} of class \code{"\linkS4class{hedonic.function}"}.
}
\value{
  \code{is.applicable.hf} returns \code{TRUE} if \code{data} is suitable for \code{hf} and \code{FALSE} otherwise.
}
\author{ Michael Beer \email{r-hepi@michael.beer.name} }
\seealso{ \code{\link{hedonic.function-class}}, \code{\link{hedonic.function}} }
\keyword{ methods }%
