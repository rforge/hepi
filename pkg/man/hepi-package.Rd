\name{hepi-package}
\alias{hepi-package}
\docType{package}
\title{
Functions for Estimating Hedonic Elementary Price Indices
}
\description{
This package provides a set of functions for estimating
hedonic functions and hedonic elementary price indices.
}
\details{
\tabular{ll}{
Package: \tab hepi\cr
Date: \tab 2008-11-10\cr
Version: \tab 1.0-1\cr
License: \tab GPL (version 2 or later) See file LICENCE.\cr
Depends: \tab R (>= 2.2.0)\cr
Suggests: \tab MASS, spdep\cr
LazyLoad: \tab yes\cr
SaveImage: \tab yes\cr
URL: \tab http://www.michael.beer.name/r-hepi\cr
}

Index:
\describe{
\item{\code{\link{hedonic.function-class}}}{Class \code{"hedonic.function"}}
\item{\code{\link{hedonic.function}}}{Constructor for a \code{"hedonic.function"} Object}
\item{\code{\link{analyse.in.hf}}}{Evaluate Expressions Within the Environment of a Hedonic Function}
\item{\code{\link{is.applicable.hf}}}{Test Applicability of Hedonic Function to Data}
\item{\code{\link{build.hf.lm}}}{Hedonic Function Based on a Linear Model}
\item{\code{\link{build.hf.lm.split}}}{Hedonic Function Based on a List of Linear Models}
\item{\code{\link{hepi}}}{Bilateral Hedonic Elementary Price Indices}
}
}
\author{Michael Beer \email{r-hepi@michael.beer.name}}
}
\references{   
  Beer, M. (2007) \emph{Hedonic Elementary Price Indices: Axiomatic Foundation and Estimation Techniques.}
  PhD thesis, University of Fribourg Switzerland, \url{http://www.michael.beer.name/phdthesis}.
 }
\keyword{ package }

