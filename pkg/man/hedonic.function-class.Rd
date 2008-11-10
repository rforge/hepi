\name{hedonic.function-class}
\docType{class}
\alias{hedonic.function-class}
\title{Class "hedonic.function" }
\description{ Hedonic functions for the use in price statistics. }
\section{Objects from the Class}{
  A hedonic function predicts a price value for an item based on a list of
  characteristics expressions.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{ An object of class \code{"\linkS4class{function}"}
  having at least the named argument \code{data} representing the data frame for
  which new prices are to be estimated. }
    \item{\code{characteristics.names}:}{A vector containing the names of all variables
  that must be available in \code{data} in order to provide valid price predictions. }
    \item{\code{call}:}{An object of class \code{"\linkS4class{call}"}, describing how \code{.Data} was created. }
    \item{\code{description}:}{A textual description of the hedonic function \code{.Data}. }
  }
}
\section{Extends}{
Class \code{"function"}, from data part.
Class \code{"OptionalFunction"}, by class \code{"function"}.
Class \code{"PossibleMethod"}, by class \code{"function"}.
}
\section{Methods}{
No methods defined with class "hedonic.function" in the signature.
}
\author{ Michael Beer \email{r-hepi@michael.beer.name} }
\seealso{
  \code{\link{hedonic.function}}
}
\keyword{classes}
