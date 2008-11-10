\name{hedonic.function}
\alias{hedonic.function}
\title{ Constructor for a "hedonic.function" Object }
\description{
  This function constructs objects of class \code{"\linkS4class{hedonic.function}"}.
}
\usage{
hedonic.function(hf, characteristics.names, env = NULL, 
                 call = match.call(), description = "")
}
\arguments{
  \item{hf}{ The hedonic function. }
  \item{characteristics.names}{ A vector of characteristics names. }
  \item{env}{ An optional list of objects needed internally for the evaluation of \code{hf}. }
  \item{call}{ An object of class \code{"\linkS4class{call}"}, describing how \code{hf} was created. }
  \item{description}{ A textual description of the hedonic function \code{hf}. }
}
\details{
  A hedonic function predicts a price value for an item based on a list of
  characteristics expressions.  \code{hf} is thus an object of class \code{"\linkS4class{function}"}
  having at least the named argument \code{data} representing the data frame for which new prices are to be estimated.

  \code{characteristics.names} should contain the names of all variables
  that must be available in \code{data} in order to provide valid price predictions.
  It can be used to check whether the returned hedonic function is applicable for
  a given data set.
}
\value{
  An object of class \code{"\linkS4class{hedonic.function}"}.
}
\author{ Michael Beer \email{r-hepi@michael.beer.name} }
\seealso{ \code{\link{hedonic.function-class}} }
\examples{
## Build hedonic function from training data set
build.hf.loglin <- function(traindata) {
    model <- lm(log(price) ~ ., traindata)

    hf <- function(data) {
        predict(model, newdata = data)
    }

    hedonic.function(
        hf = hf,
        characteristics.names = all.vars(formula(model)),
        env = list(model = model),
        call = match.call(),
        description = "Semi-logarithmic hedonic function."
    )
}
}
\keyword{ methods }%

