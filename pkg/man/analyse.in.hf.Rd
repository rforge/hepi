\name{analyse.in.hf}
\alias{analyse.in.hf}
\title{ Evaluate Expressions Within the Environment of a Hedonic Function }
\description{
  This function evaluates \R expressions within the environment of an object of class \code{"\linkS4class{hedonic.function}"}.
}
\usage{
analyse.in.hf(expr, hf)
}
\arguments{
  \item{expr}{ the expression to be evaluated }
  \item{hf}{ an object of class \code{"\linkS4class{hedonic.function}"} }
}
\details{
  Hedonic functions are often predictions from a regression model. This function allows 
  to easily access the elements in the environment of the hedonic function for further analysis.
}
\author{ Michael Beer \email{r-hepi@michael.beer.name} }
\seealso{ \code{\link{hedonic.function-class}}, \code{\link{hedonic.function}} }
\examples{
data(boston, package="spdep")

hf0 <- build.hf.lm(
    learndata = boston.c, 
    full.formula = log(MEDV) ~ CRIM + ZN + INDUS + CHAS + 
      I(NOX^2) + I(RM^2) + AGE + log(DIS) + log(RAD) + TAX + 
      PTRATIO + B + log(LSTAT), 
    backtrans = exp, 
    rm.infl =  FALSE, 
    description = NULL, 
    return.row.labels = FALSE, 
    allow.variable.selection = FALSE)

analyse.in.hf(ls(), hf0)
analyse.in.hf(summary(learndata), hf0)
analyse.in.hf(summary(hf.model), hf0)
}
\keyword{ methods }%
