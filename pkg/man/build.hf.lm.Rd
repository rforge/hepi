\name{build.hf.lm}
\alias{build.hf.lm}
\title{ Hedonic Function Based on a Linear Model }
\description{
  This function estimates hedonic functions based on a linear regression model
  for a given dataset.
}
\usage{
build.hf.lm(learndata, full.formula, min.formula, 
            backtrans = I, rm.infl = TRUE, 
            description = NULL, return.row.labels = FALSE, 
            allow.variable.selection = TRUE)
}
\arguments{
  \item{learndata}{ A \code{data.frame} containing the training data set. }
  \item{full.formula}{ The formula of the full linear model. See Details. }
  \item{min.formula}{ If variable selection is wanted, the formula of the minimal linear model. See Details. }
  \item{backtrans}{ A backtransformation function applied to all predictions. See Details. }
  \item{rm.infl}{ A logical value indicating whether influential observations should be removed. }
  \item{description}{ A character string describing the hedonic function. }
  \item{return.row.labels}{ A logical value indicating whether the row labels of the cleaned training data should be returned.  }
  \item{allow.variable.selection}{ A logical value indicating whether variable selection should be carried out. }
}
\details{
  This function estimates a hedonic function based on a linear regression model.
  An appropriate model formula must be given in \code{full.formula}.  (See \code{\link[stats]{lm}}
  for more details about specifying formulae.)
  
  The function given in \code{backtrans} is used to backtransform any predicted
  value using the linear model and defauls to the identity function \code{I}. 
  If, for example, \code{log(price)} stands on the left-hand side of the model
  formula, any predicted value needs to be transformed with the exponential function
  to a valid price.  This can be accomplished by indicating \code{backtrans = exp}.
  
  If \code{rm.infl} is \code{TRUE}, influential observations having
  \deqn{\mathrm{DFFITS}_i>2\sqrt{\frac{K}{N}}}{|DFFITS_i| > 2 * sqrt(K/N)}
  with \eqn{K} being the number of exogenous variables and \eqn{N} the dimension
  of the learning data set are removed before fitting the final model.  There, the
  \eqn{\mathrm{DFFITS}_i}{DFFITS_i} values are calculated based on the residuals
  of a first fit of a linear model using the model formula \code{full.formula}.
  
  If \code{allow.variable.selection} is \code{TRUE}, a stepwise model selection
  based on exact AIC is carried out (see \code{\link[MASS]{stepAIC}} for more details).  In this case,
  \code{full.formula} acts as upper and \code{min.formula} as lower limit of the search
  algorithm.  If \code{allow.variable.selection} is \code{FALSE} the hedonic function
  is estimated using exactly the formula given in \code{full.formula}.
  
  In \code{description}, a character string describing the hedonic function may be given
  which is saved within the returned \code{"\linkS4class{hedonic.function}"} object.
}
\value{
  If \code{return.row.labels == FALSE}, the function returns a \code{"\linkS4class{hedonic.function}"} object representing the fitted regression model.
  
  If \code{return.row.labels == TRUE}, the function returns a list with following elements:
  \item{hf }{The resulting \code{"\linkS4class{hedonic.function}"} object.}
  \item{row.labels }{ A vector containing the row labels of the cleaned training data set. }
}
\references{   
  Beer, M. (2007) \emph{Hedonic Elementary Price Indices: Axiomatic Foundation and Estimation Techniques.}
  PhD thesis, University of Fribourg Switzerland, \url{http://www.michael.beer.name/phdthesis}.
 }
\author{ Michael Beer \email{r-hepi@michael.beer.name} }
\seealso{ \code{\link{build.hf.lm.split}} }
\examples{
data(boston, package = "spdep")

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

is.applicable.hf(hf0, boston.c)
summary(hf0(boston.c))

plot(boston.c$MEDV, hf0(boston.c), xlab = "Observed", ylab = "Predicted")
abline(0,1)

hf1 <- build.hf.lm(
    learndata = boston.c, 
    full.formula = log(MEDV) ~ CRIM + ZN + INDUS + CHAS + 
      I(NOX^2) + I(RM^2) + AGE + log(DIS) + log(RAD) + TAX + 
      PTRATIO + B + log(LSTAT), 
    min.formula = log(MEDV) ~ 1, 
    backtrans = exp, 
    rm.infl =  FALSE, 
    description = NULL, 
    return.row.labels = FALSE, 
    allow.variable.selection = TRUE)
summary(hf1(boston.c))

}
\keyword{ regression }%

