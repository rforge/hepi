\name{build.hf.lm.split}
\alias{build.hf.lm.split}
\title{ Hedonic Function Based on a List of Linear Models }
\description{
  This function estimates hedonic functions based on linear regression models
  for a given dataset.  Individual sub-models are fit to subsets of the data 
  split up according to a factor variable.
}
\usage{
build.hf.lm.split(learndata, split.var, full.formula, min.formula,
    backtrans = I, rm.infl = FALSE, description = NULL, 
    return.row.labels = FALSE, allow.variable.selection = TRUE, 
    use.overall.hf = TRUE, split.threshold = 100)
}
\arguments{
  \item{learndata}{ A \code{data.frame} containing the training data set. }
  \item{split.var}{ The name of the factor variable used to split the data into subsets. See Details. }
  \item{full.formula}{ The formula of the full linear model. See Details. }
  \item{min.formula}{ If variable selection is wanted, the formula of the minimal linear model. See Details. }
  \item{backtrans}{ A backtransformation function applied to all predictions. See Details. }
  \item{rm.infl}{ A logical value indicating whether influential observations should be removed. }
  \item{description}{ A character string describing the hedonic function. }
  \item{return.row.labels}{ A logical value indicating whether the row labels of the cleaned training data should be returned.  }
  \item{allow.variable.selection}{ A logical value indicating whether variable selection should be carried out. }
  \item{use.overall.hf}{ A logical value indicating whether an overall model should be fit to the whole data set.  }
  \item{split.threshold}{ The minimal number of observations required for fitting any sub-model. }
}
\details{
  This function estimates a hedonic function based on linear regression models. 
  In contrast to \code{\link{build.hf.lm}}, however, individual linear models are fit to 
  several subsets of the data. These subsets are determined through a factor variable named
  \code{split.var} which needs to be contained in \code{learndata}. The minimal size a
  subset needs to have in order to fit a linear model is given by \code{split.threshold}. 
  
  If \code{use.overall.hf} is \code{TRUE}, an overall model for the whole data set is fit
  and stored additionally in order to predict prices for characteristics vectors belonging
  to categories of \code{split.var} where less than \code{split.threshold} observations are
  available in the learning data set.
  
  See the documentation of \code{\link{build.hf.lm}} for an explanation of the other 
  arguments of the function. Removal of influential observations and variable selection, if required,
  is carried out for each sub-model individually.
  
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
\seealso{ \code{\link{build.hf.lm}} }
\examples{
data(boston, package = "spdep")

hf0 <- build.hf.lm.split(
    learndata = boston.c,
    split.threshold = 15, 
    split.var = "TOWN",
    full.formula = log(MEDV) ~ CRIM + ZN + INDUS + CHAS + 
      I(NOX^2) + I(RM^2) + AGE + log(DIS) + log(RAD) + TAX + 
      PTRATIO + B + log(LSTAT), 
    backtrans = exp, 
    rm.infl = FALSE, 
    description = NULL, 
    return.row.labels = FALSE, 
    allow.variable.selection = FALSE)

is.applicable.hf(hf0, boston.c)
summary(hf0(boston.c))

plot(boston.c$MEDV, hf0(boston.c), xlab = "Observed", ylab = "Predicted")
abline(0,1)

hf1 <- build.hf.lm.split(
    learndata = boston.c, 
    split.var = "TOWN",
    split.threshold = 15, 
    full.formula = log(MEDV) ~ CRIM + ZN + INDUS + CHAS + 
      I(NOX^2) + I(RM^2) + AGE + log(DIS) + log(RAD) + TAX + 
      PTRATIO + B + log(LSTAT), 
    min.formula = log(MEDV) ~ 1, 
    backtrans = exp, 
    rm.infl = FALSE, 
    description = NULL, 
    return.row.labels = FALSE, 
    allow.variable.selection = TRUE)
summary(hf1(boston.c))
}
\keyword{ regression }%
