\name{hepi}
\alias{hepi}
\alias{hepi.jevons}
\alias{hepi.carli}
\alias{hepi.dutot}
\alias{hepi.hcarli}
\alias{hepi.hdutot}
\title{Bilateral Hedonic Elementary Price Indices}
\description{
  This function estimates bilateral hedonic elementary price indices
  based on estimations of the hedonic function
  in the base and the current period as well as a reference sample of 
  quality characteristics.
}
\usage{
hepi(hf0, hf1, M, 
     type = c("jevons", "dutot", "carli", "hdutot", "hcarli"), 
     na.rm = TRUE, debug = FALSE)

hepi.jevons(hf0, hf1, M)
hepi.carli(hf0, hf1, M)
hepi.dutot(hf0, hf1, M)
hepi.hcarli(hf0, hf1, M)
hepi.hdutot(hf0, hf1, M)
}
\arguments{
  \item{hf0}{ The base period hedonic function estimate. Must be of class \code{"\linkS4class{hedonic.function}"}. }
  \item{hf1}{ The current period hedonic function estimate. Must be of class \code{"\linkS4class{hedonic.function}"}. }
  \item{M}{ The reference sample. }
  \item{type}{ The type of the index estimator(s) to be used. Can be a vector if estimates of several types are requested.}
  \item{na.rm}{ A logical value indicating whether observations containing \code{NA} values should be stripped before the computation proceeds. }
  \item{debug}{ A logical value indicating whether predicted prices should be returned for debugging purposes. }
}
\details{
    This function yields an estimate of a bilateral hedonic price index.  Inputs are the two estimated
    hedonic functions \code{hf0} and \code{hf1} of the base and current period respectively.  Both of these
    must be of class \code{"\linkS4class{hedonic.function}"}.  (See \code{\link{hedonic.function}} for a constructor 
    of a \code{hedonic.function} object.)
    
    The third input is the reference sample \code{M} to be used for the estimation of the index.  This
    is usually a data frame containing \eqn{N}{N} characteristics vectors to which both hedonic functions 
    are applicable.
    
    The \code{type} argument lets one choose the index formula to be used (and yet the index to be estimated).  
    Currently, we implemented five alternative estimators, namely the
    \deqn{\sqrt[n]{\prod_{i=1}^N\frac{\hat{h}^1(m_n)}{\hat{h}^0(m_n)}}\qquad\mbox{(Jevons),}}{Jevons}
    \deqn{\frac{\sum_{n=1}^N
            \hat{h}^1(m_n)}{\sum_{n=1}^N
            \hat{h}^0(m_n)}\qquad\mbox{(Dutot),}}{Dutot}
    \deqn{\frac{1}{N}\sum_{n=1}^N\frac{
            \hat{h}^1 (m_n)}{\hat{h}^0(m_n)}\qquad\mbox{(Carli),}}{Carli}
    \deqn{\frac{\left(\sum_{n=1}^N
            \bigl(\hat{h}^1(m_n)\bigr)^{-1}\right)^{-1}}{\left(\sum_{n=1}^N
            \bigl(\hat{h}^0(m_n)\bigr)^{-1}\right)^{-1}}\qquad\mbox{(Harmonic Dutot) and}}{Harmonic Dutot and}
    \deqn{\left(\frac{1}{N}\sum_{n=1}^N\left(\frac{ \hat{h}^1
            (m_n)}{\hat{h}^0(m_n)}\right)^{\!\!-1}\right)^{\!\!-1}\qquad\mbox{(Harmonic Carli)}}{Harmonic Carli}
    formulae.
    Details can be found in Chapter 6 of the reference mentioned below.
}
\value{
    If \code{debug == FALSE}, this function returns a vector with the same length as \code{type} containing the index estimates. 
    They are returned in the same order as given by \code{type}.
    
    If \code{debug == TRUE}, this function returns a list with the following entries
  \item{index }{The vector of index estimates as above.}
  \item{p0hat }{The vector of predicted prices \eqn{\hat{\boldsymbol{p}}^0=\hat{h}^0(\boldsymbol{M})}{p^0=h^0(M)} in the base period.}
  \item{p1hat }{The vector of predicted prices \eqn{\hat{\boldsymbol{p}}^1=\hat{h}^1(\boldsymbol{M})}{p^1=h^1(M)} in the current period.}
  \item{ratios }{The vector of price ratios \eqn{\hat{p}^1_n/\hat{p}^0_n}{p_n^1/p_n^0} \eqn{(n=1,\dots,N)}.}
}
\references{   
  Beer, M. (2007) \emph{Hedonic Elementary Price Indices: Axiomatic Foundation and Estimation Techniques.}
  PhD thesis, University of Fribourg Switzerland, \url{http://www.michael.beer.name/phdthesis}.
 }
\author{ Michael Beer \email{r-hepi@michael.beer.name} }
\keyword{ manip }%
