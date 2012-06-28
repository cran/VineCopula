\name{BiCopTau2Par}           
\alias{BiCopTau2Par}            

\title{Parameter of a bivariate copula for a given Kendall's tau value}

\description{
This function computes the parameter of a one parameter bivariate copula for a given value of Kendall's tau.
}

\usage{
BiCopTau2Par(family, tau)
}

\arguments{
  \item{tau}{Kendall's tau value (numeric in [-1,1]).}
  \item{family}{An integer defining the bivariate copula family:\cr
		\code{0} = independence copula \cr
    \code{1} = Gaussian copula \cr
    \code{3} = Clayton copula \cr
    \code{4} = Gumbel copula \cr
    \code{5} = Frank copula \cr
    \code{6} = Joe copula \cr 
		\code{13} = rotated Clayton copula (180 degrees; ``survival Clayton'') \cr
		\code{14} = rotated Gumbel copula (180 degrees; ``survival Gumbel'') \cr
		\code{16} = rotated Joe copula (180 degrees; ``survival Joe'') \cr
		\code{23} = rotated Clayton copula (90 degrees) \cr
		\code{24} = rotated Gumbel copula (90 degrees) \cr
		\code{26} = rotated Joe copula (90 degrees) \cr
		\code{33} = rotated Clayton copula (270 degrees) \cr
		\code{34} = rotated Gumbel copula (270 degrees) \cr
		\code{36} = rotated Joe copula (270 degrees)\cr
    Note that two parameter bivariate copula families cannot be used.}
}

\value{
Parameter corresponding to the bivariate copula family and the value of Kendall's tau (\eqn{\tau}).
\tabular{ll}{
No. \tab Parameter \cr
\code{1, 2} \tab \eqn{\sin(\tau \frac{\pi}{2})}{sin(\tau \pi/2)} \cr
\code{3, 13} \tab \eqn{\max(0,2\frac{\tau}{1-\tau})}{max(0,2\tau/(1-\tau))} \cr
\code{4, 14} \tab \eqn{\max(1,\frac{1}{1-\tau})}{max(1,1/(1-\tau))} \cr
\code{5} \tab no closed form expression (numerical inversion) \cr
\code{6, 16} \tab no closed form expression (numerical inversion) \cr
\code{23, 33} \tab \eqn{\max(0,2\frac{\tau}{1+\tau})}{max(0,2\tau/(1+\tau))} \cr
\code{24, 34} \tab \eqn{\min(-1,-\frac{1}{1+\tau})}{min(-1,-1/(1+\tau))} \cr
\code{26, 36} \tab no closed form expression (numerical inversion)
}

}

\author{Jakob Stoeber, Eike Brechmann}

\references{
Joe, H. (1997).
Multivariate Models and Dependence Concepts.
Chapman and Hall, London.

Czado, C., U. Schepsmeier, and A. Min (2012).
Maximum likelihood estimation of mixed C-vines with application to exchange rates.
Statistical Modelling, 12(3), 229-255.
}

\seealso{\code{\link{BiCopTau2Par}}}

\examples{
## Example 1: Gaussian copula
tt1 = BiCopTau2Par(1,0.5)

# transform back
BiCopPar2Tau(1,tt1)


## Example 2: Clayton copula
BiCopTau2Par(3,0.4)
}

