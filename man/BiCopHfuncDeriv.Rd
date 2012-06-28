\name{BiCopHfuncDeriv}           
\alias{BiCopHfuncDeriv}

\title{Derivatives of the h-function of a bivariate copula}

\description{
This function evaluates the derivative of a given conditional parametric bivariate copula (h-function) with respect to its parameter(s) or one of its arguments.
}

\usage{
BiCopHfuncDeriv(u1, u2, family, par, par2=0, deriv="par")
}

\arguments{
  \item{u1,u2}{Numeric vectors of equal length with values in [0,1].}
  \item{family}{An integer defining the bivariate copula family: \cr
		\code{0} = independence copula \cr
	        \code{1} = Gaussian copula \cr
	        \code{2} = Student t copula (t-copula) \cr
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
		\code{36} = rotated Joe copula (270 degrees) \cr
		}
  \item{par}{Copula parameter.}
  \item{par2}{Second parameter for bivariate t-copula; default: \code{par2 = 0}.}
  \item{deriv}{Derivative argument \cr
  	\code{"par"} = derivative with respect to the first parameter (default)\cr
  	\code{"par2"} = derivative with respect to the second parameter (only available for the t-copula) \cr
  	\code{"u2"} = derivative with respect to the second argument \code{u2} \cr
  	}
}

\value{
A numeric vector of the conditional bivariate copula derivative with respect to \code{deriv} evaluated at \code{u1} and \code{u2} with parameter(s) \code{par} and \code{par2}.
}

\references{
Schepsmeier, U. and J. Stoeber (2012).
Derivatives and Fisher information of bivariate copulas.
Submitted for publication.
\url{http://mediatum.ub.tum.de/node?id=1106541}.
}

\author{Ulf Schepsmeier}

\seealso{\code{\link{RVineGrad}}, \code{\link{RVineHessian}}, \code{\link{BiCopDeriv2}}, \code{\link{BiCopDeriv2}}, \code{\link{BiCopHfuncDeriv}}}

\examples{
# simulate from a bivariate t-copula
simdata = BiCopSim(300,2,-0.7,par2=4)

# derivative of the conditional bivariate t-copula 
# with respect to the first parameter
u1 = simdata[,1]
u2 = simdata[,2]
BiCopHfuncDeriv(u1,u2,2,-0.7,par2=4, deriv="par")
}

