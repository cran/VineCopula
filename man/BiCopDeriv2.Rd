\name{BiCopDeriv2}           
\alias{BiCopDeriv2}

\title{Second Derivatives of a Bivariate Copula Density}

\description{
This function evaluates the second derivative of a given parametric bivariate copula density with respect to its parameter(s) and/or its arguments.
}

\usage{
BiCopDeriv2(u1, u2, family, par, par2 = 0, deriv = "par")
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
  	\code{"par"} = second derivative with respect to the first parameter (default)\cr
  	\code{"par2"} = second derivative with respect to the second parameter (only available for the t-copula) \cr
  	\code{"u1"} = second derivative with respect to the first argument \code{u1} \cr
  	\code{"u2"} = second derivative with respect to the second argument \code{u2} \cr
	\code{"par1par2"} = second derivative with respect to the first and second parameter (only available for the t-copula) \cr
	\code{"par1u1"} = second derivative with respect to the first parameter and the first argument \cr
	\code{"par2u1"} = second derivative with respect to the second parameter and the first argument (only available for the t-copula) \cr
	\code{"par1u2"} = second derivative with respect to the first parameter and the second argument \cr
	\code{"par2u2"} = second derivative with respect to the second parameter and the second argument (only available for the t-copula) \cr
  	}
}

\value{
A numeric vector of the second bivariate copula derivative with respect to \code{deriv} evaluated at \code{u1} and \code{u2} with parameter(s) \code{par} and \code{par2}.
}

\references{
Schepsmeier, U. and J. Stoeber (2012).
Derivatives and Fisher information of bivariate copulas.
Statistical Papers.
\url{http://link.springer.com/article/10.1007/s00362-013-0498-x}.
}

\author{Ulf Schepsmeier, Jakob Stoeber}

\seealso{\code{\link{RVineGrad}}, \code{\link{RVineHessian}}, \code{\link{BiCopDeriv}}, \code{\link{BiCopHfuncDeriv}}}

\examples{
# simulate from a bivariate t-copula
simdata <- BiCopSim(300, 2, -0.7, par2 = 4)

# second derivative of the bivariate t-copula with respect to the first parameter
u1 <- simdata[,1]
u2 <- simdata[,2]
BiCopDeriv2(u1, u2, 2, -0.7, par2 = 4, deriv = "par")
}

