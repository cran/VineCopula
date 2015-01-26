\name{BiCopHfunc}           
\alias{BiCopHfunc}

\title{Conditional Distribution Function of a Bivariate Copula}

\description{
This function evaluates the conditional distribution function (h-function) of a given parametric bivariate copula.
}

\usage{
BiCopHfunc(u1, u2, family, par, par2 = 0)
}

\arguments{
  \item{u1, u2}{Numeric vectors of equal length with values in [0,1].}
  \item{family}{An integer defining the bivariate copula family: \cr
		\code{0} = independence copula \cr
		\code{1} = Gaussian copula \cr
	        \code{2} = Student t copula (t-copula) \cr
	        \code{3} = Clayton copula \cr
	        \code{4} = Gumbel copula \cr
	        \code{5} = Frank copula \cr
	        \code{6} = Joe copula \cr 
		\code{7} = BB1 copula \cr
		\code{8} = BB6 copula \cr
		\code{9} = BB7 copula \cr
		\code{10} = BB8 copula \cr
		\code{13} = rotated Clayton copula (180 degrees; ``survival Clayton'') \cr
		\code{14} = rotated Gumbel copula (180 degrees; ``survival Gumbel'') \cr
		\code{16} = rotated Joe copula (180 degrees; ``survival Joe'') \cr 
		\code{17} = rotated BB1 copula (180 degrees; ``survival BB1'')\cr
		\code{18} = rotated BB6 copula (180 degrees; ``survival BB6'')\cr
		\code{19} = rotated BB7 copula (180 degrees; ``survival BB7'')\cr
		\code{20} = rotated BB8 copula (180 degrees; ``survival BB8'')\cr
		\code{23} = rotated Clayton copula (90 degrees) \cr
		\code{24} = rotated Gumbel copula (90 degrees) \cr
		\code{26} = rotated Joe copula (90 degrees) \cr
		\code{27} = rotated BB1 copula (90 degrees) \cr
		\code{28} = rotated BB6 copula (90 degrees) \cr
		\code{29} = rotated BB7 copula (90 degrees) \cr
		\code{30} = rotated BB8 copula (90 degrees) \cr
		\code{33} = rotated Clayton copula (270 degrees) \cr
		\code{34} = rotated Gumbel copula (270 degrees) \cr
		\code{36} = rotated Joe copula (270 degrees) \cr
		\code{37} = rotated BB1 copula (270 degrees) \cr
		\code{38} = rotated BB6 copula (270 degrees) \cr
		\code{39} = rotated BB7 copula (270 degrees) \cr
		\code{40} = rotated BB8 copula (270 degrees) \cr
    \code{104} = Tawn type 1 copula \cr
    \code{114} = rotated Tawn type 1 copula (180 degrees) \cr
    \code{124} = rotated Tawn type 1 copula (90 degrees)  \cr
    \code{134} = rotated Tawn type 1 copula (270 degrees) \cr
    \code{204} = Tawn type 2 copula  \cr
    \code{214} = rotated Tawn type 2 copula (180 degrees) \cr
    \code{224} = rotated Tawn type 2 copula (90 degrees)  \cr
    \code{234} = rotated Tawn type 2 copula (270 degrees) \cr
		}
  \item{par}{Copula parameter.}
  \item{par2}{Second parameter for bivariate copulas with two parameters (t, BB1, BB6, BB7, BB8, Tawn type 1 and type 2; default: \code{par2 = 0}).}
}

\details{
The h-function is defined as the conditional distribution function of a bivariate copula, i.e.,
\deqn{
h(u|v,\boldsymbol{\theta}) := F(u|v) =
\frac{\partial C(u,v)}{\partial v},
}{
h(u|v,\theta) := F(u|v) =
\partial C(u,v) / \partial v,
}
where \eqn{C} is a bivariate copula distribution function with parameter(s) \eqn{\boldsymbol{\theta}}{\theta}.  
For more details see Aas et al. (2009).
}


\value{
\item{hfunc1}{Numeric vector of the conditional distribution function (h-function) evaluated at \code{u2} given \code{u1}, i.e., \eqn{h(\code{u2}|\code{u1},\boldsymbol{\theta})}{h(u2|u1,\theta)}.}
\item{hfunc2}{Numeric vector of the conditional distribution function (h-function) evaluated at \code{u1} given \code{u2}, i.e., \eqn{h(\code{u1}|\code{u2},\boldsymbol{\theta})}{h(u1|u2,\theta)}.}
}

\references{ 
Aas, K., C. Czado, A. Frigessi, and H. Bakken (2009).
Pair-copula constructions of multiple dependence.
Insurance: Mathematics and Economics 44 (2), 182-198.
}


\author{Ulf Schepsmeier}

\seealso{\code{\link{BiCopPDF}}, \code{\link{BiCopCDF}}, \code{\link{RVineLogLik}}, \code{\link{RVineSeqEst}}}

\examples{
# load data set
data(daxreturns)

# h-functions of the Gaussian copula
h1 <- BiCopHfunc(daxreturns[,2], daxreturns[,1], 1, 0.5)
\dontshow{
h1
}
}
