\name{BiCopPar2TailDep}       
\alias{BiCopPar2TailDep} 

\title{Tail dependence coefficients of a bivariate copula}

\description{
This function computes the theoretical tail dependence coefficients of a bivariate copula for given parameter values. 
}

\usage{
BiCopPar2TailDep(family, par, par2=0)
}

\arguments{
  \item{family}{An integer defining the bivariate copula family:\cr
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
  \item{par2}{Second parameter for the two parameter t-, BB1, BB6, BB7, BB8, Tawn type 1 and type 2 copulas (default: \code{par2 = 0}).}
}                                                                                  

\value{
  \item{lower}{Lower tail dependence coefficient of the given bivariate copula family \eqn{C}:
  \deqn{
  \lambda_L = \lim_{u\searrow 0}\frac{C(u,u)}{u}
  }{
	\lambda_L = lim_{u->0} C(u,u)/u
  }
  }
  \item{upper}{Upper tail dependence coefficient of the given bivariate copula family \eqn{C}:
  \deqn{
  \lambda_U = \lim_{u\nearrow 1}\frac{1-2u+C(u,u)}{1-u}
  }{
	\lambda_U = lim_{u->1}(1-2u+C(u,u))/(1-u)
  }
  }
Lower and upper tail dependence coefficients for bivariate copula families and parameters
(\eqn{\theta} for one parameter families and the first parameter of the t-copula with \eqn{\nu} degrees of freedom,
\eqn{\theta} and \eqn{\delta} for the two parameter BB1, BB6, BB7 and BB8 copulas) are given in the following table. 
\tabular{lll}{
No. \tab Lower tail dependence \tab Upper tail dependence \cr
\code{1} \tab - \tab - \cr
\code{2} \tab \eqn{2t_{\nu+1}\left(-\sqrt{\nu+1}\sqrt{\frac{1-\theta}{1+\theta}}\right)}{2t_{\nu+1}(-\sqrt{\nu+1}\sqrt{(1-\theta)/(1+\theta)})} \tab \eqn{2t_{\nu+1}\left(-\sqrt{\nu+1}\sqrt{\frac{1-\theta}{1+\theta}}\right)}{2t_{\nu+1}(-\sqrt{\nu+1}\sqrt{(1-\theta)/(1+\theta)})} \cr
\code{3} \tab \eqn{2^{-1/\theta}} \tab - \cr
\code{4} \tab - \tab \eqn{2-2^{1/\theta}} \cr
\code{5} \tab - \tab - \cr
\code{6} \tab - \tab \eqn{2-2^{1/\theta}} \cr
\code{7} \tab \eqn{2^{-1/(\theta\delta)}} \tab \eqn{2-2^{1/\delta}} \cr
\code{8} \tab - \tab \eqn{2-2^{1/(\theta\delta)}} \cr
\code{9} \tab \eqn{2^{-1/\delta}} \tab \eqn{2-2^{1/\theta}} \cr
\code{10} \tab - \tab \eqn{2-2^{1/\theta}} if \eqn{\delta=1} otherwise 0 \cr
\code{13} \tab - \tab \eqn{2^{-1/\theta}} \cr
\code{14} \tab \eqn{2-2^{1/\theta}} \tab - \cr
\code{16} \tab \eqn{2-2^{1/\theta}} \tab - \cr
\code{17} \tab \eqn{2-2^{1/\delta}} \tab \eqn{2^{-1/(\theta\delta)}} \cr
\code{18} \tab \eqn{2-2^{1/(\theta\delta)}} \tab - \cr
\code{19} \tab \eqn{2-2^{1/\theta}} \tab \eqn{2^{-1/\delta}} \cr
\code{20} \tab \eqn{2-2^{1/\theta}} if \eqn{\delta=1} otherwise 0 \tab - \cr
\code{23, 33} \tab - \tab - \cr
\code{24, 34} \tab - \tab - \cr
\code{26, 36} \tab - \tab - \cr
\code{27, 37} \tab - \tab - \cr
\code{28, 38} \tab - \tab - \cr
\code{29, 39} \tab - \tab - \cr
\code{30, 40} \tab - \tab - \cr
\code{104,204} \tab - \tab \eqn{\delta+1-(\delta^{\theta}+1)^{1/\theta}} \cr
\code{114, 214} \tab \eqn{1+\delta-(\delta^{\theta}+1)^{1/\theta}} \tab - \cr
\code{124, 224} \tab - \tab - \cr
\code{134, 234} \tab - \tab - \cr
}
}

\author{Eike Brechmann}

\references{
Joe, H. (1997).
Multivariate Models and Dependence Concepts.
Chapman and Hall, London.
}

\seealso{\code{\link{BiCopPar2Tau}}}

\examples{
## Example 1: Gaussian copula
BiCopPar2TailDep(1,0.7)

## Example 2: t copula
BiCopPar2TailDep(2,0.7,4)
}

