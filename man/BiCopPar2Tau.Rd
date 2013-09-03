\name{BiCopPar2Tau}      
\alias{BiCopPar2Tau} 

\title{Kendall's tau value of a bivariate copula}

\description{
This function computes the theoretical Kendall's tau value of a bivariate copula for given parameter values. 
}

\usage{
BiCopPar2Tau(family, par, par2=0)
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
  \item{par2}{Second parameter for the two parameter t-, BB1, BB6, BB7, BB8, Tawn type 1 and type 2 copulas (default: \code{par2 = 0}).
    Note that the degrees of freedom parameter of the t-copula does not need to be set,
    because the theoretical Kendall's tau value of the t-copula is independent of this choice.}
}                                                                                  

\value{
Theoretical value of Kendall's tau corresponding to the bivariate copula family and parameter(s)
(\eqn{\theta} for one parameter families and the first parameter of the t-copula,
\eqn{\theta} and \eqn{\delta} for the two parameter BB1, BB6, BB7 and BB8 copulas). 
\tabular{ll}{
No. \tab Kendall's tau \cr
\code{1, 2} \tab \eqn{\frac{2}{\pi}\arcsin(\theta)}{2 / \pi arcsin(\theta)} \cr
\code{3, 13} \tab \eqn{\frac{\theta}{\theta+2}}{\theta / (\theta+2)} \cr
\code{4, 14} \tab \eqn{1-\frac{1}{\theta}}{1-1/\theta} \cr
\code{5} \tab \eqn{1-\frac{4}{\theta}+4\frac{D_1(\theta)}{\theta}}{1-4/\theta + 4 D_1(\theta)/\theta} \cr
	 \tab with \eqn{D_1(\theta)=\int_0^\theta \frac{x/\theta}{\exp(x)-1}dx}{D_1(\theta)=\int_0^\theta (x/\theta)/(exp(x)-1)dx} (Debye function) \cr
\code{6, 16} \tab \eqn{1+\frac{4}{\theta^2}\int_0^1 x\log(x)(1-x)^{2(1-\theta)/\theta}dx}{1+4/\theta^2\int_0^1 x\log(x)(1-x)^{2(1-\theta)/\theta}dx} \cr
\code{7, 17} \tab \eqn{1-\frac{2}{\delta(\theta+2)}}{1-2/(\delta(\theta+2))} \cr
\code{8, 18} \tab \eqn{1+4\int_0^1 -\log(-(1-t)^\theta+1)(1-t-(1-t)^{-\theta}+(1-t)^{-\theta}t)/(\delta\theta) dt} \cr
\code{9, 19} \tab \eqn{1+4\int_0^1 ( (1-(1-t)^{\theta})^{-\delta} - )/( -\theta\delta(1-t)^{\theta-1}(1-(1-t)^{\theta})^{-\delta-1} ) dt} \cr
\code{10, 20} \tab \eqn{1+4\int_0^1 -\log \left( ((1-t\delta)^\theta-1)/((1-\delta)^\theta-1) \right) } \cr
	 \tab \eqn{* (1-t\delta-(1-t\delta)^{-\theta}+(1-t\delta)^{-\theta}t\delta)/(\theta\delta) dt} \cr
\code{23, 33} \tab \eqn{\frac{\theta}{2-\theta}}{\theta/(2-\theta)} \cr
\code{24, 34} \tab \eqn{-1-\frac{1}{\theta}}{-1-1/\theta} \cr
\code{26, 36} \tab \eqn{-1-\frac{4}{\theta^2}\int_0^1 x\log(x)(1-x)^{-2(1+\theta)/\theta}dx}{-1-4/\theta^2\int_0^1 x\log(x)(1-x)^{-2(1+\theta)/\theta}dx} \cr
\code{27, 37} \tab \eqn{-1-\frac{2}{\delta(2-\theta)}}{1-2/(\delta(\theta+2))} \cr
\code{28, 38} \tab \eqn{-1-4\int_0^1 -\log(-(1-t)^{-\theta}+1)(1-t-(1-t)^{\theta}+(1-t)^{\theta}t)/(\delta\theta) dt} \cr
\code{29, 39} \tab \eqn{-1-4\int_0^1 ( (1-(1-t)^{-\theta})^{\delta} - )/( -\theta\delta(1-t)^{-\theta-1}(1-(1-t)^{-\theta})^{\delta-1} ) dt} \cr
\code{30, 40} \tab \eqn{-1-4\int_0^1 -\log \left( ((1+t\delta)^{-\theta}-1)/((1+\delta)^{-\theta}-1) \right)} \cr 
	 \tab \eqn{* (1+t\delta-(1+t\delta)^{\theta}-(1+t\delta)^{\theta}t\delta)/(\theta\delta) dt} \cr
   \code{104,114} \tab \eqn{\int_0^1 \frac{t(1-t)A^{\prime\prime}(t)}{A(t)}dt} \cr
   \tab with \eqn{A(t) = (1-\delta)t+[(\delta(1-t))^{\theta}+t^{\theta}]^{1/\theta}} \cr
   \code{204,214} \tab \eqn{\int_0^1 \frac{t(1-t)A^{\prime\prime}(t)}{A(t)}dt} \cr
   \tab with \eqn{A(t) = (1-\delta)(1-t)+[(1-t)^{-\theta}+(\delta t)^{-\theta}]^{-1/\theta}} \cr
   \code{124,134} \tab \eqn{-\int_0^1 \frac{t(1-t)A^{\prime\prime}(t)}{A(t)}dt} \cr
   \tab with \eqn{A(t) = (1-\delta)t+[(\delta(1-t))^{-\theta}+t^{-\theta}]^{-1/\theta}} \cr
   \code{224,234} \tab \eqn{-\int_0^1 \frac{t(1-t)A^{\prime\prime}(t)}{A(t)}dt} \cr
   \tab with \eqn{A(t) = (1-\delta)(1-t)+[(1-t)^{-\theta}+(\delta t)^{-\theta}]^{-1/\theta}} \cr
   
}
}

\author{Ulf Schepsmeier}

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
tt1 = BiCopPar2Tau(1,0.7)

# transform back
BiCopTau2Par(1,tt1)


## Example 2: Clayton copula
BiCopPar2Tau(3,1.3)
}

