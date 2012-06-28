\name{BiCopGofKendall}        
\alias{BiCopGofKendall}     

\title{Goodness-of-fit test based on Kendall's process for bivariate copula data}

\description{
This function performs the goodness-of-fit test based on Kendall's process for bivariate copula data.
It computes the Cramer-von Mises and Kolmogorov-Smirnov test statistics, respectively, as well as the according p-values using bootstrapping.
}

\usage{
BiCopGofKendall(u1, u2, family, B=100, level=0.05)
}

\arguments{
  \item{u1,u2}{Data vectors of equal length with values in [0,1].}
  \item{family}{An integer defining the bivariate copula family for which the test is performed: \cr
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
		\code{40} = rotated BB8 copula (270 degrees) 
		}
  \item{B}{Integer; number of bootstrap samples (default: \code{B = 100}).
    For \code{B = 0} only the the test statistics are returned.
    WARNING: If \code{B} is chosen too large, computations will take very long.}
  \item{level}{Numeric; significance level of the goodness-of-fit test (default: \code{level = 0.05}).}  
}

\value{
  \item{p.value.CvM}{P-value of the goodness-of-fit test using the Cramer-von Mises statistic \cr
  (if \code{B > 0}).}
  \item{p.value.KS}{P-value of the goodness-of-fit test using the Kolmogorov-Smirnov statistic \cr
  (if \code{B > 0}).}
  \item{statistic.CvM}{The observed Cramer-von Mises test statistic.}
  \item{statistic.KS}{The observed Kolmogorov-Smirnov test statistic.}
}

\details{
This copula goodness-of-fit test is based on Kendall's process as investigated by Genest and Rivest (1993) and Wang and Wells (2000).
For rotated copulas the input arguments are transformed and the goodness-of-fit procedure for the corresponding non-rotated copula is used. 
}

\references{
Genest, C. and L.-P. Rivest (1993).
Statistical inference procedures for bivariate Archimedean copulas.
Journal of the American Statistical Association, 88 (423), 1034-1043.

Luo J. (2011).
Stepwise estimation of D-vines with arbitrary specified copula pairs and EDA Tools.
Diploma thesis, Technische Universitaet Muenchen.\cr
\url{http://mediatum.ub.tum.de/?id=1079291}.

Wang, W. and M. T. Wells (2000).
Model selection and semiparametric inference for bivariate failure-time data.
Journal of the American Statistical Association, 95 (449), 62-72.
}

\author{Jiying Luo, Eike Brechmann}

\seealso{\code{\link{BiCopIndTest}}, \code{\link{BiCopSelect}}, \code{\link{BiCopVuongClarke}}, \code{\link{BiCopKPlot}}, \code{\link{BiCopLambda}}}

\examples{
# sample from a Gaussian copula
par1 = 3
fam1 = 3
dat1 = BiCopSim(500,fam1,par1)

# perform the goodness-of-fit test for the true copula
gof = BiCopGofKendall(dat1[,1],dat1[,2],fam1)
gof$p.value.CvM
gof$p.value.KS

# perform the goodness-of-fit test for the Frank copula
gof = BiCopGofKendall(dat1[,1],dat1[,2],5)
gof$p.value.CvM
gof$p.value.KS
}
