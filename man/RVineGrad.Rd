\name{RVineGrad}        
\alias{RVineGrad}                

\title{Gradient of the log-likelihood of an R-vine copula model}

\description{
This function calculates the gradient of the log-likelihood of a d-dimensional R-vine copula model with respect to the copula parameter and evaluates it on a given copula data set.
}

\usage{
RVineGrad(data, RVM, par=RVM$par, par2=RVM$par2, start.V=NA,
          posParams=(RVM$family>0))
}

\arguments{
  \item{data}{An N x d data matrix (with uniform margins).}
  \item{RVM}{An \code{\link{RVineMatrix}} object including the structure and the pair-copula families and parameters.}
  \item{par}{A d x d matrix with the pair-copula parameters (optional; default: \code{par = RVM$par}).}
  \item{par2}{A d x d matrix with the second parameters of pair-copula families with two parameters (optional; default: \code{par2 = RVM$par2}).}
  \item{start.V}{Transformations (h-functions and log-likelihoods of each pair-copula) of previous calculations (see output; default: \code{start.V = NA}).}
  \item{posParams}{A d x d matrix indicating which copula has to be considered in the gradient (default: \code{posParams = (RVM$family > 0)}).}
}

\note{
The gradient for R-vine copula models with two parameter Archimedean copulas, i.e. BB1, BB6, BB7, BB8 and their rotated versions.
}

\details{
The ordering of the gradient is due to the ordering of the R-vine matrix. The gradient starts at the lower right corner of the R-vine matrix and goes column by column 
to the left and up, i.e. the first entry of the gradient is the last entry of the second last column of the \code{par}-matrix followed by the last entry of the third last
column and the second last entry of this column.
If there is a copula family with two parameters, i.e. the t-copula, the derivative with respect to the second parameter is at the end of the gradient vector in order of their 
occurrence.
}

\value{
  \item{gradient}{The calculated gradient of the log-likelihood value of the R-vine copula model.}
  %\item{V}{The stored transformations (h-functions and log-likelihoods of each pair-copula) which may be used for posterior updates
  %  (three matrices: \code{direct}, \code{indirect} and \code{value}).}
}


\references{
Schepsmeier, U. and J. Stoeber (2012).
Derivatives and Fisher information of bivariate copulas.
Submitted for publication.
\url{http://mediatum.ub.tum.de/node?id=1106541}.

Stoeber, J. and U. Schepsmeier (2012).
Is there significant time-variation in multivariate dependence?.
In preparation.
\url{http://de.arxiv.org/abs/1205.4841}.

Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka (2011).
Selecting and estimating regular vine copulae and application to financial returns.
Submitted for publication.
\url{http://mediatum.ub.tum.de/node?id=1079277}
}

\author{Ulf Schepsmeier, Jakob Stoeber}

\seealso{\code{\link{BiCopDeriv}}, \code{\link{BiCopDeriv2}}, \code{\link{BiCopHfuncDeriv}}, \code{\link{BiCopHfuncDeriv2}}, \cr
\code{\link{RVineMatrix}}, \code{\link{RVineMLE}}, \code{\link{RVineHessian}}}

\examples{
# define 5-dimensional R-vine tree structure matrix
Matrix = c(5,2,3,1,4,0,2,3,4,1,0,0,3,4,1,0,0,0,4,1,0,0,0,0,1)
Matrix = matrix(Matrix,5,5)

# define R-vine pair-copula family matrix
family = c(0,1,3,4,4,0,0,3,4,1,0,0,0,4,1,0,0,0,0,3,0,0,0,0,0)
family = matrix(family,5,5)

# define R-vine pair-copula parameter matrix
par = c(0,0.2,0.9,1.5,3.9,0,0,1.1,1.6,0.9,0,0,0,1.9,0.5,
        0,0,0,0,4.8,0,0,0,0,0)
par = matrix(par,5,5)

# define second R-vine pair-copula parameter matrix
par2 = matrix(0,5,5)

# define RVineMatrix object
RVM = RVineMatrix(Matrix=Matrix,family=family,par=par,par2=par2,
                  names=c("V1","V2","V3","V4","V5"))

# simulate a sample of size 300 from the R-vine copula model
simdata = RVineSim(300,RVM)

# compute the gradient of the first row of the data
out2 = RVineGrad(simdata[1,],RVM)
out2$gradient

}
