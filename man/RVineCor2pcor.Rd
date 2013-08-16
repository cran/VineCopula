\name{RVineCor2pcor}
\Rdversion{1.1}
\alias{RVinePcor2cor}
\alias{RVineCor2pcor}

\title{correlations to partial correlations and vice versa for R-vines}
\description{
correlations to partial correlations and vice versa for R-vines
(C vine, D vine or general R vine) 
}

\usage{
RVineCor2pcor(RVM, corMat)
RVinePcor2cor(RVM)
}
\arguments{
  \item{RVM}{\code{\link{RVineMatrix}} defining only the R-vine structure for \code{cor2pcor} and providing as well the partial correlations for \code{pcor2cor}.}
  \item{corMat}{correlation matrix}
}
\value{
  \item{RVM}{RVineMatrix with transformed partial correlations}
  \item{cor}{correlation matrix}
}

\examples{
corMat <- matrix(c(1.00, 0.17, 0.15, 0.14, 0.13,
                   0.17, 1.00, 0.30, 0.28, 0.05,
                   0.15, 0.30, 1.00, 0.17, 0.05,
                   0.14, 0.28, 0.17, 1.00, 0.04,
                   0.13, 0.05, 0.05, 0.04, 1.00),5,5)

Matrix = matrix(c(5,2,3,1,4,
                  0,2,3,4,1,
                  0,0,3,4,1,
                  0,0,0,4,1,
                  0,0,0,0,1),5,5)
family = matrix(1,5,5)

par = matrix(c(0,0.2,0.9,0.5,0.8,
               0,  0,0.1,0.6,0.9,
               0,  0,  0,0.7,0.5,
               0,  0,  0,  0,0.8,
               0,  0,  0,  0,  0),5,5)

# define RVineMatrix object
RVM = RVineMatrix(Matrix,family,par)

# adjust the un-ordered RVine
newRVM <- RVineCor2pcor(RVM, corMat)
round(cor(qnorm(RVineSim(1000, newRVM)))-corMat, 2)

# normalise the RVine
normRVM <- RVineMatrixNormalize(RVM)

# adjust the normalised RVine
newNormRVM <- RVineCor2pcor(normRVM, corMat)

# newRVM and newNormRVM are the same vine using only different naming:
newNormRVM$par - newRVM$par

# the variable now do have a different ordering in the correlation matrix
newNormCor <- cor(qnorm(RVineSim(1000, newNormRVM)))
round(newNormCor,2)

# permuted, they meet the initial correlation matrix up to +/- 0.01
round(newNormCor[c(1,4,3,2,5),c(1,4,3,2,5)]-corMat, 2)

# re-order names of the normalised RVine generating a new RVine
normRVM2 <- normRVM
normRVM2$names <- c("V1", "V2", "V3", "V4", "V5")

# adjust the normalised RVine
newNormRVM2 <- RVineCor2pcor(normRVM2, corMat)
# check whether the parameters are different beyond permutation (that's why 
# permutation does not work)
newNormRVM2$par
newRVM$par

# adjust the normalised RVine
newNormRVM2 <- RVineCor2pcor(normRVM2, corMat[c(1,4,3,2,5),c(1,4,3,2,5)])
# check whether the parameters are now identical
round(newNormRVM2$par - newRVM$par,2)

# back and forth
RVinePcor2cor(RVineCor2pcor(RVM, corMat))-corMat
RVinePcor2cor(RVineCor2pcor(normRVM, corMat))-corMat
RVinePcor2cor(RVineCor2pcor(normRVM2, corMat))-corMat
}

\keyword{vine}
\keyword{partial correlation}
