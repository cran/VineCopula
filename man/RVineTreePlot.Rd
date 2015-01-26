\name{RVineTreePlot}      
\alias{RVineTreePlot}

\title{Visualisation of R-Vine Tree Structure}

\description{
This function plots one or all trees of a given R-vine copula model.
}

\usage{
RVineTreePlot(data = NULL, RVM, method = "mle", max.df = 30,
              max.BB = list(BB1=c(5,6),BB6=c(6,6),BB7=c(5,6),BB8=c(6,1)),
              tree = "ALL", edge.labels = c("family"), P = NULL)
}

\arguments{
  \item{data}{An N x d data matrix (with uniform margins), default: \code{data = NULL}.}
  \item{RVM}{An \code{\link{RVineMatrix}} object including the structure and the pair-copula families and parameters.}  
  \item{method}{Character indicating the estimation method:
    either maximum likelihood estimation (\code{method = "mle"}; default) or inversion of Kendall's tau (\code{method = "itau"}).}  
  \item{max.df}{Numeric; upper bound for the estimation of the degrees of freedom parameter of the t-copula
    (default: \code{max.df = 30}; for more details see \code{\link{BiCopEst}}).}
  \item{max.BB}{List; upper bounds for the estimation of the two parameters (in absolute values) of the BB1, BB6, BB7 and BB8 copulas \cr
    (default: \code{max.BB = list(BB1=c(5,6),BB6=c(6,6),BB7=c(5,6),BB8=c(6,1))}).}
  \item{tree}{Number of the tree to be plotted or \code{tree = "ALL"} (default) to plot all trees.}
  \item{edge.labels}{Vector of edge labels. Possible choices:\cr
    \code{FALSE}: no edge labels\cr
    \code{"family"}: pair-copula families (default)\cr
    \code{"par"}: pair-copula parameters\cr
    \code{"par2"}: second pair-copula parameters\cr
    \code{"theotau"}: theoretical Kendall's tau values corresponding to pair-copula families and parameters (see \code{\link{BiCopPar2Tau}})\cr
    \code{"emptau"}: empirical Kendall's tau values (only if data is provided!)\cr
	  \code{"pair"}: indices of (conditioned) pair of variables identified by the edges}
  \item{P}{A list of matrices with two columns for the x-y-coordinates of the nodes in the plot(s) (optional; default: \code{P = NULL}).}
}

\note{
The function computes the positions of the nodes automatically with the Fruchterman-Reingold algorithm (see \code{\link{plot.igraph}} 
for a detailed description).
If one would like to set the positions manually, one has to specify a list of matrices \code{P} in the argument list.
A good starting point may be to run the function \code{\link{RVineTreePlot}} and manipulate the returning matrix P.

If data is provided, the parameters of the R-vine copula model are estimated sequentially using \cr
\code{\link{RVineSeqEst}}/\code{\link{BiCopEst}}.
Then the edge width is chosen according to the empirical Kendall's tau values. Otherwise theoretical values are used.
}

\value{
A list of matrices \code{P} with two columns for the x-y-coordinates of the nodes in the plot(s).
}

\author{Eike Brechmann}

\seealso{\code{\link{BiCopName}}}

\examples{
# define 5-dimensional R-vine tree structure matrix
Matrix <- c(5, 2, 3, 1, 4,
            0, 2, 3, 4, 1,
            0, 0, 3, 4, 1,
            0, 0, 0, 4, 1,
            0, 0, 0, 0, 1)
Matrix <- matrix(Matrix, 5, 5)

# define R-vine pair-copula family matrix
family <- c(0, 1, 3, 4, 4,
            0, 0, 3, 4, 1,
            0, 0, 0, 4, 1,
            0, 0, 0, 0, 3,
            0, 0, 0, 0, 0)
family <- matrix(family, 5, 5)

# define R-vine pair-copula parameter matrix
par <- c(0, 0.2, 0.9, 1.5, 3.9,
         0, 0, 1.1, 1.6, 0.9,
         0, 0, 0, 1.9, 0.5,
         0, 0, 0, 0, 4.8,
         0, 0, 0, 0, 0)
par <- matrix(par, 5, 5)

# define second R-vine pair-copula parameter matrix
par2 <- matrix(0, 5, 5)

# define RVineMatrix object
RVM <- RVineMatrix(Matrix = Matrix, family = family,
                   par = par, par2 = par2,
                   names = c("V1", "V2", "V3", "V4", "V5"))

# set random seed for testing
set.seed(123)

# plot all trees with pair-copula families and
# theoretical Kendall's tau values as edge labels
P <- RVineTreePlot(data = NULL, RVM = RVM, tree = "ALL",
                   edge.labels = c("family","theotau"),
                   P = NULL)

# manipulate the first matrix of x-y-coordinates
P[[1]][1,] = P[[1]][1,]*2  

# re-set random seed for testing
set.seed(123)
# plot only the first tree with new coordinates
P <- RVineTreePlot(data = NULL, RVM = RVM, tree = 1,
                   edge.labels = FALSE, P = P)
}

