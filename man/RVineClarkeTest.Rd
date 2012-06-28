\name{RVineClarkeTest}     
\alias{RVineClarkeTest}                

\title{Clarke test comparing two R-vine copula models}

\description{
This function performs a Clarke test between two d-dimensional R-vine copula models as specified by their \code{\link{RVineMatrix}} objects.
}

\usage{
RVineClarkeTest(data, RVM1, RVM2)
}

\arguments{
  \item{data}{An N x d data matrix (with uniform margins).} 
  \item{RVM1, RVM2}{\code{\link{RVineMatrix}} objects of models 1 and 2.}
}

\value{
  \item{statistic, statistic.Akaike, statistic.Schwarz}{Test statistics without correction, with Akaike correction and with Schwarz correction.}
  \item{p.value, p.value.Akaike, p.value.Schwarz}{P-values of tests without correction, with Akaike correction and with Schwarz correction.}
}

\details{
The test proposed by Clarke (2007) allows to compare non-nested models.
For this let \eqn{c_1} and \eqn{c_2} be two competing vine copulas in terms of their densities
and with estimated parameter sets \eqn{\hat{\boldsymbol{\theta}}_1}{\theta_1} and \eqn{\hat{\boldsymbol{\theta}}_2}{\theta_2}.
The null hypothesis of statistical indistinguishability of the two models is 
\deqn{
H_0: P(m_i > 0) = 0.5\ \forall i=1,..,N,
}{
H_0: P(m_i > 0) = 0.5 forall i=1,..,N,
}
where \eqn{m_i:=\log\left[\frac{c_1(\boldsymbol{u}_i|\hat{\boldsymbol{\theta}}_1)}{c_2(\boldsymbol{u}_i|\hat{\boldsymbol{\theta}}_2)}\right]}{m_i:=log[ c_1(u_i|\theta_1) / c_2(u_i|\theta_2) ]}
for observations \eqn{\boldsymbol{u}_i,\ i=1,...,N}{u_i, i=1,...,N}.

Since under statistical equivalence of the two models the log likelihood ratios of the single observations
are uniformly distributed around zero and in expectation \eqn{50\%} of the log likelihood ratios greater than zero,
the tets statistic
\deqn{
\texttt{statistic} := B = \sum_{i=1}^N \mathbf{1}_{(0,\infty)}(m_i),
}{
statistic := B = \sum_{i=1}^N 1_{(0,\infty)}(m_i),
}
where \eqn{\mathbf{1}}{1} is the indicator function, is distributed Binomial with parameters \eqn{N} and \eqn{p=0.5},
and critical values can easily be obtained.
Model 1 is interpreted as statistically equivalent to model 2 if \eqn{B} is not significantly different from the expected value \eqn{Np = \frac{N}{2}}{np=N/2}.

Like AIC and BIC, the Clarke test statistic may be corrected for the number of parameters used in the models.
There are two possible corrections; 
the Akaike and the Schwarz corrections, which correspond to the penalty terms in the AIC and the BIC, respectively. 
}

\references{
Clarke, K. A. (2007).
A Simple Distribution-Free Test for Nonnested Model Selection.
Political Analysis, 15, 347-363.
}

\author{Jeffrey Dissmann, Eike Brechmann}

\seealso{\code{\link{RVineVuongTest}}, \code{\link{RVineAIC}}, \code{\link{RVineBIC}}}


\examples{
# define first 5-dimensional RVineMatrix object
Matrix1 = c(5,2,3,1,4,0,2,3,4,1,0,0,3,4,1,0,0,0,4,1,0,0,0,0,1)
Matrix1 = matrix(Matrix1,5,5)

family1 = c(0,1,3,4,4,0,0,3,4,1,0,0,0,4,1,0,0,0,0,3,0,0,0,0,0)
family1 = matrix(family1,5,5)

par1 = c(0,0.2,0.9,1.5,3.9,0,0,1.1,1.6,0.9,0,0,0,1.9,0.5,
         0,0,0,0,4.8,0,0,0,0,0)
par1 = matrix(par1,5,5)

RVM1 = RVineMatrix(Matrix=Matrix1,family=family1,par=par1,
                   par2=matrix(0,5,5),names=c("V1","V2","V3","V4","V5"))
                  
# define second 5-dimensional RVineMatrix object
Matrix2 = c(5,4,3,2,1,0,4,3,2,1,0,0,3,2,1,0,0,0,2,1,0,0,0,0,1)
Matrix2 = matrix(Matrix2,5,5)

family2 = c(0,3,1,3,2,0,0,1,5,3,0,0,0,2,3,0,0,0,0,1,0,0,0,0,0)
family2 = matrix(family2,5,5)

par2 = c(0,0.8,0.3,1.5,0.8,0,0,-0.4,1.6,1.2,0,0,0,-0.4,1.5,
         0,0,0,0,0.6,0,0,0,0,0)
par2 = matrix(par2,5,5)

nu2 = c(0,0,0,0,5,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0)
nu2 = matrix(nu2,5,5)

RVM2 = RVineMatrix(Matrix=Matrix2,family=family2,par=par2,par2=nu2,
                   names=c("V1","V2","V3","V4","V5"))                  

# simulate a sample of size 300 from the first R-vine copula model
simdata = RVineSim(300,RVM1)

# compare the two models based on this sample
clarke = RVineClarkeTest(simdata,RVM1,RVM2)
clarke$statistic
clarke$statistic.Schwarz
clarke$p.value
clarke$p.value.Schwarz
}

