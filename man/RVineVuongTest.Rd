\name{RVineVuongTest}         
\alias{RVineVuongTest}                

\title{Vuong test comparing two R-vine copula models}

\description{
This function performs a Vuong test between two d-dimensional R-vine copula models as specified by their \code{\link{RVineMatrix}} objects.
}

\usage{
RVineVuongTest(data, RVM1, RVM2)
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
The likelihood-ratio based test proposed by Vuong (1989) can be used for comparing non-nested models.
For this let \eqn{c_1} and \eqn{c_2} be two competing vine copulas in terms of their densities
and with estimated parameter sets \eqn{\hat{\boldsymbol{\theta}}_1}{\theta_1} and \eqn{\hat{\boldsymbol{\theta}}_2}{\theta_2}.
We then compute the standardized sum, \eqn{\nu}, of the log differences of their 
pointwise likelihoods \eqn{m_i:=\log\left[\frac{c_1(\boldsymbol{u}_i|\hat{\boldsymbol{\theta}}_1)}{c_2(\boldsymbol{u}_i|\hat{\boldsymbol{\theta}}_2)}\right]}{m_i:=log[c_1(u_i|\theta_1) / c_2(u_i|\theta_2) ]}
for observations \eqn{\boldsymbol{u}_i\in[0,1],\ i=1,...,N}{u_i in [0,1],i=1,...,N} , i.e.,
\deqn{
\texttt{statistic} := \nu = \frac{\frac{1}{n}\sum_{i=1}^N m_i}{\sqrt{\sum_{i=1}^N\left(m_i - \bar{m} \right)^2}}.
}{
statistic := \nu = (1/n\sum_{i=1}^N m_i) / ((\sum_{i=1}^N (m_i - \bar{m} )^2)^0.5).
}
Vuong (1989) shows that \eqn{\nu} is asymptotically standard normal.
According to the null-hypothesis
\deqn{
H_0: E[m_i] = 0\ \forall i=1,...,N,
}{
H_0: E[m_i] = 0 forall i=1,...,N,
}
we hence prefer vine model 1 to vine model 2 at level \eqn{\alpha} if
\deqn{
\nu>\Phi^{-1}\left(1-\frac{\alpha}{2}\right),
}{
\nu > \Phi^{-1}(1-\alpha/2),
}
where \eqn{\Phi^{-1}} denotes the inverse of the standard normal distribution function.
If \eqn{\nu<-\Phi^{-1}\left(1-\frac{\alpha}{2}\right)}{\nu<-\Phi^{-1}(1-\alpha/2)} we choose model 2. 
If, however, \eqn{|\nu|\leq\Phi^{-1}\left(1-\frac{\alpha}{2}\right)}{|\nu| <= \Phi^{-1}(1-\alpha/2)}, no decision among the models is possible.

Like AIC and BIC, the Vuong test statistic may be corrected for the number of parameters used in the models.
There are two possible corrections; 
the Akaike and the Schwarz corrections, which correspond to the penalty terms in the AIC and the BIC, respectively. 
}


\references{
Vuong, Q. H. (1989).
Ratio tests for model selection and non-nested hypotheses.
Econometrica 57 (2), 307-333.
}

\author{Jeffrey Dissmann, Eike Brechmann}

\seealso{\code{\link{RVineClarkeTest}}, \code{\link{RVineAIC}}, \code{\link{RVineBIC}}}


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
vuong = RVineVuongTest(simdata,RVM1,RVM2)
vuong$statistic
vuong$statistic.Schwarz
vuong$p.value
vuong$p.value.Schwarz
}

