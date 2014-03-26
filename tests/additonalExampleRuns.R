## switch for testing the following time consuming examples
docheck <- FALSE

if(docheck){
  
## tests from excluded examples
library(VineCopula)

## Not run: 
# chi-plots for bivariate Gaussian copula data
n = 500
tau = 0.5

# simulate copula data
fam = 1  
theta = BiCopTau2Par(fam,tau)
set.seed(666)
dat = BiCopSim(n,fam,theta)	

# create chi-plots
par(mfrow=c(1,3))
BiCopChiPlot(dat[,1],dat[,2],xlim=c(-1,1),ylim=c(-1,1),
             main="General chi-plot")
BiCopChiPlot(dat[,1],dat[,2],mode="lower",xlim=c(-1,1),
             ylim=c(-1,1),main="Lower chi-plot")
BiCopChiPlot(dat[,1],dat[,2],mode="upper",xlim=c(-1,1),
             ylim=c(-1,1),main="Upper chi-plot")

# simulate from a bivariate Clayton copula
set.seed(666)
simdata = BiCopSim(300,3,2)
u1 = simdata[,1]
u2 = simdata[,2]

# perform White's goodness-of-fit test for the true copula
BiCopGofTest(u1,u2,family=3)

# perform Kendall's goodness-of-fit test for the Frank copula
BiCopGofTest(u1,u2,family=5)

# perform Kendall's goodness-of-fit test for the true copula
gof = BiCopGofTest(u1,u2,family=3,method="kendall")
gof$p.value.CvM
gof$p.value.KS

# perform Kendall's goodness-of-fit test for the Frank copula
gof = BiCopGofTest(u1,u2,family=5,method="kendall")
gof$p.value.CvM
gof$p.value.KS
  
# Not run: 
# Gaussian and Clayton copulas
n = 500
tau = 0.5

# simulate from Gaussian copula
fam1 = 1  
theta1 = BiCopTau2Par(fam1,tau)
set.seed(666)
dat1 = BiCopSim(n,fam1,theta1)	

# simulate from Clayton copula
fam2 = 3
theta2 = BiCopTau2Par(fam2,tau)
set.seed(666)
dat2 = BiCopSim(n,fam2,theta2)

# create K-plots
par(mfrow=c(1,2))
BiCopKPlot(dat1[,1],dat1[,2],main="Gaussian copula")
BiCopKPlot(dat2[,1],dat2[,2],main="Clayton copula")
# End(Not run)
  
# Not run: 
# Clayton and rotated Clayton copulas
n = 1000
tau = 0.5

# simulate from Clayton copula
fam = 3  
theta = BiCopTau2Par(fam,tau)
set.seed(666)
dat = BiCopSim(n,fam,theta)

# create lambda-function plots
par(mfrow=c(1,3))
BiCopLambda(dat[,1],dat[,2])	# empirical lambda-function	
BiCopLambda(family=fam,par=theta)	# theoretical lambda-function
BiCopLambda(dat[,1],dat[,2],family=fam,par=theta)	# both

# simulate from rotated Clayton copula (90 degrees)
fam = 23  
theta = BiCopTau2Par(fam,-tau)
set.seed(666)
dat = BiCopSim(n,fam,theta)

# rotate the data to standard Clayton copula data
rot_dat = 1-dat[,1]

par(mfrow=c(1,3))
BiCopLambda(rot_dat,dat[,2])  # empirical lambda-function	
BiCopLambda(family=3,par=-theta)	# theoretical lambda-function
BiCopLambda(rot_dat,dat[,2],family=3,par=-theta)	# both
# End(Not run)
  
# Not run: 
## Example 3: empirical data
data(daxreturns)
cop3 = BiCopSelect(daxreturns[,1],daxreturns[,4],
       familyset=c(1:10,13,14,16,23,24,26))
cop3$family
cop3$par
cop3$par2
# End(Not run)
  
# Not run: 
# simulate from a t-copula
set.seed(666)
dat = BiCopSim(500,2,0.7,5)

# apply the test for families 1-10
vcgof = BiCopVuongClarke(dat[,1],dat[,2],familyset=c(1:10))

# display the Vuong test scores
vcgof[1,]
# End(Not run)
  
# Not run: 
# select the R-vine structure, families and parameters
RVM = RVineStructureSelect(daxreturns[,1:5],c(1:6))
RVM$Matrix
RVM$par
RVM$par2

# select the C-vine structure, families and parameters
CVM = RVineStructureSelect(daxreturns[,1:5],c(1:6),type="CVine")
CVM$Matrix
CVM$par
CVM$par2

# compare the two models based on the data
clarke = RVineClarkeTest(daxreturns[,1:5],RVM,CVM)
clarke$statistic
clarke$statistic.Schwarz
clarke$p.value
clarke$p.value.Schwarz
# End(Not run)
  
# Not run: 
# White test with asymptotic p-value
RVineGofTest(daxreturns[,1:5], RVM, B=0)

# ECP2 test with Cramer-von-Mises test statistic and a bootstrap with 200 replications 
# for the calculation of the p-value
RVineGofTest(daxreturns[,1:5], RVM, method="ECP2", statistic="CvM", B=200)
# End(Not run)
  
  ## Not run: 
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
set.seed(666)
simdata = RVineSim(300,RVM)

# compute the MLE
mle = RVineMLE(simdata,RVM,grad=TRUE)
mle$RVM
# End(Not run)

##TODO shorten this test, takes too long
# # Not run: 
# RVM = RVineStructureSelect(daxreturns,c(1:6),progress=TRUE)
# # End(Not run)
# 
# # specify a C-vine copula model with only Clayton, Gumbel and Frank copulas
# # Not run: 
# CVM = RVineStructureSelect(daxreturns,c(3,4,5),"CVine")
# # End(Not run)
# # determine the order of the nodes in a D-vine using the package TSP
# # Not run: 
# library(TSP)
# d = dim(daxreturns)[2]
# M = 1 - abs(TauMatrix(daxreturns))
# hamilton = insert_dummy(TSP(M),label="cut")
# sol = solve_TSP(hamilton,method="repetitive_nn")
# order = cut_tour(sol,"cut")
# DVM = D2RVine(order,family=rep(0,d*(d-1)/2),par=rep(0,d*(d-1)/2))
# RVineCopSelect(daxreturns,c(1:6),DVM$Matrix)
# End(Not run)
  
# Not run: 
RVM = RVineStructureSelect(daxreturns[,1:5],c(1:6))
CVM = RVineStructureSelect(daxreturns[,1:5],c(1:6),type="CVine")

# compare the two models based on the data
vuong = RVineVuongTest(daxreturns[,1:5],RVM,CVM)
vuong$statistic
vuong$statistic.Schwarz
vuong$p.value
vuong$p.value.Schwarz
# End(Not run)
  
# Not run: 
library(copula)
vine <- vineCopula(4L,"CVine")

set.seed(666)
rCopula(500,vine)

# End(Not run)

}