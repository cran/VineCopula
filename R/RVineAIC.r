RVineAIC <-function(data,RVM,par=RVM$par,par2=RVM$par2){
	
  if(is.vector(data)){
    data = t(as.matrix(data))
  }else{
    data=as.matrix(data)
  }
  if(any(data>1) || any(data<0)) stop("Data has be in the interval [0,1].")
	d=dim(data)[2]
	T=dim(data)[1]
	n<-d
	N<-T
	if(n != dim(RVM)) stop("Dimensions of 'data' and 'RVM' do not match.")
  if(!is(RVM,"RVineMatrix")) stop("'RVM' has to be an RVineMatrix object.")
  
	npar = sum(RVM$family >= 1, na.rm=TRUE) + sum(RVM$family %in% c(2,7:10,17:20,27:30,37:40),na.rm=TRUE)
  npar_pair = (RVM$family>=1)+(RVM$family%in%c(2,7:10,17:20,27:30,37:40))

  like = RVineLogLik(data,RVM)

  AIC = -2*like$loglik + 2*npar
  pair.AIC = -2*like$V$value + 2*npar_pair
  
  return(list(AIC=AIC,pair.AIC=pair.AIC))
}

RVineBIC <-function(data,RVM,par=RVM$par,par2=RVM$par2){
	
  if(is.vector(data)){
    data = t(as.matrix(data))
  }else{
    data=as.matrix(data)
  }
	d=dim(data)[2]
	T=dim(data)[1]
	n<-d
	N<-T
	if(n != dim(RVM)) stop("Dimensions of 'data' and 'RVM' do not match.")
  if(is(RVM) != "RVineMatrix") stop("'RVM' has to be an RVineMatrix object.")
  
	npar = sum(RVM$family >= 1, na.rm=TRUE) + sum(RVM$family %in% c(2,7:10,17:20,27:30,37:40),na.rm=TRUE)
  npar_pair = (RVM$family>=1)+(RVM$family%in%c(2,7:10,17:20,27:30,37:40))

  like = RVineLogLik(data,RVM)

  BIC = -2*like$loglik + log(T)*npar
  pair.BIC = -2*like$V$value + log(T)*npar_pair
  
  return(list(BIC=BIC,pair.BIC=pair.BIC))	
}

