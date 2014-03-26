RVineCopSelect <- function(data,familyset=NA,Matrix,selectioncrit="AIC",indeptest=FALSE,level=0.05,trunclevel=NA) {

	n = dim(data)[2]
	N = nrow(data)

	if(dim(Matrix)[1]!=dim(Matrix)[2]) 
    stop("Structure matrix has to be quadratic.")
	if(max(Matrix)>dim(Matrix)[1]) 
    stop("Error in the structure matrix.")

	if(N<2) 
    stop("Number of observations has to be at least 2.")
	if(n<2) 
    stop("Dimension has to be at least 2.")
	if(any(data>1) || any(data<0)) 
    stop("Data has be in the interval [0,1].")
  if(!is.na(familyset[1])) {
      if(any(!(familyset %in% c(0,1:10,13,14,16:20,23,24,26:30,33,34,36:40,104,114,124,134,204,214,224,234)))) 
        stop("Copula family not implemented.")
  }
  if(selectioncrit != "AIC" && selectioncrit != "BIC") 
    stop("Selection criterion not implemented.")
  if(level < 0 & level > 1) 
    stop("Significance level has to be between 0 and 1.")

  if(is.na(trunclevel)) 
    trunclevel = n

  types = familyset
	if(trunclevel == 0) 
    types = 0
  
	M = Matrix

  Mold = M

	o = diag(M)
	M = reorderRVineMatrix(M)

  data = data[,o[length(o):1]]

  MaxMat = createMaxMat(M)
  CondDistr = neededCondDistr(M)

  Types = matrix(0,n,n)

  Params = matrix(0,n,n)
  Params2 = matrix(0,n,n)

  V = list()
  V$direct = array(NA,dim=c(n,n,N))
  V$indirect = array(NA,dim=c(n,n,N))
  V$direct[n,,] = t(data[,n:1])

	for(i in (n-1):1){
		for(k in n:(i+1)){

			m = MaxMat[k,i]

			zr1 = V$direct[k,i,]

			if(m == M[k,i]){
				zr2 = V$direct[k,(n-m+1),]
			}else{
				zr2 = V$indirect[k,(n-m+1),]
			}

      if(n+1-k > trunclevel){
        outcop = BiCopSelect(zr2,zr1,0,selectioncrit,indeptest,level)
      }else{
        #outcop = BiCopSelect(zr1,zr2,types,selectioncrit,indeptest,level)
        outcop = BiCopSelect(zr2,zr1,types,selectioncrit,indeptest,level)
      }

      Types[k,i] = outcop$family
      Params[k,i] = outcop$par
      Params2[k,i] = outcop$par2

			if(CondDistr$direct[k-1,i])
        #V$direct[k-1,i,] = outcop$CondOn.2
				V$direct[k-1,i,] = .C("Hfunc1",as.integer(Types[k,i]),as.integer(N),as.double(zr1),as.double(zr2),as.double(Params[k,i]),as.double(Params2[k,i]),as.double(rep(0,N)),PACKAGE='VineCopula')[[7]]
			if(CondDistr$indirect[k-1,i])
        #V$indirect[k-1,i,] = outcop$CondOn.1
				V$indirect[k-1,i,] = .C("Hfunc2",as.integer(Types[k,i]),as.integer(N),as.double(zr2),as.double(zr1),as.double(Params[k,i]),as.double(Params2[k,i]),as.double(rep(0,N)),PACKAGE='VineCopula')[[7]]
		}
	}

  varnames = paste("V", 1:n, sep="")
 
 print(Types)
  
  RVM = RVineMatrix(Mold, family=Types, par=Params, par2=Params2, names=varnames)

	return(RVM)
}
