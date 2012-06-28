RVineMatrix  = function(Matrix,family=array(0,dim=dim(Matrix)),par=array(NA,dim=dim(Matrix)),par2=array(NA,dim=dim(Matrix)),names=NULL)
{

	Matrix[is.na(Matrix)]=0
	family[is.na(family)]=0
	family[upper.tri(family,diag=T)]=0
	par[is.na(par)]=0
	par[upper.tri(par,diag=T)]=0
	par2[is.na(par2)]=0
	par2[upper.tri(par2,diag=T)]=0
	
	if(dim(Matrix)[1]!=dim(Matrix)[2]) stop("Structure matrix has to be quadratic.")
	if(any(par!=NA) & dim(par)[1]!=dim(par)[2]) stop("Parameter matrix has to be quadratic.")
	if(any(par2!=NA) & dim(par2)[1]!=dim(par2)[2]) stop("Second parameter matrix has to be quadratic.")
	if(any(family!=0) & dim(family)[1]!=dim(family)[2]) stop("Copula family matrix has to be quadratic.")
	if(max(Matrix)>dim(Matrix)[1]) stop("Error in the structure matrix.")
	if(any(!(family %in% c(0,1:10,13,14,16:20,23,24,26:30,33,34,36:40,43,44)))) stop("Copula family not implemented.")
	if(length(names)>0 & length(names)!=dim(Matrix)[1]) stop("Length of the vector 'names' is not correct.")
	
	if(!all(par %in% c(0,NA)))
	{
	for(i in 2:dim(Matrix)[1]){
	for(j in 1:(i-1)){
    	if((family[i,j]==1 || family[i,j]==2) && abs(par[i,j])>=1) stop("The parameter of the Gaussian and t-copula has to be in the interval (-1,1).")
	if(family[i,j]==2 && par2[i,j]<=2) stop("The degrees of freedom parameter of the t-copula has to be larger than 2.")    	
	if((family[i,j]==3 || family[i,j]==13) && par[i,j]<=0) stop("The parameter of the Clayton copula has to be positive.")
    	if((family[i,j]==4 || family[i,j]==14) && par[i,j]<1) stop("The parameter of the Gumbel copula has to be in the interval [1,oo).")
    	if((family[i,j]==6 || family[i,j]==16) && par[i,j]<=1) stop("The parameter of the Joe copula has to be in the interval (1,oo).")	
    	if(family[i,j]==5 && par[i,j]==0) stop("The parameter of the Frank copula has to be unequal to 0.")
    	if((family[i,j]==7 || family[i,j]==17) && par[i,j]<=0) stop("The first parameter of the BB1 copula has to be positive.")
    	if((family[i,j]==7 || family[i,j]==17) && par2[i,j]<1) stop("The second parameter of the BB1 copula has to be in the interval [1,oo).")
	if((family[i,j]==8 || family[i,j]==18) && par[i,j]<=0) stop("The first parameter of the BB6 copula has to be in the interval [1,oo).")
	if((family[i,j]==8 || family[i,j]==18) && par2[i,j]<1) stop("The second parameter of the BB6 copula has to be in the interval [1,oo).")
    	if((family[i,j]==9 || family[i,j]==19) && par[i,j]<1) stop("The first parameter of the BB7 copula has to be in the interval [1,oo).")
    	if((family[i,j]==9 || family[i,j]==19) && par2[i,j]<=0) stop("The second parameter of the BB7 copula has to be positive.")
	if((family[i,j]==10 || family[i,j]==20) && par[i,j]<1) stop("The first parameter of the BB8 copula has to be in the interval [1,oo).")
	if((family[i,j]==10 || family[i,j]==20) && (par2[i,j]<=0 || par2[i,j]>1)) stop("The second parameter of the BB8 copula has to be in the interval (0,1].")    	
        if((family[i,j]==23 || family[i,j]==33) && par[i,j]>=0) stop("The parameter of the rotated Clayton copula has to be negative.")
    	if((family[i,j]==24 || family[i,j]==34) && par[i,j]>-1) stop("The parameter of the rotated Gumbel copula has to be in the interval (-oo,-1].")
    	if((family[i,j]==26 || family[i,j]==36) && par[i,j]>=-1) stop("The parameter of the rotated Joe copula has to be in the interval (-oo,-1).")
	if((family[i,j]==27 || family[i,j]==37) && par[i,j]>=0) stop("The first parameter of the rotated BB1 copula has to be negative.")
	if((family[i,j]==27 || family[i,j]==37) && par2[i,j]>-1) stop("The second parameter of the rotated BB1 copula has to be in the interval (-oo,-1].")
	if((family[i,j]==28 || family[i,j]==38) && par[i,j]>=0) stop("The first parameter of the rotated BB6 copula has to be in the interval (-oo,-1].")
	if((family[i,j]==28 || family[i,j]==38) && par2[i,j]>-1) stop("The second parameter of the rotated BB6 copula has to be in the interval (-oo,-1].")
	if((family[i,j]==29 || family[i,j]==39) && par[i,j]>-1) stop("The first parameter of the rotated BB7 copula has to be in the interval (-oo,-1].")
	if((family[i,j]==29 || family[i,j]==39) && par2[i,j]>=0) stop("The second parameter of the rotated BB7 copula has to be negative.")
	if((family[i,j]==30 || family[i,j]==40) && par[i,j]>-1) stop("The first parameter of the rotated BB8 copula has to be in the interval (-oo,-1].")
	if((family[i,j]==30 || family[i,j]==40) && (par2[i,j]>=0 || par2[i,j]<(-1))) stop("The second parameter of the rotated BB8 copula has to be in the interval [-1,0).")
  	  }
	  }
	}

  MaxMat=createMaxMat(Matrix)
  CondDistr=neededCondDistr(Matrix)
					
	RVM = list(Matrix = Matrix, family = family, par = par, par2 = par2, names = names, MaxMat = MaxMat, CondDistr = CondDistr)
			   
	class(RVM) = "RVineMatrix"
	return(RVM)	
}

normalizeRVineMatrix = function(RVM){
	
	oldOrder = diag(RVM$Matrix)
	Matrix = reorderRVineMatrix(RVM$Matrix)
	
	return(RVineMatrix(Matrix, RVM$family, RVM$par, RVM$par2, names = rev(RVM$names[oldOrder])))
}

reorderRVineMatrix = function(Matrix){
	oldOrder = diag(Matrix)
	
	O = apply(t(1:nrow(Matrix)),2,"==", Matrix)
	
	for(i in 1:nrow(Matrix)){
		Matrix[O[,oldOrder[i]]] = nrow(Matrix)-i+1
	}
	
	return(Matrix)
}

dim.RVineMatrix = function(x){
	RVine=x
	return(dim(RVine$Matrix)[1])
	NextMethod("dim")
}

print.RVineMatrix = function(x, ...){
	RVine=x
	message("R-vine matrix:")
	print(RVine$Matrix, ...)
	
	#Falls namen diese auch ausgeben
	if(!is.null(RVine$names)){
		message("")
		message("Where")
		for(i in 1:length(RVine$names)){
			message(i," <-> ",RVine$names[[i]])
		}
	}
	#NextMethod("print")
}

createMaxMat = function(Matrix){

  if(dim(Matrix)[1]!=dim(Matrix)[2]) stop("Structure matrix has to be quadratic.")
	
	MaxMat = reorderRVineMatrix(Matrix)
	
	n = nrow(MaxMat)
	
	for(j in 1:(n-1)){
		for(i in (n-1):j){
			MaxMat[i,j] = max(MaxMat[i:(i+1),j])
		}
	}
	
	tMaxMat = MaxMat
	tMaxMat[is.na(tMaxMat)] = 0
	
	oldSort = diag(Matrix)
	oldSort = oldSort[n:1]
	
	for(i in 1:n){
		MaxMat[tMaxMat == i] = oldSort[i]
	}
	
	return(MaxMat)
}

neededCondDistr = function(Vine){
	
	if(dim(Vine)[1]!=dim(Vine)[2]) stop("Structure matrix has to be quadratic.")
	
	Vine = reorderRVineMatrix(Vine)
	
	MaxMat = createMaxMat(Vine)
		
	d = nrow(Vine)
	
	M = list()
	M$direct = matrix(FALSE,d,d)
	M$indirect = matrix(FALSE,d,d)
	
	M$direct[2:d,1] = TRUE
	
	for(i in 2:(d-1)){
		v = d-i+1
		
		bw = as.matrix(MaxMat[i:d,1:(i-1)]) == v
		
		direct = Vine[i:d,1:(i-1)] == v
		
		M$indirect[i:d,i] = apply(as.matrix(bw & (!direct)),1,any)

		M$direct[i:d,i] = TRUE

		M$direct[i,i] = any(as.matrix(bw)[1,] & as.matrix(direct)[1,])
	}
	
	return(M)
}

as.RVineMatrix = function(RVine){

	n = length(RVine$Tree)+1
	con = list()
	names = V(RVine$Tree[[1]])$name
	
	conditionedSets = NULL
	corresppondingParams = list()
	corresppondingTypes = list()
	
	conditionedSets[[n-1]][[1]] = (E(RVine$Tree[[n-1]])$conditionedSet)
	for(k in 1:(n-2)){
		conditionedSets[[k]] = E(RVine$Tree[[k]])$conditionedSet
		corresppondingParams[[k]] = as.list(E(RVine$Tree[[k]])$Copula.param)
		corresppondingTypes[[k]] = as.list(E(RVine$Tree[[k]])$Copula.type)
	}
	corresppondingParams[[n-1]] = list()
	corresppondingParams[[n-1]][[1]] = (E(RVine$Tree[[n-1]])$Copula.param)
	corresppondingTypes[[n-1]] = as.list(E(RVine$Tree[[n-1]])$Copula.type)
	
	Param = array(dim=c(n,n))
	Params2 = array(0,dim=c(n,n))
	Type = array(dim=c(n,n))
	M = matrix(NA,n,n)
	
	for(k in 1:(n-1)){
		w = conditionedSets[[n-k]][[1]][1]
		
		M[k,k] = w
		M[(k+1),k] = conditionedSets[[n-k]][[1]][2]
		
		Param[(k+1),k] = corresppondingParams[[n-k]][[1]][1]
		Params2[(k+1),k] = corresppondingParams[[n-k]][[1]][2]
		
		Type[(k+1),k] = corresppondingTypes[[n-k]][[1]]
		
		if(k == (n-1)){
			M[(k+1),(k+1)] = conditionedSets[[n-k]][[1]][2]
		}else{
			for(i in (k+2):n){
				for(j in 1:length(conditionedSets[[n-i+1]])){
					cs = conditionedSets[[n-i+1]][[j]]
					if(cs[1] == w){
						M[i,k] = cs[2]
						break
					} else if(cs[2] == w){
						M[i,k] = cs[1]
						break
					}
				}
				Param[i,k] = corresppondingParams[[n-i+1]][[j]][1]
				Params2[i,k] = corresppondingParams[[n-i+1]][[j]][2]
				Type[i,k] = corresppondingTypes[[n-i+1]][[j]]
				
				conditionedSets[[n-i+1]][[j]] = NULL
				corresppondingParams[[n-i+1]][[j]] = NULL
				corresppondingTypes[[n-i+1]][[j]] = NULL
			}
		}
		
	}
	
	M = M+1 
	M[is.na(M)]=0
	Type[is.na(Type)]=0
	
	return(RVineMatrix(M, family = Type, par = Param, par2 = Params2, names = names))
	
}
