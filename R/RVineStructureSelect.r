RVineStructureSelect = function(data,familyset=NA,type=0,selectioncrit="AIC",indeptest=FALSE,level=0.05,trunclevel=NA,progress=FALSE,weights=NA){

  if(type == 0) type = "RVine"
  else if(type == 1) type = "CVine"
  if(type != "RVine" & type != "CVine") stop("Vine model not implemented.")

	n = dim(data)[2]
	d = dim(data)[1]
	
	if(dim(data)[1]<2) stop("Number of observations has to be at least 2.")
  if(d<2) stop("Dimension has to be at least 2.")
  if(any(data>1) || any(data<0)) stop("Data has be in the interval [0,1].")
	
  if(!is.na(familyset[1])) for(i in 1:length(familyset)) if(!(familyset[i] %in% c(0,1:10,13,14,16:20,23,24,26:30,33,34,36:40, 104,114,124,134,204,214,224,234))) stop("Copula family not implemented.")  
  if(selectioncrit != "AIC" && selectioncrit != "BIC") stop("Selection criterion not implemented.")
  if(level < 0 & level > 1) stop("Significance level has to be between 0 and 1.")
  	
	if(is.null(colnames(data))) colnames(data) = paste("V",1:n,sep="") 

  if(is.na(trunclevel)) trunclevel = d

	RVine = list(Tree = NULL, Graph=NULL)

  if(trunclevel == 0) familyset = 0
	
	g = initializeFirstGraph(data,weights)
	mst = findMaximumTauTree(g,mode=type)
	VineTree = fit.FirstTreeCopulas(mst,data,familyset,selectioncrit,indeptest,level,weights=weights)
	
	RVine$Tree[[1]] = VineTree
	RVine$Graph[[1]] = g
	oldVineGraph  = VineTree
  
	
	for(i in 2:(n-1)){

    if(trunclevel == i-1) familyset = 0
	
		g = buildNextGraph(VineTree,weights)
		mst = findMaximumTauTree(g,mode=type)
		
		VineTree = fit.TreeCopulas(mst, VineTree,familyset,selectioncrit,indeptest,level,progress,weights=weights)
		
		RVine$Tree[[i]] = VineTree
		RVine$Graph[[i]] = g
	}
	
	return(as.RVM(RVine))
}

initializeFirstGraph <- function(data.univ,weights)
{

	#C = cor(data.univ,method="kendall")
	q=dim(data.univ)[2]
	C=matrix(rep(1,q*q), ncol=q)

	for(i in 1:(q-1))
	{
		for(j in (i+1):q)
		{
			tau=fasttau(data.univ[,i],data.univ[,j],weights)
			C[i,j]=tau
			C[j,i]=tau
		}
	}

	rownames(C)=colnames(C)=colnames(data.univ)

	g = graph.adjacency(C, mode="lower",weighted=TRUE,diag=FALSE)
	

	E(g)$tau = E(g)$weight
	
	E(g)$name = paste(get.edgelist(g)[,1],get.edgelist(g)[,2],sep=",")
	
	for(i in 1:ecount(g)){
		E(g)$conditionedSet[[i]] = get.edges(g,i)
	}
	return(g)
}

findMaximumTauTree <- function(g,mode="RVine")
{

	if(mode == "RVine")
	{
		return(minimum.spanning.tree(g, weights=1-abs(E(g)$weight)))
	}
	else if(mode == "CVine")
	{
		M = abs(get.adjacency(g,attr="weight",sparse=0))
		sumtaus = rowSums(M)
		root = which.max(sumtaus)
		
		Ecken = get.edges(g,1:ecount(g))
		pos = Ecken[,2]== root | Ecken[,1]== root
		
		mst = delete.edges(g, E(g)[!pos])
		
		return(mst)
	}
}

fit.FirstTreeCopulas <- function(mst,data.univ,type,copulaSelectionBy,testForIndependence,testForIndependence.level,weights=NA)
{
	
	d = ecount(mst)
	
	parameterForACopula = list()
	
	for(i in 1:d)
	{
		parameterForACopula[[i]] = list()
		
		a = get.edges(mst,i)
		
		parameterForACopula[[i]]$zr1 = data.univ[,a[1]]
		parameterForACopula[[i]]$zr2 = data.univ[,a[2]]
		
		E(mst)[i]$Copula.Data.1 =  list(data.univ[,a[1]])
		E(mst)[i]$Copula.Data.2 =  list(data.univ[,a[2]])
		
		if(is.null(V(mst)[a[1]]$name))
			E(mst)[i]$Copula.CondName.1 = a[1]
		else
			E(mst)[i]$Copula.CondName.1 = V(mst)[a[1]]$name
		
		if(is.null(V(mst)[a[2]]$name))
			E(mst)[i]$Copula.CondName.2 = a[2]
		else
			E(mst)[i]$Copula.CondName.2 = V(mst)[a[2]]$name
		
		if(is.null(V(mst)[a[1]]$name) || is.null(V(mst)[a[2]]$name))
			E(mst)[i]$Copula.Name = paste(a[1],a[2],sep=" , ")
		else
			E(mst)[i]$Copula.Name = paste(V(mst)[a[1]]$name,V(mst)[a[2]]$name,sep=" , ")
	}

	outForACopula = lapply(X = parameterForACopula, FUN=wrapper_fit.ACopula, type,copulaSelectionBy,testForIndependence,testForIndependence.level,weights)
	
	for(i in 1:d)
	{
		E(mst)$Copula.param[[i]] = c(outForACopula[[i]]$par,outForACopula[[i]]$par2)
		E(mst)[i]$Copula.type = outForACopula[[i]]$family
    E(mst)[i]$Copula.out = list(outForACopula[[i]])
		
		E(mst)[i]$Copula.CondData.1 <- list(outForACopula[[i]]$CondOn.1)
		E(mst)[i]$Copula.CondData.2 <- list(outForACopula[[i]]$CondOn.2)
	}
	
	return(mst)
}

fit.TreeCopulas <- function(mst, oldVineGraph, type,copulaSelectionBy,testForIndependence,testForIndependence.level,progress,weights=NA)
{
	d = ecount(mst)
	
	parameterForACopula = list()
	
	for(i in 1:d)
	{
		parameterForACopula[[i]] = list()
		
		con = get.edge(mst,i)
		
		temp = get.edges(oldVineGraph,con)
			
		if((temp[1,1] == temp[2,1])|| (temp[1,2] == temp[2,1]))
		{
			same = temp[2,1]
		}
		else
		{
			if((temp[1,1] == temp[2,2]) || (temp[1,2] == temp[2,2]))
			{
				same = temp[2,2]
			}
		}
		
		other1 = temp[1,temp[1,] != same]
		other2 = temp[2,temp[2,] != same]
	
		if(temp[1,1] == same){
			zr1 = E(oldVineGraph)[con[1]]$Copula.CondData.2
			n1 = E(oldVineGraph)[con[1]]$Copula.CondName.2
		}else{
			zr1 = E(oldVineGraph)[con[1]]$Copula.CondData.1
			n1 = E(oldVineGraph)[con[1]]$Copula.CondName.1
		}
		
		if(temp[2,1] == same){
			zr2 = E(oldVineGraph)[con[2]]$Copula.CondData.2
			n2 = E(oldVineGraph)[con[2]]$Copula.CondName.2
		}else{
			zr2 = E(oldVineGraph)[con[2]]$Copula.CondData.1
			n2 = E(oldVineGraph)[con[2]]$Copula.CondName.1
		}
		if(progress == TRUE) message(n1," + ",n2," --> ", E(mst)[i]$name)
		
		
		parameterForACopula[[i]]$zr1 = zr1
		parameterForACopula[[i]]$zr2 = zr2

		E(mst)[i]$Copula.Data.1 =  list(zr1)
		E(mst)[i]$Copula.Data.2 =  list(zr2)
		
		E(mst)[i]$Copula.CondName.2 = n1
		E(mst)[i]$Copula.CondName.1 = n2
	}

	outForACopula = lapply(X = parameterForACopula, FUN=wrapper_fit.ACopula, type,copulaSelectionBy,testForIndependence,testForIndependence.level,weights)
	
	for(i in 1:d)
	{
		E(mst)$Copula.param[[i]] = c(outForACopula[[i]]$par,outForACopula[[i]]$par2)
		E(mst)[i]$Copula.type = outForACopula[[i]]$family
		E(mst)[i]$Copula.out = list(outForACopula[[i]])
		
		E(mst)[i]$Copula.CondData.2 <- list(outForACopula[[i]]$CondOn.1)
		E(mst)[i]$Copula.CondData.1 <- list(outForACopula[[i]]$CondOn.2)
	}
	
	return(mst)
}	

buildNextGraph <- function(oldVineGraph,weights=NA)
{

	EL = get.edgelist(oldVineGraph)
	d = ecount(oldVineGraph)
	
	
	g = graph.full(d)
	V(g)$name = E(oldVineGraph)$name
	V(g)$conditionedSet = E(oldVineGraph)$conditionedSet

	if(!is.null(E(oldVineGraph)$conditioningSet)){
		V(g)$conditioningSet = E(oldVineGraph)$conditioningSet
	}
	
	for(i in 1:ecount(g)){
		
		con = get.edge(g,i)
		
		temp = get.edges(oldVineGraph,con)
		
		ok = FALSE
		
		if((temp[1,1] == temp[2,1])|| (temp[1,2] == temp[2,1])){
			ok = TRUE
			same = temp[2,1]
		}else{if((temp[1,1] == temp[2,2]) || (temp[1,2] == temp[2,2])){
				ok = TRUE
				same = temp[2,2]
			}}
		
		if(ok){
			other1 = temp[1,temp[1,] != same]
			other2 = temp[2,temp[2,] != same]
		
			if(temp[1,1] == same){
				zr1 = E(oldVineGraph)[con[1]]$Copula.CondData.2
			}else{
				zr1 = E(oldVineGraph)[con[1]]$Copula.CondData.1
			}

			if(temp[2,1] == same){
				zr2 = E(oldVineGraph)[con[2]]$Copula.CondData.2
			}else{
				zr2 = E(oldVineGraph)[con[2]]$Copula.CondData.1
			}
			
			keine_nas = !(is.na(zr1) | is.na(zr2))
			#E(g)[i]$weight = cor(x=zr1[keine_nas],y=zr2[keine_nas], method="kendall")
			E(g)[i]$weight = fasttau(zr1[keine_nas],zr2[keine_nas],weights)
			
			name.node1 = strsplit( V(g)[con[1]]$name,split=" *[,|] *")[[1]]
			name.node2 = strsplit( V(g)[con[2]]$name,split=" *[,|] *")[[1]]
			
			schnitt = c()
			
			for(j in 1:length(name.node1))
			{
				for(k in 1:length(name.node2))
				{
					if(name.node1[j] == name.node2[k])
					{
						schnitt = c(schnitt,name.node1[j])
						name.node1[j] = ""
						name.node2[k] = ""
						break
					}
				}
			}
			
			differenz = c()
			for(j in 1:length(name.node1)){
				if(name.node1[j] != ""){
					differenz = c(differenz, name.node1[j])
				}
			}
			for(j in 1:length(name.node2)){
				if(name.node2[j] != ""){
					differenz = c(differenz, name.node2[j])
				}
			}
			
			E(g)[i]$name = paste(
					paste(differenz, collapse= ","),
					paste(schnitt, collapse= ","),
					sep= " | ")
			
			l1 = c(V(g)[con[1]]$conditionedSet,V(g)[con[1]]$conditioningSet)
			l2 = c(V(g)[con[2]]$conditionedSet,V(g)[con[2]]$conditioningSet)
			
			out = intern_SchnittDifferenz(l1,l2)
			
			suppressWarnings({E(g)$conditionedSet[i] = list(out$differenz)})
			suppressWarnings({E(g)$conditioningSet[i]  = list(out$schnitt)})
		}
		
		E(g)[i]$todel = !ok
	}
	
	E(g)$tau = E(g)$weight
	
	g = delete.edges(g, E(g)[E(g)$todel])
	
	return(g)
}

wrapper_fit.ACopula <- function(parameterForACopula,type,...)
{
	return(fit.ACopula(parameterForACopula$zr1,parameterForACopula$zr2,type,...))
}

intern_SchnittDifferenz = function(liste1,liste2){
	out = list()	
	out$schnitt = c()
	out$differenz = c()
	
	for(j in 1:length(liste1)){
		for(k in 1:length(liste2)){
			if(!is.na(liste2[k]) && liste1[j] == liste2[k]){
				out$schnitt = c(out$schnitt, liste1[j])
				liste1[j] = NA
				liste2[k] = NA
				break
			}
		}
	}
	
	for(j in 1:length(liste1)){
		if(!is.na(liste1[j])){
			out$differenz = c(out$differenz, liste1[j])
		}
	}
	for(j in 1:length(liste2)){
		if(!is.na(liste2[j])){
			out$differenz = c(out$differenz, liste2[j])
		}
	}
	
	return(out)
}

fit.ACopula <- function(u1,u2,familyset=NA,selectioncrit="AIC",indeptest=FALSE,level=0.05,weights=NA) 
{

	out=BiCopSelect(u1,u2,familyset,selectioncrit,indeptest,level,weights=weights)
	if(out$family%in%c(23,24,26:30,124,224))
	{
		out$family=out$family+10
	}
	else if(out$family%in%c(33,34,36:40,134,234))
	{
		out$family=out$family-10
	}
      out$CondOn.1 = .C("Hfunc1",as.integer(out$family),as.integer(length(u1)),as.double(u1),as.double(u2),as.double(out$par),as.double(out$par2),as.double(rep(0,length(u1))),PACKAGE='VineCopula')[[7]]
      out$CondOn.2 = .C("Hfunc2",as.integer(out$family),as.integer(length(u1)),as.double(u2),as.double(u1),as.double(out$par),as.double(out$par2),as.double(rep(0,length(u1))),PACKAGE='VineCopula')[[7]]



  return(out)

}

as.RVM = function(RVine){

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
					cty = corresppondingTypes[[n-i+1]][[j]]
					if(cs[1] == w){
						M[i,k] = cs[2]
					  Type[i,k] = cty #Mathias 21.3.
						break
					} else if(cs[2] == w){
						M[i,k] = cs[1]
						if(any(cty == c(23,24,26))){ Type[i,k] = cty+10} #Mathias 21.3.
						if(any(cty == c(33,34,36))){ Type[i,k] = cty-10} #Mathias 21.3.
						if(!any(cty == c(23,24,26,33,34,36))){ Type[i,k] = cty} #Mathias 21.3.
						break
					}
				}
				Param[i,k] = corresppondingParams[[n-i+1]][[j]][1]
				Params2[i,k] = corresppondingParams[[n-i+1]][[j]][2]
# changed Mathias 21.3.				Type[i,k] = corresppondingTypes[[n-i+1]][[j]]

				conditionedSets[[n-i+1]][[j]] = NULL
				corresppondingParams[[n-i+1]][[j]] = NULL
				corresppondingTypes[[n-i+1]][[j]] = NULL
			}
		}

	}

	M = M#+1
	M[is.na(M)]=0
	Type[is.na(Type)]=0

	return(RVineMatrix(M, family = Type, par = Param, par2 = Params2, names = names))

}
