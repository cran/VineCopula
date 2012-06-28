###########################################
# TauMatrix
#
# Input:
# data		data matrix
#
# Output:
# ktauMarix	Kendall's tau matrix
############################################

TauMatrix <- function(data)
{
	data=as.matrix(data)
	if(any(data>1) || any(data<0)) stop("Data has be in the interval [0,1].")
	d=dim(data)[2]
	N=dim(data)[1]

	ktau=rep(0,d*(d-1)/2)

	out <- .C("ktau_matrix",
		as.double(data),
		as.integer(d),
		as.integer(N),
		as.double(ktau),
		PACKAGE = 'VineCopula')
	
	ktau=out[[4]]

	ktauMatrix=matrix(1,d,d)
	k=1
	for(i in 1:(d-1))
	{
		for(j in (i+1):d)
		{
			ktauMatrix[i,j]=ktau[k]
			ktauMatrix[j,i]=ktau[k]
			k=k+1
		}
	}
	if(!is.null(colnames(data)))
	{
		rownames(ktauMatrix)=colnames(ktauMatrix)=colnames(data)
	}

return(ktauMatrix)
}
