BiCopPar2Beta <- function(family,par,par2=0)
{
	blomBeta=4*BiCopCDF(0.5,0.5,family,par,par2)-1
	
	return(blomBeta)
}