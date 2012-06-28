BiCopSim = function(N,family,par,par2=0)
{

	if(!any(family %in% c(2,7:10,17:20,27:30,37:40)))
        tmp = .C("pcc",
        as.integer(N),
        as.integer(2),
        as.integer(family),
        as.integer(1),
        as.double(par),
        as.double(0),
        as.double(rep(0,N*2)),
        PACKAGE='VineCopula')[[7]]
    else tmp = .C("pcc",
         as.integer(N),
         as.integer(2),
         as.integer(family),
         as.integer(1),
         as.double(par),
         as.double(par2),
         as.double(rep(0,N*2)),
         PACKAGE='VineCopula')[[7]]
    U <- matrix(tmp,ncol=2)
    
  return(U)

}
