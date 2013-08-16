# wrapping C functions to be used in spcopula

RHfunc1 <- function(fam, n, u, param) {
  .C("Hfunc1", as.integer(fam), as.integer(n), as.double(u[,2]), as.double(u[,1]), 
     as.double(param[1]), as.double(param[2]), as.double(rep(0, n)), 
     PACKAGE = "VineCopula")
}

RHfunc2 <- function(fam, n, u, param) {
  .C("Hfunc2", as.integer(fam), as.integer(n), as.double(u[,1]), as.double(u[,2]), 
     as.double(param[1]), as.double(param[2]), as.double(rep(0, n)), 
     PACKAGE = "VineCopula")
}
    
RLL_mod_separate <- function(fam, n, u, param) {
  .C("LL_mod_seperate", as.integer(fam), as.integer(n), as.double(u[,1]), 
     as.double(u[,2]), as.double(param[1]), as.double(param[2]), 
     as.double(rep(0, n)), PACKAGE = "VineCopula")
}

RarchCDF <- function(fam, n, u, param) {
  .C("archCDF", as.double(u[,1]), as.double(u[,2]), as.integer(n), as.double(param),
     as.integer(fam), as.double(rep(0, n)), PACKAGE = "VineCopula")
}

Rpcc <- function(fam, n, param) {
  .C("pcc", as.integer(n), as.integer(2), as.integer(fam), as.integer(1), 
     as.double(param[1]), as.double(param[2]), as.double(rep(0, n * 2)), 
     PACKAGE = "VineCopula")
}