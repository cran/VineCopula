\name{copulaFromFamilyIndex}
\alias{copulaFromFamilyIndex}

\title{
Construction of a copula object from a VineCopula family index
}
\description{
A VineCopula family index along with its parameters is used to construct a corresponding \code{\linkS4class{copula}} object. 
}
\usage{
copulaFromFamilyIndex(family, par, par2 = 0)
}

\arguments{
  \item{family}{a family inde aS DEFINED IN \code{\link{VineCopula-package}}}
  \item{par}{its first parameter value}
  \item{par2}{if present, its second parameter}
}
\value{An object inherting \code{\linkS4class{copula}} corresponding to the specific family.}

\author{Benedikt Graeler}

\examples{
# normalCopula with parameter 0.5
copulaFromFamilyIndex(1,0.5)

# rotated Tawn T2 copula with parameters 
copulaFromFamilyIndex(224, -2, 0.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
