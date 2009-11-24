\name{Supplementary functions}
\alias{determinant.numeric}
\alias{map}
\alias{loglikelihood.mModel}
\alias{DEprobs}
\title{
Set of supplementary functions for bgmm package.
}
\description{
Set of supplementary functions for bgmm package.
}
\usage{
 determinant.numeric(x, logarithm = TRUE, ...)
 
 map(B) 
 
 loglikelihood.mModel(model, X) 
 
 DEprobs(model, verbose=FALSE)
}

\arguments{
  \item{x}{single number}
  \item{X}{a data.frame with the unlabeled observations, the rows correspond to the observations and the columns to the dimensions of the data}
  \item{B}{a beliefs matrix with the distribution of beliefs for the
  labeled observations. }
  \item{model}{an object of the class mModel}
  \item{logarithm, verbose, \dots}{these arguments will be ignored  }
}
\references{
  http://bgmm.molgen.mpg.de
}
\author{
Przemyslaw Biecek
}
\examples{
data(genotypes)

map(genotypes$B)

}
