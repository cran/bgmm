\name{getModelStructure}
\alias{getModelStructure}
\title{
Model structure description.
}
\description{
This function creates object which describes constraints on model structure.
}
\usage{
getModelStructure(mean = "D", between = "D", within = "D", cov = "D")
}
\arguments{
  Each argument is single character (E, equal) or (D, different). If set to E then given property should be constrained.  By default all are set to D.

  \item{mean}{forces equality of the means among the components,}
  \item{between}{forces equality of the covariance matrices among the components, }
  \item{within}{forces equality of variances within each covariance matrix i to some constant vi and equality of covariances to some constant wi}
  \item{cov}{forces equality of covariances within each covariance matrix i to 0,}
}
\value{
List of four elements.
}
\references{
  http://bgmm.molgen.mpg.de
}
\author{
Ewa Szczurek
}
\examples{
 getModelStructure()
 getModelStructure(mean="E")
}


