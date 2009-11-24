\name{chooseModels}
\alias{chooseModels}
\alias{chooseOptimal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subset of models.
}
\description{
The function chooseModels extracts sublist of models which match requirements for number of components and for models structure. The function chooseOptimal returns the best model with respect to given model selection criteria. 
}
\usage{
chooseModels(models, kList = NULL, struct = NULL)

chooseOptimal(models, penalty=2)
}

\arguments{
  \item{models}{Object of the class \code{\link{mModelList}}
  }
  \item{kList}{The vector which specifies requested numbers of Gaussian components. }
  \item{struct}{The vector which specifies names of requested structures. }
  \item{penalty}{The penalty for GIC criteria. This parameter can be a single number or one of the strings c("BIC", "AIC").   }
}
\details{
  The function chooseModels() extracts a sublist of models from the \code{models} argument.
  The returned sublist is composed of models that simultaneously satisfy restrictions for number of Gaussian components defined by \code{kList} and restrictions for model structure defined by \code{struct}.
  If argument \code{kList} is set to NULL then no restrictions for number of components is applied, same with \code{struct} argument.

  The function chooseOptimal() returns the single model of the class \code{\link{mModel}} which has best (smallest) GIC score.
}
\value{
The object of the class \code{\link{mModelList}} or \code{\link{mModel}}.
}
\references{
  http://bgmm.molgen.mpg.de
}
\author{
Przemyslaw Biecek
}
\examples{
 simulated = simulateData(d=2, k=3, n=100, m=50, cov="0", within="E", n.labels=2)
 
 models3 = beliefList(X=simulated$X, knowns=simulated$knowns, B=simulated$B,
               kList=2:4, mean="D", within="D")
 plotGIC(models3, penalty="BIC")

## Do not run
## could take more than one minute
# simulated = simulateData(d=2, k=3, n=300, m=60, cov="0", within="E", n.labels=2)
# 
# models3 = beliefList(X=simulated$X, knowns=simulated$knowns, B=simulated$B,
#               kList=2:7, mean="D")
# plotGIC(models3, penalty="BIC")
#
# models4 = chooseModels(models3, kList=2:5, struct=c("DDDD","DDED","DDE0"))
# plot(models4)
# plotGIC(models4, penalty="BIC")
#
# model4 = chooseOptimal(models3, "BIC")
# plot(model4)
}


