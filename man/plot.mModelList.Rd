\name{plot.mModelList}
\alias{plot.mModelList}
\title{
Plots a graphical visualization for every model from list of models.
}
\description{
The function \code{plot.mModelList()} creates a grid of panels and then plots fitted models in consecutive panels. The \code{plot.mModel()} function is used to plot single model.
}
\usage{
plot.mModelList(x, ...)
}

\arguments{
  \item{x}{An object of the class \code{mModelList},
  }
  \item{\dots}{Graphical arguments that will be passed to underlying \code{plot()} function.
  }
}
\details{
The argument \code{x} specifies a list of models. If these models differ both in number of components and in model structures, 
then in the resulting grid of panels columns will correspond to different model structures while rows will correspond to different number of components. 

If considered models differ only in number of components or in model structures then the grid of panels will be as close to square as possible and consecutive panels will present consecutive models from model list.
}
\references{
  http://bgmm.molgen.mpg.de
}
\seealso{\code{\link{plot.mModel}}}
\author{
Przemyslaw Biecek
}
\examples{
 simulated = simulateData(d=2, k=3, n=100, m=60, cov="0", within="E", n.labels=2)
 models1=mModelList(X=simulated$X, knowns=simulated$knowns, B=simulated$B, 
             kList=3:4, mean=c("D","E"), between="D", within="D", 
             cov="0", funct=belief)
 plot(models1)

## Do not run
## could take more than one minute
# simulated = simulateData(d=2, k=3, n=300, m=60, cov="0", within="E", n.labels=2)
# 
# models1=mModelList(X=simulated$X, knowns=simulated$knowns, B=simulated$B, 
#             kList=3, mean=c("D","E"), between=c("D","E"), within=c("D","E"), 
#             cov=c("D","0"), funct=belief)
# plot(models1)
#
# models2 = beliefList(X=simulated$X, knowns=simulated$knowns, B=simulated$B,
#                kList=2:7, mean="D", between="D", within="E", cov="0")
# plot(models2)
# 
# models3 = beliefList(X=simulated$X, knowns=simulated$knowns, B=simulated$B,
#               kList=2:7, mean="D")
# plot(models3)
}


