\name{plotGIC}
\alias{plotGIC}
\title{
GIC scores plot.
}
\description{
The function \code{plotGIC()} plots GIC scores for collection of models.
}
\usage{
plotGIC(models, penalty = 2, plot.it = TRUE)
}

\arguments{
  \item{models}{an object of the class \code{mModelList} or matrix of GIC scores.
  }
  \item{penalty}{penalty for GIC criteria. This parameter can be a single number or one of the strings c("BIC", "AIC").     
  }
  \item{plot.it}{logical value which specify does the chart with GIC scores should be plotted.
  }
}
\details{
The function \code{plotGIC()} calculates GIC scores for each model from the \code{models} list and then 
plots dotchart with GIC scores for models specified by the \code{models} argument. 
The function \code{plotGIC()} returns a matrix with calculated GIC scores. This matrix or its submatrix can be used in next call the \code{plotGIC()} function as \code{models} argument.

The best model, i.e. model with smallest GIC scores is marked by star. 
}
\value{
The matrix with GIC scores calculated for list of models specified in \code{models} argument. 
Columns in resulting matrix corresponds to different number of components while rows corresponds to different structures.
}
\references{
  http://bgmm.molgen.mpg.de
}
\author{
Przemyslaw Biecek
}
\examples{
 simulated = simulateData(d=2, k=3, n=100, m=60, cov="0", within="E", n.labels=2)
 models1=mModelList(X=simulated$X, knowns=simulated$knowns, B=simulated$B, 
             kList=3:4, mean=c("D","E"), between="D", within="D", 
             cov="0", funct=belief)
 plotGIC(models1, penalty="BIC")

## Do not run
## could take more than one minute
# simulated = simulateData(d=2, k=3, n=300, m=60, cov="0", within="E", n.labels=2)
# 
# models1=mModelList(X=simulated$X, knowns=simulated$knowns, B=simulated$B, 
#             kList=3, mean=c("D","E"), between=c("D","E"), within=c("D","E"), 
#             cov=c("D","0"), funct=belief)
# plotGIC(models1, penalty="BIC")
# 
# models2 = beliefList(X=simulated$X, knowns=simulated$knowns, B=simulated$B,
#                kList=2:7, mean="D", between="D", within="E", cov="0")
# plotGIC(models2, penalty="BIC")
# 
# models3 = beliefList(X=simulated$X, knowns=simulated$knowns, B=simulated$B,
#               kList=2:7, mean="D")
# plotGIC(models3, penalty="BIC")
}


