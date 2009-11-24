\name{mModelList}
\alias{mModelList}
\alias{beliefList}
\alias{softList}
\alias{semisupervisedList}
\alias{supervisedList}
\alias{unsupervisedList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gaussian mixture model fitting for collection of models.
}
\description{
These functions fit collection of models of one particular variant. 
Models may differ in number of Gaussian components or in model structures.
}
\usage{
mModelList(X, knowns, B = NULL, P = NULL, class = NULL, kList = ncol(B), 
    init.params = NULL, stop.likelihood.change = 10^-5, stop.max.nsteps = 100, 
    trace = FALSE, mean = c("D", "E"), between = c("D", "E"), within = c("D", 
        "E"), cov = c("D", "0"), funct = belief, clusterAssigment = TRUE,
		all.possible.permutations = FALSE)
        
beliefList(..., funct=belief)

softList(..., funct=soft)

semisupervisedList(..., funct=semisupervised)

unsupervisedList(X, kList = 2, ...)
}

\arguments{
  \item{X}{a data.frame with the unlabeled observations, the rows correspond to the observations and the columns to the dimensions of the data }
  \item{knowns}{a data.frame with the labeled observations, the rows correspond to the observations and the columns to the dimensions of the data  }
  \item{B}{a beliefs matrix with the distribution of beliefs for the
  labeled observations. The rows of B correspond to the labeled
  observations and are assumed to be given in the same order as in the
  argument knowns. The columns correspond to the model components. If matrix B is provided then the number of columns has to be not greater than min(kList)   }
  \item{P}{a  matrix  with plausibilities, i.e., weights for the prior
  probabilities for the labeled observations. If matrix P is provided then the number of columns has to be not greater than min(kList)  }
  \item{class}{a vector of classes for the labeled observations. If vector class is provided then the number of unique values has to be not greater than min(kList)  }
  \item{kList}{vector or list with numbers of components to fit, by default equal to the number of columns of B  }
  \item{init.params}{initial values for model parameter that will be passed to funct function  }
  \item{stop.likelihood.change, stop.max.nsteps, trace}{parameters for the EM  algorithm that will be passed to funct function    }
  \item{mean, between, within, cov}{four vectors which define model structures to fit. If mean="E" then only models with constraint over mean will be considered but if mean=c("E", "D") than both model with constraints over mu parameter and model without constraints will be fitted }
  \item{funct}{function which fit variant of Gaussian mixture models, one of the: belief, soft, semisupervised or unsupervised functions
  }
  \item{clusterAssigment}{parameter for the function init.model.params() which specifies if cluster labels should be agreed with knowns labels  }
 \item{\dots}{arguments that will be passed to function \code{funct}.
  }
  \item{all.possible.permutations}{If TRUE then all possible permutations of initial parameters are considered, the resulting model is the model with highest likelihood.}
}
\details{
  Arguments \code{kList} and \code{mean}, \code{between}, \code{within}, \code{cov} define list of models to fit. All combinations of specified model sizes and model structures will be considered. List of fitted models will be returned as result. 
  
  The argument  \code{funct} defines which variant of Gaussian mixture models should be used for fitting. One can use 
  wrappers \code{beliefList()}, \code{softList()}, \code{semisupervisedList()}, \code{unsupervisedList()} with prespecified argument \code{funct}.
}
\value{
An object of the class mModelList, with the following slots:
\item{models}{ a list of models, each of the class \code{\link{mModel}}}
\item{loglikelihoods}{ a vector with log likelihoods for models from list \code{models}}
\item{names}{ a vector with names for models from list \code{models}} 
\item{params}{ a vector with number of parameters of models from list \code{models}}
\item{kList}{as kList input argument}
}
\references{
  http://bgmm.molgen.mpg.de
}
\author{
Przemyslaw Biecek
}
\seealso{\code{\link{mModel}}, \code{\link{getModelStructure}}}
\examples{
 simulated = simulateData(d=2, k=3, n=100, m=60, cov="0", within="E", n.labels=2)
 
 models1=mModelList(X=simulated$X, knowns=simulated$knowns, B=simulated$B, 
             kList=3:4, mean=c("D","E"), between="D", within="D", 
             cov="0", funct=belief)

 plot(models1)
 plotGIC(models1, penalty="BIC")

## Do not run
## could take more than one minute
# simulated = simulateData(d=2, k=3, n=300, m=60, cov="0", within="E", n.labels=2)
# 
# models1=mModelList(X=simulated$X, knowns=simulated$knowns, B=simulated$B, 
#             kList=3, mean=c("D","E"), between=c("D","E"), within=c("D","E"), 
#             cov=c("D","0"), funct=belief)
# plot(models1)
# plotGIC(models1, penalty="BIC")
# 
# models2 = beliefList(X=simulated$X, knowns=simulated$knowns, B=simulated$B,
#                kList=2:7, mean="D", between="D", within="E", cov="0")
# plot(models2)
# plotGIC(models2, penalty="BIC")
# 
# models3 = beliefList(X=simulated$X, knowns=simulated$knowns, B=simulated$B,
#               kList=2:7, mean="D")
# plotGIC(models3, penalty="BIC")
}


