\name{mModel}
\alias{mModel}
\alias{belief}
\alias{soft}
\alias{semisupervised}
\alias{supervised}
\alias{unsupervised}
\title{
Gaussian mixture model fitting.
}
\description{
These functions fit different variants of Gaussian mixture models. These variants 
differ in the amount of knowledge incorporated into the fitting procedure.
}
\usage{
belief(X, knowns, B = NULL, k = ifelse(!is.null(B), ncol(B), 
    ifelse(!is.null(P), ncol(P), length(unique(class)))), P = NULL, 
    class = map(B), init.params = init.model.params(X, knowns, 
        B = B, P = P, class = class, k = k), model.structure = getModelStructure(), 
    stop.likelihood.change = 10^-5, stop.max.nsteps = 100, trace = FALSE, 
    b.min = 0.025, clusterAssigment = TRUE, all.possible.permutations=FALSE)
    
soft(X, knowns, P = NULL, k = ifelse(!is.null(P), ncol(P), 
    ifelse(!is.null(B), ncol(B), length(unique(class)))), B = NULL, 
    class = NULL, init.params = init.model.params(X, knowns, 
        class = class, B = P, k = k, clusterAssigment = clusterAssigment), 
    model.structure = getModelStructure(), stop.likelihood.change = 10^-5, 
    stop.max.nsteps = 100, trace = FALSE, clusterAssigment = TRUE, b.min = 0.025, 
	all.possible.permutations=FALSE)    
    
semisupervised(X, knowns, class = NULL, k = ifelse(!is.null(class), 
    length(unique(class)), ifelse(!is.null(B), ncol(B), ncol(P))), 
    B = NULL, P = NULL, ..., clusterAssigment = TRUE, 
	all.possible.permutations=FALSE)    
    
supervised(knowns, class = NULL, k = length(unique(class)), B = NULL, P = NULL,
     clusterAssigment = TRUE, model.structure = getModelStructure())

unsupervised(X, k, init.params=init.model.params(X, knowns=NULL, k=k, clusterAssigment=FALSE), 
      model.structure=getModelStructure(), stop.likelihood.change=10^-5, stop.max.nsteps=100, 
      trace=FALSE, clusterAssigment=FALSE) 
}

\arguments{
  \item{X}{A data.frame with unlabeled observations. Rows correspond to  observations and  columns to  dimensions of the data}
  \item{knowns}{A data.frame with  labeled observations. Rows correspond to  observations and columns to  dimensions of the data}
  \item{B}{A beliefs matrix with the distribution of beliefs for 
  labeled observations. Rows of B correspond to  labeled
  observations and are assumed to be given in the same order as in the
  argument knowns. Columns correspond to model components.  }
  \item{P}{A  matrix  with plausibilities, i.e., weights for prior
  probabilities for labeled observations. }
  \item{class}{A vector of classes for labeled observations}
  \item{k}{Number of components, by default equal to the number of columns of B}
  \item{init.params}{Initial values for model parameter estimates
  (means, variances and mixing proportions), bu default derived with the
  use of the \code{\link{init.model.params}} function}
  \item{stop.likelihood.change, stop.max.nsteps}{Parameters for the EM
  algorithms defining the stop criteria, i.e., the minimum required improvement in loglikelihood and the maximum number of steps}
  \item{trace}{If trace=TRUE then loglikelihoods for every step of EM algorithm will be presented
  }
  \item{model.structure}{An object returned by the \code{\link{getModelStructure}} function, which defines the constraints for the parameters of the model to be fitted.}
  \item{b.min}{This argument will be passed to function \code{\link{init.model.params}}  }
  \item{\dots}{These arguments will be passed to function \code{\link{init.model.params}}    }
  \item{clusterAssigment}{Parameter for the function init.model.params() which specifies if cluster labels should be agreed with knowns labels  }
  \item{all.possible.permutations}{If TRUE then all possible permutations of initial parameters are considered, the resulting model is the model with highest likelihood.}
}
\details{
  In the belief() function, if the argument B is not provided, it is
  first set to the argument P. If the argument P is not
  provided, B is derived from the class argument by function get.simple.beliefs() 
  which assigns 1-(k-1)*b.min to the component corresponding to class argument and 
  b.min to other components.

  In the soft() function, if the argument P is not provided, it is
  first set to the argument B. If the argument B is not
  provided, P is derived from the class argument as in the belief()
  function.

  In the supervised() function,  if the argument class is not provided,
  it is derived from argument B or P, taking the class of each
  observation as its most believed or plausible component.
  
  The number of columns in the beliefs matrix B or in the matrix of
  plausibilities P may be lower than the number of model components
  defined by the argument k.  Such situation corresponds to the scenario
  when the user does not know any examples of observations from some
  component. In other words, this component is not used as a label for
  any observation, and thus can be omitted from the beliefs matrix. An
  equivalent would be to include a column for this component and fill it
  with beliefs/plausibilities equal to 0.

  Slots in the returned object are listed in the next section.
  The value returned by the functions differs accordingly: the belief() function returns an object with a slot B. 
  Instead of B, soft() returns an object with a slot P, and supervised() and semisupervised() with a slot class. 
  The object returned by the supervised() function does not have a slot X.
}
\value{
An object of the class mModel, with the following slots:
\item{pi}{ a vector with the fitted mixing proportions}
\item{mu}{ a matrix with the means' vectors, fitted for all components}
\item{cvar}{ a three-dimensional matrix with the covariance matrices, fitted for all components} 
\item{X}{ the unlabeled observations}
\item{knowns}{the labeled observations}
\item{B}{the beliefs matrix}
\item{n}{ the number of all observations}
\item{m}{ the number of the unlabeled observations}
\item{k}{the number of fitted model components}
\item{d}{the data dimension}
\item{likelihood}{log-likelihood of the fitted model}
\item{n.steps}{the number of steps performed by the EM algorithm}
\item{model.structure}{the set of constraints kept during the fitting process}
}
\references{
  http://bgmm.molgen.mpg.de
}
\author{
Przemyslaw Biecek
}
\examples{
data(genotypes)

modelSupervised = supervised(knowns=genotypes$knowns, 
            class=genotypes$labels)
plot(modelSupervised)

modelSemiSupervised = semisupervised(X=genotypes$X, 
            knowns=genotypes$knowns, class = genotypes$labels)
plot(modelSemiSupervised)

modelBelief = belief(X=genotypes$X, 
            knowns=genotypes$knowns, B=genotypes$B)
plot(modelBelief)

modelSoft = soft(X=genotypes$X, 
            knowns=genotypes$knowns, P=genotypes$B)
plot(modelSoft)

modelUnSupervised = unsupervised(X=genotypes$X, k=3)
plot(modelUnSupervised)

}
