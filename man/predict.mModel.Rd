\name{predict.mModel}
\alias{predict.mModel}
\title{
Predictions for fitted Gaussian component model
}
\description{
For every row in the X matrix the posterior probability of belonging to class i is calculated.
}
\usage{
predict.mModel(object, X, knowns = NULL, B = NULL, P = NULL, ...)
}
\arguments{
  \item{object}{an object of the class \code{mModel},   }
  \item{X}{matrix or data.frame with number of columns equal to object$d.   }
  \item{knowns}{a data.frame or matrix with the labeled observations, is specified then B or P need to be specified}
  \item{P}{a  matrix  with plausibilities for object 'knowns' }
  \item{B}{a  matrix  with beliefs for object 'knowns' }
  \item{\dots}{all remaining arguments will be neglected.   }
}
\details{
  The matrix tij of posterior probabilities is calculated as normalized products of priors pi's and density of model components in values specified be rows of matrix \code{X}.

  If arguments \code{knowns} and \code{B} are specified, then the priors's for objects in \code{knowns} are replaced by belief matrix \code{B}.
  If arguments \code{knowns} and \code{P} are specified, then the priors's for objects in \code{knowns} are multiplied by plausibility matrix \code{P}.
}
\value{
An list with the following elements:
\item{tij}{A matrix with number of rows equal number of rows in matrix \code{X} and number of columns equal to number of components in model 
    defined by argument object. Values in this matrix are the posterior probabilities that observation i belongs to component j.}
\item{class}{Classification obtained by MAP rule derived from posterior distribution.}
}
\references{
  http://bgmm.molgen.mpg.de
}
\author{
Przemyslaw Biecek
}
\seealso{\code{\link{belief}}}
\examples{
 data(genotypes)

 modelSoft = soft(X=genotypes$X, knowns=genotypes$knowns, P=genotypes$B)

 preds = predict(modelSoft, X = genotypes$X)
 str(preds)
}


