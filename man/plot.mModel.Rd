\name{plot.mModel}
\alias{plot.mModel}
\title{
Plots a graphical visualization of Gaussian component models for single model for one-dimensional or two-dimensional data.
}
\description{
This function will work only for one od two dimensinal data.

The generic function \code{plot} is used to visualize data set and fitted Gaussian components. 
The observations without labels are denoted by black points, the labeled observations are marked by different colors and different symbols. The fitted Gaussian components are represented by ellipses for two-dimensional case and by densities in one dimensional case.
}
\usage{
plot.mModel(x, ...)
}

\arguments{
  \item{x}{An object of the class \code{mModel}.
  }
  \item{\dots}{Graphical arguments that will be passed to underlying \code{plot()} function.
  }
}
\details{
  For one dimensional data the width of density corresponds to sd for fitted Gaussian component. Fitted means are marked by veritical dashed lines.
  For two dimensional data ellipses represents covariances for corresponding model components.   
}
\references{
  http://bgmm.molgen.mpg.de
}
\author{
Przemyslaw Biecek
}
\examples{
 data(genotypes)
 modelSupervised = supervised(knowns=genotypes$knowns, class=genotypes$labels)
 plot(modelSupervised)

 # semi-supervised modeling
 modelSemiSupervised = semisupervised(X=genotypes$X, knowns=genotypes$knowns,
                           class = genotypes$labels)
 plot(modelSemiSupervised)

 # belief-based modeling
 modelBelief = belief(X=genotypes$X, knowns=genotypes$knowns, B=genotypes$B)
 plot(modelBelief)

 # soft-label modeling
 modelSoft = soft(X=genotypes$X, knowns=genotypes$knowns, P=genotypes$B)
 plot(modelSoft)

 # unsupervised modeling
 modelUnSupervised = unsupervised(X=genotypes$X, k=3)
 plot(modelUnSupervised)
}


