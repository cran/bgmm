\name{plot.mModel}
\alias{plot.mModel}
\title{
Plotting a graphical visualization of a Gaussian model or a list of models for one-dimensional or two-dimensional data
}
\description{
The generic function \code{plot} is used to visualize the data set and Gaussian model components fitted to this data. On the resulting plot the observations without labels are presented with  black points, whereas the labeled observations are marked by different colors and different symbols. This function works only for one or two dimensional data. The fitted Gaussian components are represented by ellipses int the two-dimensional case and by densities in the one dimensional case.
}
\usage{
plot.mModel(x, ...)
}

\arguments{
  \item{x}{an object of the class \code{mModel}.
  }
  \item{\dots}{graphical arguments that are passed to the underlying \code{plot()} function.
  }
}
\details{
  For one dimensional data the width of the density corresponds to standard deviation of the fitted Gaussian component. Fitted means are marked by veritical dashed lines.
  
  For two dimensional data ellipses represents covariances for the corresponding model components.   
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
