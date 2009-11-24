\name{init.model.params}
\alias{init.model.params}
\alias{init.model.params.knowns}
\title{
Initiation for model parameters 
}
\description{
Methods for initiation of the EM algorithm. Two initiation algorithms are implemented. One which take into account only labeled cases (siutable for datasets with high percentage of labeled cases) and one which does not take any information information about observations labeling.
}
\usage{
init.model.params(X = NULL, knowns = NULL, class = NULL, 
    k = length(unique(class)), method = "all", B = P, P = NULL, 
    clusterAssigment = TRUE)
}
\arguments{
  \item{X}{A data.frame with the unlabeled observations, the rows correspond to the observations and the columns to the dimensions of the data}
  \item{knowns}{A data.frame with the labeled observations, the rows correspond to the observations and the columns to the dimensions of the data}
  \item{class, B, P}{Class is a vector of labels for known observations. If not specified then it will be derived from B or P matrices by MAP rule,  }
  \item{k}{Desired number of components,  }
  \item{method}{Method for parameter initialization, one of following c("knowns","all"), see details.  }
  \item{clusterAssigment}{This argument specify if labels from initial clustering should be agreed with initial labels of labeled cases. Should be set to FALSE for unsupervised learning while there is no labeled cases and FALSE for method="knowns" while then it is no need for cluster agreeing.
  }
}
\details{
If method='knowns' then the initialization is based only on the labeled cases. The initial model parameters for each component are estimated in one step from the observations assigned to this component, 
      as in the fully supervised learning. 

If method='all' (default) then the initialization is based on all observations. In this case, to obtain the initial set of model components, we start by clustering the data using k-means algorithm (repeated 10 times to get stable results). 
      For one dimensional data the clusters are identified by dividing the data into intervals of equal sizes. 
      Each resulting cluster corresponds to one initial model component. 
      
For the partially and semi-supervised methods, correspondence of labels from clustering algorithm and labels for known cases rises a technical problem: 
      the cluster corresponding to component y should be as close as possible to the set of labeled observations with label y. 
      
      Note that for the unsupervised modeling this problem is irrelevant and any cluster may be used to initialize any component.
      
      To label the clusters with model components, we use a greedy heuristic. The heuristic calculates weighted distances between all possible pairs of cluster centers and sets of observations grouped by their labels. 
      In each step, the pair with a minimal distance is chosen. For the chosen pair, the cluster is first labeled the same as the set of observations. 
      Next, this pair is removed and the heuristic choice repeats. 
      
      If  method='all' and data is one dimensional then the initial clustering is calculated by dividing data into k equal subsets of observations, subset are separated by empirical quantiles c(1/2k, 3/2k, 5/2k, ..., (2k-1)/2k).}
\value{
A list with the following elements:
\item{pi}{ A vector of length k with the initial values for mixing proportions.}
\item{mu}{ A matrix with the means' vectors with the initial values for k components.}
\item{cvar}{ A three-dimensional matrix with the covariance matrices with the initial values for k components.} 
}
\references{
  http://bgmm.molgen.mpg.de
}
\author{
Przemyslaw Biecek
}
\examples{
 data(genotypes)
 initial.params = init.model.params(X=genotypes$X, knowns=genotypes$knowns,
									 class = genotypes$labels)
 str(initial.params)
}


