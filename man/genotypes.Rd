\name{genotypes}
\alias{genotypes}
\docType{data}
\title{
Fluorescence signals corresponding to a given allele for 333 SNPs.
}
\description{
This dataset describes 333 SNPs. Rows in the slots \code{X} and \code{knowns} 
correspond to different SNPs, while values in columns represents an intensity of 
the fluorescence signal corresponding to the selected alleles of a given SNP. 

To prepare the input for the partially and semi-supervised modeling methods 15 SNPs was 
separated by choosing at random five SNPs per each cluster. For the partially supervised 
methods, the belief/plausibility values of the most probable type reflect the high 
certainty of the labels and are set to 0.95, and of the remaining types are equal 0.025. 
The remaining 318 SNPs are kept unlabeled.
}
\usage{data(genotypes)}
\format{
  X       : matrix 318x2
  knowns  : matrix 15x2
  B       : matrix 15x3
  labels  : vector 15
}
\details{
  \code{genotypes} A list, where slot X corresponds to fluorescence signals of unlabeled SNPs, slot knowns corresponds to fluorescence signal of lablled SNPs, slot B corresponds to belief matrix while labels contains true labels for labeled SNPs.
}

\references{
   Takitoh, S.  Fujii, S.  Mase, Y.  Takasaki, J.  Yamazaki, T.  Ohnishi, Y.  Yanagisawa, M.  Nakamura, Y.  Kamatani, N., Accurate automated clustering of two-dimensional data for single-nucleotide polymorphism genotyping by a combination of clustering methods: evaluation by large-scale real data, Bioinformatics (2007) Vol. 23, 408--413.
}
\examples{
library(bgmm)
data(gnotypes)
}
\keyword{datasets}
