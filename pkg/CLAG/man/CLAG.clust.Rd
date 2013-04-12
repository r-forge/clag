\name{CLAG.clust}
\alias{CLAG.clust}
\title{Run cluster analysis with CLAG}
\description{
  This function computes clusters with the CLAG algorithm.
  CLAG calculates a global distribution for all the values in the
  matrix, which implies variables are comparable magnitudes.
  To understand what the parameters are, please refer to the reference
  article (see below).
}
\usage{
CLAG.clust <- function(DATA, delta=0.05, threshold=0, analysisType=1, verbose=FALSE, keepTempFiles=FALSE)
}
\arguments{
\item{DATA}{DATA is either a matrix (or data frame) or a list composed of
 - a matrix M
 - a vector colIds with column ids (same length as ncol(M))
 - a vector rowIds with row ids (same length as nrow(M)).
}
\item{delta}{Delta is a parameter which influences how CLAG decides whether two values are near or not. Examples: 0.05, 0.1, 0.2}
\item{threshold}{Threshold to apply to the environment score
(if analysisType is 1 or 2) or to both symmetricity and environment
scores (if analysisType is 3). Examples: 0, 0.5, 1}
\item{analysisType}{Set to 1 for a general matrix where rows are observations and columns are real variables. Set to 2 for the case if the variables are boolean. Set to 3 if the variables are real and are themselves a (strict or not) subset of the observations AND you want to take in account the symmetric score. In case columns are a strict subset, you specify the mapping  
}
}
\examples{
DATA <- CLAG.loadExampleData("DIM128")
PCA <- prcomp(DATA)
RES <- CLAG.clust(DATA, delta=0.05, threshold=0.5, verbose=TRUE)
plot(PCA$x[,1], PCA$x[,2], col=rainbow(RES$nclusters)[RES$cluster])
}
\references{
CLAG: an unsupervised non hierarchical clustering algorithm handling biological data, Linda Dib, Alessandra Carbone, BMC Bioinformatics 2012, 13:194
}
\keyword{misc}
