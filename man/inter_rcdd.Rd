% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functional.betapart.core.pairwise.R
\name{inter_rcdd}
\alias{inter_rcdd}
\title{Internal function to compute convexhull volume}
\usage{
inter_rcdd(set1, set2, qhull.opt = "FA", conv2 = function(...) NA)
}
\arguments{
\item{set1}{A matrix of coordinates}

\item{set2}{A matrix of coordinates}

\item{qhull.opt}{Qhull options, see \url{http://www.qhull.org/html/qh-optq.htm}}

\item{conv2}{A function applyed if the convexhull function crashes}
}
\value{
A volume corresponding to the intersection of the two hypervolumes
}
\description{
Estimation of the convexhull volume of the intersection of two hypervolumes based on rcdd functions
}
\examples{
\dontrun{mat1 <- matrix(runif(30), 10)
mat2 <- matrix(runif(30), 10)
inter_rcdd(mat1, mat2)}
}
\seealso{
\code{\link{inter_geom}}
}
