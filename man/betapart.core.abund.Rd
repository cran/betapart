\encoding{utf8}
\name{betapart.core.abund}
\alias{betapart.core.abund}

\title{
Core calculations of betapart for abundance-based dissimilarity measures
}
\description{
Computes the basic quantities needed for computing the abundance-based multiple-site 
dissimilarity measures and pairwise dissimilarity matrices.
}
\usage{
betapart.core.abund(x)
}

\arguments{
\item{x}{ data frame, where rows are sites and columns are species} 
}

\value{
The function returns an object of class \code{betapart.abund} with the following elements:
\item{multiple.shared.abund}{ the multiple-site intersection component in terms of abundances (AM)}
\item{pair.shared.abund}{ a matrix containing the agreement in abundance between pairs of sites (A)}
\item{min.not.shared.abund}{ a matrix containing the minimum disagreement in abundance between pairs of sites: min(B,C)}
\item{max.not.shared.abund}{ a matrix containing the maximum disagreement in abundance between pairs of sites between pairs of sites: max(B,C)}
\item{pair.not.shared.abund}{ a matrix containing the total disagreement in abundance between pairs of sites: B+C}
}

\references{
Baselga, A. 2013. Separating the two components of abundance-based dissimilarity: balanced changes in 
abundance vs. abundance gradients. Methods in Ecology and Evolution, 4: 552–557

Legendre, P. 2014. Interpreting the replacement and richness difference components of beta diversity. 
Global Ecology and Biogeography, 23: 1324–1334

Baselga, A. 2017. Partitioning abundance-based multiple-site dissimilarity into components: balanced variation in abundance and abundance gradients. Methods in Ecology and Evolution, 8: 799-808
}

\author{
Andrés Baselga
}


\seealso{
\code{\link{beta.multi.abund}}, \code{\link{beta.pair.abund}}, \code{\link{beta.sample.abund}}, \code{\link{betapart.core}}
}
\examples{
require(vegan)
data(BCI)
core.BCI<-betapart.core.abund(BCI)
pair.BCI<-beta.pair.abund(core.BCI)
multi.BCI<-beta.multi.abund(core.BCI)
}