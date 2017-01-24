\encoding{utf8}
\name{betapart.core}
\alias{betapart.core}

\title{
Core calculations of betapart
}
\description{
Computes the basic quantities needed for computing the multiple-site beta diversity measures
and pairwise dissimilarity matrices.
}
\usage{
betapart.core(x)
}

\arguments{
\item{x}{ data frame, where rows are sites and columns are species} 
}

\value{
The function returns an object of class \code{betapart} with the following elements:
\item{sumSi}{ the sum of the species richness values of all sites}
\item{St}{ the total richness in the dataset}
\item{a}{ the multiple-site analog of the shared species term}
\item{shared}{ a matrix containing the number of species shared between pairs of sites}
\item{not.shared}{ a matrix containing the number of species not shared between pairs of sites: b, c}
\item{sum.not.shared}{ a matrix containing the total number of species not shared between pairs of sites: b+c}
\item{max.not.shared}{ a matrix containing the total maximum number of species not shared between pairs of sites: max(b,c)}
\item{min.not.shared}{ a matrix containing the total minimum number of species not shared between pairs of sites: min(b,c)}
}

\references{
Baselga, A. 2010. Partitioning the turnover and nestedness components of beta diversity. Global Ecology and Biogeography 19:134-143

Baselga, A. 2012. The relationship between species replacement, dissimilarity derived from nestedness, and nestedness. 
Global Ecology and Biogeography 21, 1223-1232
}

\author{
Andrés Baselga and David Orme
}


\seealso{
\code{\link{beta.multi}}, \code{\link{beta.pair}}, \code{\link{beta.sample}}, \code{\link{beta.temp}},
}
\examples{
data(ceram.s)
ceram.core.s<-betapart.core(ceram.s)
ceram.dist.jac<-beta.pair(ceram.core.s, index.family="jac")
ceram.dist.sor<-beta.pair(ceram.core.s, index.family="sor")
ceram.multi.jac<-beta.multi(ceram.core.s, index.family="jac")
ceram.multi.sor<-beta.multi(ceram.core.s, index.family="sor")

}