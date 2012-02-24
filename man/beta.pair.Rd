\name{beta.pair}
\alias{beta.pair}
\encoding{utf8}

\title{
Pair-wise dissimilarities
}
\description{
Computes 3 distance matrices accounting for the spatial turnover and nestedness components 
of beta diversity, and the sum of both values
}
\usage{
beta.pair(x, index.family="sorensen")
}

\arguments{
\item{x}{ data frame, where rows are sites and columns are species. Alternatively \code{x} can 
be a \code{betapart} object derived from the \code{betapart.core} function}
\item{index.family}{ family of dissimilarity indices, partial match of \code{"sorensen"} or \code{"jaccard"}.}
}

\value{
The function returns a list with three dissimilarity matrices.

For \code{index.family="sorensen"} the three matrices are:
\item{beta.sim}{ \code{dist} object, dissimilarity matrix accounting for spatial turnover, measured as Simpson pair-wise dissimilarity}
\item{beta.sne}{ \code{dist} object, dissimilarity matrix accounting for nestedness-resultant dissimilarity, measured as the nestedness-fraction of Sorensen pair-wise dissimilarity}
\item{beta.sor}{ \code{dist} object, dissimilarity matrix accounting for beta diversity, measured as Sorensen pair-wise dissimilarity}


For \code{index.family="jaccard"} the three matrices are:
\item{beta.jtu}{ \code{dist} object, dissimilarity matrix accounting for spatial turnover, measured as the turnover-fraction of Jaccard pair-wise dissimilarity}
\item{beta.jne}{ \code{dist} object, dissimilarity matrix accounting for nestedness-resultant dissimilarity, measured as the nestedness-fraction of Jaccard pair-wise dissimilarity}
\item{beta.jac}{ \code{dist} object, dissimilarity matrix accounting for beta diversity, measured as Jaccard pair-wise dissimilarity}
}

\references{
Baselga, A. 2010. Partitioning the turnover and nestedness components of beta diversity. Global Ecology and Biogeography 19:134-143

Baselga, A. 2012. The relationship between species replacement, dissimilarity derived from nestedness, and nestedness. 
Global Ecology and Biogeography, in press
}

\author{
Andrés Baselga and David Orme
}


\seealso{
\code{\link{beta.multi}}, \code{\link{beta.sample}}, \code{\link{betapart.core}}, \code{\link{beta.temp}}
}
\examples{
data(ceram.s)
ceram.dist<-beta.pair(ceram.s, index.family="jac")
}

