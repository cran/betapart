\name{beta.multi}
\alias{beta.multi}
\encoding{utf8}

\title{
Multiple-site dissimilarities
}
\description{
Computes 3 multiple-site dissimilarities accounting for the spatial turnover 
and the nestedness components of beta diversity, and the sum of both values
}
\usage{
beta.multi(x, index.family="sorensen")
}

\arguments{
\item{x}{ data frame, where rows are sites and columns are species. Alternatively \code{x} can 
be a \code{betapart} object derived from the \code{betapart.core} function}
\item{index.family}{ family of dissimilarity indices, partial match of \code{"sorensen"} or \code{"jaccard"}.}
}

\value{
The function returns a list with the three multiple site dissimilarity values.

For \code{index.family="sorensen"} the three indices are:
\item{beta.SIM}{ value of the turnover component, measured as Simpson dissimilarity}
\item{beta.SNE}{ value of the nestedness component, measured as Nestedness-resultant fraction of Sorensen dissimilarity}
\item{beta.SOR}{ value of the overall beta diversity, measured as Sorensen dissimilarity}


For \code{index.family="jaccard"} the three indices are:
\item{beta.JTU}{ value of the turnover component, measured as turnover fraction of Jaccard dissimilarity}
\item{beta.JNE}{ value of the nestedness component, measured as Nestedness-resultant fraction of Jaccard dissimilarity}
\item{beta.JAC}{ value of the overall beta diversity, measured as Jaccard dissimilarity}
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
\code{\link{beta.pair}}, \code{\link{beta.sample}}, \code{\link{betapart.core}}, \code{\link{beta.temp}}
}
\examples{
data(ceram.s)
ceram.beta<-beta.multi(ceram.s, index.family="sor")
}



