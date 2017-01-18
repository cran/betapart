\name{beta.temp}
\alias{beta.temp}
\encoding{utf8}

\title{
Temporal change in community composition
}
\description{
Computes the dissimilarity for each locality between time 1 and time 2, considering 
the turnover and nestedness components of temporal change, and the sum of both values (overall change)
}
\usage{
beta.temp(x, y, index.family="sorensen")
}

\arguments{
\item{x}{ data frame for time 1, where rows are sites and columns are species}
\item{y}{ data frame for time 2, where rows are sites and columns are species. 
\code{x} and \code{y} must contain exactly the same sites and species}
\item{index.family}{ family of dissimilarity indices, partial match of \code{"sorensen"} or \code{"jaccard"}.}
}

\value{
The function returns a data frame where rows are sites and columns are pairwise dissimilarity 
values between cell composition in time 1 and time 2. For \code{index.family="sorensen"} the indices are 
beta.sim, beta.sne, and beta.sor. For \code{index.family="jaccard"} the indices are beta.jtu, beta.sne, and beta.jac. 
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
\code{\link{beta.multi}}, \code{\link{beta.pair}}, \code{\link{beta.sample}}, \code{\link{betapart.core}}, 
}
\examples{
data(bbsData)
bbs.t <- beta.temp(bbs1980, bbs2000, index.family="sor")
}