\encoding{utf8}
\name{beta.multi.abund}
\alias{beta.multi.abund}

\title{
Abundance-based multiple-site dissimilarities
}
\description{
Computes 3 multiple-site dissimilarities accounting for the (i) balanced variation 
and (ii) abundance gradient components of dissimilarity, and the sum of both values (i.e. total abundance-based dissimilarity)
}
\usage{
beta.multi.abund(x, index.family="bray")
}

\arguments{
\item{x}{ data frame, where rows are sites and columns are species. Alternatively \code{x} can 
be a \code{betapart.abund} object derived from the \code{betapart.core.abund} function}
\item{index.family}{ family of dissimilarity indices, partial match of \code{"bray"} or \code{"ruzicka"}.}
}

\value{
The function returns a list with the three multiple site dissimilarity values.

For \code{index.family="bray"} the three indices are:
\item{beta.BRAY.BAL}{ value of the balanced variation component of Bray-Curtis multiple-site dissimilarity}
\item{beta.BRAY.GRA}{ value of the abundance-gradient component of Bray-Curtis multiple-site dissimilarity}
\item{beta.BRAY}{ value of the overall dissimilarity, measured as Bray-Curtis multiple-site dissimilarity}


For \code{index.family="ruzicka"} the three indices are:
\item{beta.RUZ.BAL}{ value of the balanced variation component of Ruzicka multiple-site dissimilarity}
\item{beta.RUZ.GRA}{ value of the abundance-gradient component of Ruzicka multiple-site dissimilarity}
\item{beta.RUZ}{ value of the overall dissimilarity, measured as Ruzicka multiple-site dissimilarity}
}

\references{
Baselga, A. 2017. Partitioning abundance-based multiple-site dissimilarity into components: balanced variation in abundance and abundance gradients. Methods in Ecology and Evolution 8: 799-808
}

\author{
Andrés Baselga
}


\seealso{
\code{\link{beta.pair.abund}}, \code{\link{beta.sample.abund}}, \code{\link{betapart.core.abund}}, \code{\link{beta.multi}}
}
\examples{
require(vegan)
data(BCI)
beta.multi.abund(BCI, index.family="bray")
}