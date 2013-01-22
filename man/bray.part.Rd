\name{bray.part}
\alias{bray.part}
\encoding{utf8}

\title{
Partitioning pair-wise Bray-Curtis dissimilarities
}
\description{
Computes 3 distance matrices accounting for the balanced variation and abundance gradient components 
of Bray-Curtis dissimilarity, and the sum of both values (i.e. Bray-Curtis dissimilarity)
}
\usage{
bray.part(x)
}

\arguments{
\item{x}{ data frame of species abundances, where rows are sites and columns are species.}
}

\value{
The function returns a list with three dissimilarity matrices.

\item{bray.bal}{ \code{dist} object, dissimilarity matrix accounting for the dissimilarity derived from 
balanced variation in abundance between sites}
\item{bray.gra}{ \code{dist} object, dissimilarity matrix accounting for the dissimilarity derived from
unidirectional abundance gradients}
\item{bray}{ \code{dist} object, dissimilarity matrix accounting for total abundance-based dissimilarity 
between sites, measured as the Bray-Curtis index}
}

\references{
Baselga, A. in press. Separating the two components of abundance-based dissimilarity: balanced changes in 
abundance vs. abundance gradients. Methods in Ecology and Evolution. DOI: 10.1111/2041-210X.12029
}


\author{
Andrés Baselga
}


\seealso{
\code{\link{beta.pair}}
}
\examples{
require(vegan)
data(BCI)
BCI.matrices<-bray.part(BCI)
}

