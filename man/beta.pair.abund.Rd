\name{beta.pair.abund}
\alias{beta.pair.abund}
\encoding{utf8}

\title{
Abundance-based pair-wise dissimilarities
}
\description{
Computes 3 distance matrices accounting for the (i) balanced variation in abundances,
(ii) abundance gradients, and (iii) total dissimilarity (i.e. the sum of both components).
}
\usage{
beta.pair.abund(x, index.family = "bray")
}
\arguments{
  \item{x}{
data frame, where rows are sites and columns are species. Alternatively \code{x} can be a \code{betapart.abund} object derived from the \code{betapart.core.abund} function
}
  \item{index.family}{
family of dissimilarity indices, partial match of \code{"bray"} or \code{"ruzicka"}.
}
}
\value{
The function returns a list with three dissimilarity matrices.
For \code{index.family="bray"} the three matrices are:
\item{beta.bray.bal}{ \code{dist} object, dissimilarity matrix accounting for the dissimilarity derived from balanced variation in abundance between sites}
\item{beta.bray.gra}{ \code{dist} object, dissimilarity matrix accounting for the dissimilarity derived from unidirectional abundance gradients}
\item{beta.bray}{ \code{dist} object, dissimilarity matrix accounting for total abundance-based dissimilarity between sites, measured as the Bray-Curtis index}

For \code{index.family="ruzicka"} the three matrices are:
\item{beta.ruz.bal}{ \code{dist} object, dissimilarity matrix accounting for the dissimilarity derived from balanced variation in abundance between sites}
\item{beta.ruz.gra}{ \code{dist} object, dissimilarity matrix accounting for the dissimilarity derived from unidirectional abundance gradients}
\item{beta.ruz}{ \code{dist} object, dissimilarity matrix accounting for total abundance-based dissimilarity between sites, measured as the Ruzicka index}
}
\references{
Baselga, A. 2013. Separating the two components of abundance-based dissimilarity: balanced changes in abundance vs. abundance gradients. Methods in Ecology and Evolution 4: 552–557

Legendre, P. 2014. Interpreting the replacement and richness difference components of beta diversity. Global Ecology and Biogeography, 23: 1324–1334
}
\author{
Andrés Baselga
}


\seealso{
\code{\link{beta.multi.abund}}, \code{\link{beta.sample.abund}}, \code{\link{betapart.core.abund}}, \code{\link{beta.pair}}
}
\examples{
require(vegan)
data(BCI)
BCI.pair<-beta.pair.abund(BCI, index.family="bray")
}

