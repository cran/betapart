\encoding{utf8}
\name{beta.sample.abund}
\alias{beta.sample.abund}

\title{
Resampling abundance-based multiple-site dissimilarity for n sites
}
\description{
Resamples the 3 abundance-based multiple-site dissimilarities (balanced variation fraction,abundance-gradient fraction, and overall 
dissimilarity) for a subset of sites of the original data frame.
}
\usage{
beta.sample.abund(x, index.family="bray", sites = nrow(x), samples = 1)
}

\arguments{
\item{x}{ data frame, where rows are sites and columns are species}
\item{index.family}{ family of dissimilarity indices, partial match of \code{"bray"} or \code{"ruzicka"}.}
\item{sites}{ number of sites for which multiple-site dissimilarities will be computed. If not specified, default is all sites.}
\item{samples}{ number of repetitions. If not specified, default is 1.}
}

\value{
The function returns a list with a dataframe with the resampled 3 multiple-site dissimilarities 
(balanced variation fraction, abundance-gradient fraction and overall dissimilarity; see \code{\link{beta.multi.abund}}), 
a vector with the respective means and a vector with the respective standard deviation.

For \code{index.family="bray"}:
\item{sampled.values}{ dataframe containing beta.BRAY.BAL, beta.BRAY.GRA and beta.BRAY for all samples}
\item{mean.values}{ vector containing the mean values of beta.BRAY.BAL, beta.BRAY.GRA and beta.BRAY among samples}
\item{sd.values}{ vector containing the sd values of beta.BRAY.BAL, beta.BRAY.GRA and beta.BRAY among samples}


For \code{index.family="ruzicka"}:
\item{sampled.values}{ dataframe containing beta.RUZ.BAL, beta.RUZ.GRA and beta.RUZ for all samples}
\item{mean.values}{ vector containing the mean values of beta.RUZ.BAL, beta.RUZ.GRA and beta.RUZ among samples}
\item{sd.values}{ vector containing the sd values of beta.RUZ.BAL, beta.RUZ.GRA and beta.RUZ among samples}
}

\references{
Baselga, A. 2017. Partitioning abundance-based multiple-site dissimilarity into components: balanced variation in abundance and abundance gradients. Methods in Ecology and Evolution 8: 799-808
}

\author{
Andrés Baselga
}

\seealso{
\code{\link{beta.multi.abund}}, \code{\link{beta.sample}}
}
\examples{
require(vegan)
data(BCI)
beta.sample.abund(BCI, index.family="bray", sites=10, samples=100)
}