\name{beta.sample}
\alias{beta.sample}
\encoding{utf8}

\title{
Resampling multiple-site dissimilarity for n sites
}
\description{
Resamples the 3 multiple-site dissimilarities (turnover, nestedness-resultant fraction, and overall 
beta diversity) for a subset of sites of the original data frame.
}
\usage{
beta.sample(x, index.family="sorensen", sites=10, samples=100)
}

\arguments{
\item{x}{ data frame, where rows are sites and columns are species}
\item{index.family}{ family of dissimilarity indices, partial match of \code{"sorensen"} or \code{"jaccard"}.}
\item{sites}{ number of sites for which multiple-site dissimilarities will be computed. If not specified, default is all sites.}
\item{samples}{ number of repetitions. If not specified, default is 1.}
}


\value{
The function returns a list with a dataframe with the resampled 3 multiple-site dissimilarities 
(turnover fraction, nestedness-resultant fraction and overall dissimilarity; see \code{\link{beta.multi}}), 
a vector with the respective means and a vector with the respective standard deviation.

For \code{index.family="sorensen"}:
\item{sampled.values}{ dataframe containing beta.SIM, beta.SNE and beta.SOR for all samples}
\item{mean.values}{ vector containing the mean values of beta.SIM, beta.SNE and beta.SOR among samples}
\item{sd.values}{ vector containing the sd values of beta.SIM, beta.SNE and beta.SOR among samples}


For \code{index.family="jaccard"}:
\item{sampled.values}{ dataframe containing beta.JTU, beta.JNE and beta.JAC for all samples}
\item{mean.values}{ vector containing the mean values of beta.JTU, beta.JNE and beta.JAC among samples}
\item{sd.values}{ vector containing the sd values of beta.JTU, beta.JNE and beta.JAC among samples}
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
\code{\link{beta.multi}}, \code{\link{beta.sample}}, \code{\link{beta.temp}}
}
\examples{
# Read the data for Northern and Southern European cerambycids
data(ceram.s)
data(ceram.n)

# Resample 100 times the multiple-site dissimilarities
# for 10 countries.
beta.ceram.s<-beta.sample(ceram.s, index.family="sor", sites=10, samples=100)
beta.ceram.n<-beta.sample(ceram.n, index.family="sor", sites=10, samples=100)

# Plot the distributions of beta.SIM in Southern Europe (red) 
# and Northern Europe (blue)
hist(beta.ceram.s$sampled.values$beta.SIM, breaks=10, col="red", xlim=c(0,1))
hist(beta.ceram.n$sampled.values$beta.SIM, breaks=10, col="blue", add=TRUE)
}

