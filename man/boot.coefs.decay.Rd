\encoding{utf8}
\name{boot.coefs.decay}
\alias{boot.coefs.decay}

\title{
Bootstrapping the parameters of distance-decay models computed with decay.model() 
}

\description{
Takes the output of decay.model() and bootstraps the parameters (intercept and slope).
}

\usage{
boot.coefs.decay(x, R)
}

\arguments{
\item{x}{ the output of decay.model().}
\item{R}{ the number of bootstrap samples.}
}

\value{
The function returns a list with: 

\item{model.type}{  functional form of the model, either negative exponential or power law.}
\item{y.type}{  similarties or dissimilarities.}
\item{boot.coefs}{  a matrix with the coefficients bootstrapped distributions, including intercept
 values in the first column, and slope values in the second column.}
\item{original.coefs}{  model coefficients as estimated with a GLM using decay.model().}
\item{mean.boot}{  the mean of the bootstrapped distributions.}
\item{sd.boot}{  the standard deviation of the bootstrapped distributions.}
}

\references{
Gómez-Rodríguez, C. & Baselga, A. 2018. Variation among European beetle taxa in patterns of distance decay 
of similarity suggests a major role of dispersal processes. Ecography, in press
}

\author{
Andrés Baselga
}

\seealso{
\code{\link{decay.model}}
}

\examples{
require(vegan)

data(BCI)
## UTM Coordinates (in metres)
UTM.EW <- rep(seq(625754, 626654, by=100), each=5)
UTM.NS <- rep(seq(1011569,  1011969, by=100), len=50)

spat.dist<-dist(data.frame(UTM.EW, UTM.NS))

dissim.BCI<-beta.pair.abund(BCI)$beta.bray.bal

BCI.decay.pow<-decay.model(dissim.BCI, spat.dist, model.type="pow", perm=100)

boot.coefs.decay(BCI.decay.pow, 100)
}