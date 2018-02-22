\encoding{utf8}
\name{decay.model}
\alias{decay.model}

\title{
Fitting distance decay models to pair-wise assemblage similarity 
}

\description{
Fits a negative exponential or power law function (via GLM) describing (i) the decay of assemblage similarity 
with spatial (or any other) distance, or, equivalently, (ii) the increase of assemblage dissimilarity with distance.
}

\usage{
decay.model(y, x, model.type="exponential", y.type="similarities", perm=100)
}

\arguments{
\item{y}{ code{dist} object, either containing similarities or dissimilarities between pairs of assemblages.}
\item{x}{ code{dist} object, containing distances (spatial or other) between pairs of assemblages.}
\item{model.type}{ functional form of the model, either negative exponential or power law, 
partial match of \code{"exponential"} or \code{"power"}.}
\item{y.type}{ polarity of the code{dist} object (i.e. 1 means total similarity or total dissimilarity), 
partial match of \code{"similarities"} or \code{"dissimilarities"}.}
\item{perm}{ number of randomizations to assess significance.}
}

\value{
The function returns a list with: 

\item{data}{ dataframe containing distances (spatial or other) and similarities (or dissimilarities).}
\item{model}{ the fitted GLM.}
\item{model.type}{  functional form of the model, either negative exponential or power law.}
\item{y.type}{  similarties or dissimilarities.}
\item{pseudo.r.squared}{  similarties or dissimilarities.}
\item{a.intercept}{  intercept of the model, i.e. similarity or dissimilarity at distance=0.}
\item{b.slope}{  slope of the model, i.e. rate at which similarity decreases with distance, or dissimilarity 
increases with distance in a negative exponential or power law model.}
\item{p.value}{  significance of the model, as estimated from a randomization test.}
}

\references{
Nekola, J.C. & McGill, B.J. 2014. Scale dependency in the functional form of the distance decay 
relationship. Ecography, 37, 309-320

Gómez-Rodríguez, C. & Baselga, A. 2018. Variation among European beetle taxa in patterns of distance decay 
of similarity suggests a major role of dispersal processes. Ecography, in press
}

\author{
Andrés Baselga
}

\seealso{
\code{\link{beta.pair}}, \code{\link{beta.pair.abund}}
}

\examples{
require(vegan)

data(BCI)
## UTM Coordinates (in metres)
UTM.EW <- rep(seq(625754, 626654, by=100), each=5)
UTM.NS <- rep(seq(1011569,  1011969, by=100), len=50)

spat.dist<-dist(data.frame(UTM.EW, UTM.NS))

dissim.BCI<-beta.pair.abund(BCI)$beta.bray.bal

plot(spat.dist, dissim.BCI, ylim=c(0,1), xlim=c(0, max(spat.dist)))

BCI.decay.exp<-decay.model(dissim.BCI, spat.dist, model.type="exp", perm=100)

BCI.decay.pow<-decay.model(dissim.BCI, spat.dist, model.type="pow", perm=100)

plot.decay(BCI.decay.exp, col=rgb(0,0,0,0.5))
plot.decay(BCI.decay.exp, col="red", remove.dots=TRUE, add=TRUE)
plot.decay(BCI.decay.pow, col="blue", remove.dots=TRUE, add=TRUE)
}