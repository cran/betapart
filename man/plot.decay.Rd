\encoding{utf8}
\name{plot.decay}
\alias{plot.decay}

\title{
Plotting distance decay curves from models computed with decay.model() 
}

\description{
Takes the output of decay.model() and plots a distance-decay curve, either a negative exponential 
or power law function as estimated with decay.model().
}

\usage{
\method{plot}{decay}(x, xlim=c(0,max(x$data[,1])), ylim=c(0,1), 
add=FALSE, remove.dots=FALSE, col="black", pch=1, lty=1, lwd=5, cex=1, ...)

}

\arguments{
\item{x}{ the output of decay.model().}
\item{xlim}{ the range of spatial distances to be plotted, default is from 0 to the maximum distance in the data.}
\item{ylim}{ the range of assemblage similarities or dissimilarities to be plotted, default is from 0 to 1.}
\item{add}{ add to the previous plot.}
\item{remove.dots}{ remove the dots from the plot, thus retaining just the decay curve.}
\item{col}{ colour used.}
\item{pch}{ symbol used for points.}
\item{lty}{ line type.}
\item{lwd}{ line width.}
\item{cex}{ scale of text and symbols.}
\item{...}{ other parameters for plotting functions.}

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

plot(spat.dist, dissim.BCI, ylim=c(0,1), xlim=c(0, max(spat.dist)))

BCI.decay.exp<-decay.model(dissim.BCI, spat.dist, model.type="exp", perm=100)

BCI.decay.pow<-decay.model(dissim.BCI, spat.dist, model.type="pow", perm=100)

plot.decay(BCI.decay.exp, col=rgb(0,0,0,0.5))
plot.decay(BCI.decay.exp, col="red", remove.dots=TRUE, add=TRUE)
plot.decay(BCI.decay.pow, col="blue", remove.dots=TRUE, add=TRUE)
}