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
\method{plot}{decay}(x, xlim=c(0,max(x$data.x)), ylim=c(0,1), 
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
of similarity suggests a major role of dispersal processes. Ecography 41: 1825-1834
}

\author{
Andrés Baselga
}

\seealso{
\code{\link{decay.model}}
}

\examples{

# presence/absence tables for longhorn beetles of South and North Europe
data(ceram.s)
data(ceram.n)

# spatial coordinates of territories in South and North Europe
data(coords.s)
data(coords.n)

# dissimilarity matrices
ceram.s.sim<-beta.pair(ceram.s)$beta.sim
ceram.n.sim<-beta.pair(ceram.n)$beta.sim

# spatial distances in km
distgeo.s<-dist(coords.s[,1:2])
distgeo.n<-dist(coords.n[,1:2])

# Negative exponential models for the decay of similarity with spatial distance
decay.south<-decay.model(y=1-ceram.s.sim, x=distgeo.s, y.type="sim", model.type="exp")
decay.north<-decay.model(y=1-ceram.n.sim, x=distgeo.n, y.type="sim", model.type="exp")

# Plot the decay models
plot.decay(decay.south, col="red")
plot.decay(decay.north, col="blue", add=TRUE)

# Equivalent models for the increase of dissimilarity with spatial distance
increase.south<-decay.model(y=ceram.s.sim, x=distgeo.s, y.type="dissim", model.type="exp")
increase.north<-decay.model(y=ceram.n.sim, x=distgeo.n, y.type="dissim", model.type="exp")

# Plot the decay models
plot.decay(increase.south, col="red")
plot.decay(increase.north, col="blue", add=TRUE)

}