\encoding{utf8}
\name{decay.model}
\alias{decay.model}

\title{
Fitting distance decay models to pair-wise assemblage similarity 
}

\description{
Fits a nonlinear model describing (i) the decay of assemblage similarity with spatial (or any other) distance,
or, equivalently, (ii) the increase of assemblage dissimilarity with distance. Nonlinear models are fitted 
via the nls.lm function in the minpack.lm package (which uses the Levenberg-Marquardt Nonlinear Least-Squares 
Algorithm). Implemented functional forms are either the (i) negative exponential, (ii) power law, or (iii) Gompertz models. 
}

\usage{
decay.model(y, x, model.type = "exponential", y.type = "similarities", 
perm = 100, st.val = c(1, 0))
}

\arguments{
\item{y}{ \code{dist} object, either containing similarities or dissimilarities between pairs of assemblages.}
\item{x}{ \code{dist} object, containing distances (spatial or other) between pairs of assemblages.}
\item{model.type}{ functional form of the model, either negative exponential, power law, or Gompertz,
partial match of \code{"exponential"}, \code{"power"}, or \code{"gompertz"}.}
\item{y.type}{ polarity of the \code{dist} object (i.e. 1 means total similarity or total dissimilarity), 
partial match of \code{"similarities"} or \code{"dissimilarities"}.}
\item{perm}{ number of permutations to assess significance.}
\item{st.val}{ starting values for the nonlinear model.}
}

\value{
The function returns a list with: 

\item{data.y}{ original y data, \code{dist} object, either containing similarities or dissimilarities between pairs of assemblages.}
\item{data.x}{ original x data, \code{dist} object, containing distances (spatial or other) between pairs of assemblages.}
\item{model}{ the fitted nonlinear model.}
\item{model.type}{  functional form of the model, either negative exponential, power law, or Gompertz.}
\item{y.type}{  similarties or dissimilarities.}
\item{first.parameter}{  first parameter of the model. It can be either the intercept, i.e. similarity or dissimilarity at distance=0, in negative exponential or power law models, or the position parameter in Gompertz models.}
\item{second.parameter}{  slope of the model, i.e. rate at which similarity decreases with distance, or dissimilarity 
increases with distance in a negative exponential, power law or Gompertz model.}
\item{aic}{  AIC of the model.}
\item{pseudo.r.squared}{  proportion of the variation in the dependent variable that the model accounts for.}
\item{p.value}{  significance of the model, as estimated from a permutation test.}
}

\references{
Gómez-Rodríguez C, Baselga A. 2018. Variation among European beetle taxa in patterns of distance decay 
of similarity suggests a major role of dispersal processes. Ecography 41: 1825-1834

Martínez-Santalla S, Martín-Devasa R, Gómez-Rodríguez C, Crujeiras RM, Baselga A. 2022. Assessing the 
non-linear decay of community similarity: permutation and site-block resampling significance tests. 
Journal of Biogeography 49: 968-978

Martín-Devasa R, Martínez-Santalla S, Gómez-Rodríguez C, Crujeiras RM, Baselga A. 2022. Species range size shapes distance decay in community similarity. Diversity and Distributions 28: 1348-1357
}

\author{
Sara Martínez-Santalla, Ramiro Martín-Devasa, Carola Gómez-Rodríguez, Rosa M. Crujeiras, Andrés Baselga
}

\seealso{
\code{\link{beta.pair}}, \code{\link{beta.pair.abund}}, \code{\link{plot.decay}}, \code{\link{boot.coefs.decay}}, \code{\link{zdep}}
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