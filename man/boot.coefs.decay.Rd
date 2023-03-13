\encoding{utf8}
\name{boot.coefs.decay}
\alias{boot.coefs.decay}

\title{
Bootstrapping the parameters of distance-decay models computed with decay.model() 
}

\description{
Takes the output of decay.model() and bootstraps the parameters of the model (i.e. intercept and slope in negative exponential or power law models, or position parameter and slope in Gompertz models).
}

\usage{
boot.coefs.decay(m1, resamples, st.val = c(1, 0))
}

\arguments{
\item{m1}{ the output of decay.model().}
\item{resamples}{ the number of bootstrap resamples.}
\item{st.val}{ starting values for the nonlinear model.}
}

\value{
The function returns a list with: 

\item{model.type}{  functional form of the model, either negative exponential or power law.}
\item{y.type}{  similarties or dissimilarities.}
\item{boot.coefs}{  a matrix with the coefficients bootstrapped distributions, including 
 values of the first parameter (intercept or position parameter) in the first column, and values of the second parameter (slope) in the second column.}
\item{original.coefs}{  model coefficients as estimated with a nonlinear model using decay.model().}
\item{mean.boot}{  the mean of the bootstrapped distributions.}
\item{sd.boot}{  the standard deviation of the bootstrapped distributions.}
}

\references{
Gómez-Rodríguez, C. & Baselga, A. 2018. Variation among European beetle taxa in patterns of distance decay 
of similarity suggests a major role of dispersal processes. Ecography 41: 1825-1834

Martínez-Santalla S, Martín-Devasa R, Gómez-Rodríguez C, Crujeiras RM, Baselga A. 2022. Assessing the 
non-linear decay of community similarity: permutation and site-block resampling significance tests. 
Journal of Biogeography 49: 968-978
}

\author{
Sara Martínez-Santalla, Ramiro Martín-Devasa, Carola Gómez-Rodríguez, Rosa M. Crujeiras, Andrés Baselga
}

\seealso{
\code{\link{decay.model}}, \code{\link{zdep}}
}

\examples{

# presence/absence tables for longhorn beetles of South Europe
data(ceram.s)

# spatial coordinates of territories in South Europe
data(coords.s)

# dissimilarity matrix
ceram.s.sim<-beta.pair(ceram.s)$beta.sim

# spatial distances in km
distgeo.s<-dist(coords.s[,1:2])

# Negative exponential model for the decay of similarity with spatial distance
decay.south<-decay.model(y=1-ceram.s.sim, x=distgeo.s, y.type="sim", model.type="exp")

# Site-block bootstrap
boot.coefs.decay(decay.south, resamples=100)

}