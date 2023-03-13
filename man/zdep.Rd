\encoding{utf8}
\name{zdep}
\alias{zdep}

\title{
Assessing the differences between parameters of two distance-decay models computed with decay.model(). 
}

\description{
Takes two distance-decay models fitted with decay.model() and assess via block-site bootstrap whether the parameters of both models are equal. Two tests are conducted independently, one for the equality of the first parameters of both models, another for the equality of the second parameters (slopes) of both models. The null hypothesis is that parameters of both models are equal. 
}

\usage{
zdep(m1, m2, resamples, st.val = c(1, 0))
}

\arguments{
\item{m1}{ first distance decay model, the output of decay.model().}
\item{m2}{ second distance decay model, the output of decay.model().}
\item{resamples}{ the number of bootstrap resamples.}
\item{st.val}{ starting values for the nonlinear model.}
}

\value{
The function returns a dataframe with the tests statistics (Z.dep) and respective p values for both model parameters.
}

\references{
Martín-Devasa R, Martínez-Santalla S, Gómez-Rodríguez C, Crujeiras RM, Baselga A. 2022. Comparing distance-decay parameters: a novel test under pairwise dependence. Ecological Informatics 72: 101894 

Martínez-Santalla S, Martín-Devasa R, Gómez-Rodríguez C, Crujeiras RM, Baselga A. 2022. Assessing the 
non-linear decay of community similarity: permutation and site-block resampling significance tests. 
Journal of Biogeography 49: 968-978
}

\author{
Ramiro Martín-Devasa, Sara Martínez-Santalla, Carola Gómez-Rodríguez, Rosa M. Crujeiras, Andrés Baselga
}

\seealso{
\code{\link{decay.model}}, \code{\link{plot.decay}}, \code{\link{boot.coefs.decay}}
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

# Negative exponential distance decay models
decay.south<-decay.model(y=1-ceram.s.sim, x=distgeo.s, y.type="sim", model.type="exp")
decay.north<-decay.model(y=1-ceram.n.sim, x=distgeo.n, y.type="sim", model.type="exp")

# Plot the decay models
plot.decay(decay.south, col="red")
plot.decay(decay.north, col="blue", add=TRUE)

# Assess North-South difference between intercepts and slopes
zdep(decay.south, decay.north, resamples=1000)

}