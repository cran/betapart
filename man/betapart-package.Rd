\name{betapart}
\alias{betapart-package}
\alias{betapart}
\encoding{utf8}
\docType{package}

\title{Partitioning beta diversity into turnover and nestedness components}
\description{
\pkg{betapart} allows computing pair-wise dissimilarities (distance matrices),
multiple-site dissimilarities and temporal changes in comunity composition, 
separating the turnover and nestedness components of beta diversity as described in Baselga (2010). 
}
\details{
The partitioning of beta diversity can be performed for two different families of indices: 
Sorensen and Jaccard. The pairwise function \command{beta.pair} yields 3 distance matrices accounting 
for the spatial turnover and the nestedness components of beta-diversity. The third distance 
matrix accounts for the sum of both componentes, i.e. beta diversity. The multiple site function \command{beta.multi} 
yields the spatial turnover and the nestedness components of overall beta diversity, and the sum 
of both components, i.e. beta diversity. The basic calculations for all these multiple-site measures and dissimilarity 
matrices can be computed using the function \command{betapart.core}, which returns an object of class \code{betapart}. 
This is useful for large datasets as the consuming calculations are done only once, and its result can then be used 
for computing many indices. The multiple-site values can be randomly sampled a specified number of 
times for a specified number of sites using the function \command{beta.sample}. The aforementioned indices used for 
assessing spatial patterns can also be used for measuring temporal changes in community composition with the 
function \command{beta.temp}.
}
\author{
Andrés Baselga and David Orme

Maintainer: Andrés Baselga <andres.baselga@usc.es>
}
\references{
Baselga, A. 2010. Partitioning the turnover and nestedness components of beta diversity. Global Ecology and Biogeography 19:134-143

Baselga, A. 2012. The relationship between species replacement, dissimilarity derived from nestedness, and nestedness. 
Global Ecology and Biogeography, in press
}