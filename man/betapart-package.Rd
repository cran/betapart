\encoding{utf8}
\name{betapart}
\alias{betapart-package}
\alias{betapart}
\docType{package}

\title{Partitioning beta diversity into turnover and nestedness components}
\description{
\pkg{betapart} allows computing pair-wise dissimilarities (distance matrices) and multiple-site dissimilarities, separating the turnover 
and nestedness-resultant components of taxonomic (incidence and abundance based), functional and phylogenetic beta diversity.
}
\details{
The partitioning of incidence-based dissimilarity can be performed for two different families of indices: 
Sorensen and Jaccard. The pairwise function \command{beta.pair} yields 3 distance matrices accounting 
for the spatial turnover and the nestedness components of beta-diversity. The third distance 
matrix accounts for the sum of both components, i.e. total dissimilarity (a monotonic transformation 
of beta diversity). 
The multiple site function \command{beta.multi} yields the spatial turnover and the nestedness components of 
overall dissimilarity, and the sum of both components, total dissimilarity. The basic calculations for all these 
multiple-site measures and pairwise dissimilarity matrices can be computed using the function \command{betapart.core}, 
which returns an object of class \code{betapart}. This is useful for large datasets as the consuming calculations 
are done only once, and its result can then be used for computing many indices. 
The multiple-site values can be randomly sampled a specified number of times for a specified number of sites using 
the function \command{beta.sample}. 
The aforementioned indices used for assessing spatial patterns can also be used for measuring temporal changes in community composition with the 
function \command{beta.temp}.
Likewise, an analogous framework has been implemented for separating the two components of abundance-based 
dissimilarity (balanced changes in abundance vs. abundance gradients) using commands \command{beta.pair.abund}, \command{beta.multi.abund}, 
\command{betapart.core.abund}, and \command{beta.sample.abund}.
The framework has been extended for functional beta diversity with commands \command{functional.betapart.core}, 
\command{functional.beta.pair} and \command{functional.beta.multi}, and for phylogenetic beta diversity with commands \command{phylo.betapart.core}, 
\command{phylo.beta.pair} and \command{phylo.beta.multi}.
The package also allows fitting negative exponential or power law distance-decay models for assessing the relationship between assemblage 
(dis)similarity and spatial (or other) distance. \command{decay.model} will fit the distance-decay function via GLM, \command{plot.decay} will 
plot the distance-decay pattern, and \command{boot.coefs.decay} will bootstrap the paramaters of the distance-decay model.
}
\author{
Andrés Baselga, David Orme, Sébastien Villéger, Julien De Bortoli and Fabien Leprieur
}

\references{
Baselga, A. 2010. Partitioning the turnover and nestedness components of beta diversity. Global Ecology and Biogeography 19:134-143

Baselga, A. 2012. The relationship between species replacement, dissimilarity derived from nestedness, and nestedness. 
Global Ecology and Biogeography 21, 1223-1232

Baselga, A. 2013. Separating the two components of abundance-based dissimilarity: balanced changes in abundance vs. 
abundance gradients. Methods in Ecology and Evolution, 4: 552-557

Baselga, A. 2017. Partitioning abundance-based multiple-site dissimilarity into components: balanced variation in abundance and 
abundance gradients. Methods in Ecology and Evolution 8: 799-808

Baselga A, Leprieur, F. 2015. Comparing methods to separate components of beta diversity. Methods in Ecology and Evolution 6: 1069-1079

Baselga A, Orme CDL. 2012. betapart: an R package for the study of beta diversity. Methods Ecol. Evol. 3: 808-812

Gómez-Rodríguez, C. & Baselga, A. 2018. Variation among European beetle taxa in patterns of distance decay 
of similarity suggests a major role of dispersal processes. Ecography, in press

Legendre P. 2014. Interpreting the replacement and richness difference components of beta diversity. Global Ecology and Biogeography, 23: 1324–1334

Leprieur F, Albouy C, De Bortoli J, Cowman PF, Belwood DR, Mouillot D. 2012. Quantifying phylogenetic beta diversity: 
distinguishing between "true" turnover of lineages and phylogenetic diversity gradients. PLoS One 7(8): e42760 

Villéger, S. Grenouillet, G., Brosse, S. 2013. Decomposing functional beta-diversity reveals that low functional 
beta-diversity is driven by low functional turnover in European fish assemblages. Global Ecology and Biogeography, 22: 671-681

}