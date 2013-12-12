\name{betapart}
\alias{betapart-package}
\alias{betapart}
\encoding{utf8}
\docType{package}

\title{Partitioning beta diversity into turnover and nestedness components}
\description{
\pkg{betapart} allows computing pair-wise dissimilarities (distance matrices) and multiple-site dissimilarities, separating the turnover and nestedness-resultant components of taxonomic (incidence and abundance based), functional and phylogenetic beta diversity.
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
This framework has been extended for functional beta diversity with commands \command{functional.betapart.core}, 
\command{functional.beta.pair} and \command{functional.beta.multi}, and for phylogenetic beta diversity with commands \command{phylo.betapart.core}, 
\command{phylo.beta.pair} and \command{phylo.beta.multi}.
Likewise, an analogous framework has been implemented for separating the two components of Bray-Curtis abundance-based 
dissimilarity (balanced changes in abundance vs. abundance gradients) using command \command{bray.part}.
}
\author{
Andrés Baselga, David Orme, Sébastien Villeger, Julien De Bortoli and Fabien Leprieur

Maintainer: Andrés Baselga <andres.baselga@usc.es>
}
\references{
Baselga, A. 2010. Partitioning the turnover and nestedness components of beta diversity. Global Ecology and Biogeography 19:134-143

Baselga, A. 2012. The relationship between species replacement, dissimilarity derived from nestedness, and nestedness. 
Global Ecology and Biogeography 21, 1223-1232

Baselga, A. 2013. Separating the two components of abundance-based dissimilarity: balanced changes in abundance vs. 
abundance gradients. Methods in Ecology and Evolution, 4: 552-557

Baselga A, Orme CDL. 2012. betapart: an R package for the study of beta diversity. Methods Ecol. Evol. 3: 808-812 

Leprieur F, Albouy C, De Bortoli J, Cowman PF, Belwood DR, Mouillot D. 2012. Quantifying phylogenetic beta diversity: 
distinguishing between "true" turnover of lineages and phylogenetic diversity gradients. PLoS One 7(8): e42760 

Villéger, S. Grenouillet, G., Brosse, S. 2013. Decomposing functional beta-diversity reveals that low functional 
beta-diversity is driven by low functional turnover in European fish assemblages. Global Ecology and Biogeography, 22: 671-681

}