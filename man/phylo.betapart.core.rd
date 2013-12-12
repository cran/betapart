\name{phylo.betapart.core}
\alias{phylo.betapart.core}
\encoding{utf8}

\title{
Core calculations of phylogenetic dissimilarities metrics
}
\description{
Computes the basic quantities needed for computing the multiple-site phylogenetic beta diversity measures
and pairwise phylogenetic dissimilarity matrices.
}
\usage{
phylo.betapart.core(x, tree)
}

\arguments{
\item{x}{ a community matrix or data frame, where rows are sites and columns are species.}
\item{tree}{ a phylogenetic tree of class phylo with tips names identic to species names from the community matrix.}
}

\value{
The function returns a list with:
\item{sumSi}{ the sum of the phylogenetic diversity values of all sites}
\item{St}{ the total phylogenetic diversity in the dataset}
\item{shared}{ a matrix containing the phylogenetic diversity shared between pairs of sites}
\item{sum.not.shared}{ a matrix containing the total phylogenetic diversity not shared between pairs of sites: b+c}
\item{max.not.shared}{ a matrix containing the total maximum phylogenetic diversity not shared between pairs of sites: max(b,c)}
\item{min.not.shared}{ a matrix containing the total minimum phylogenetic diversity not shared between pairs of sites: min(b,c)}
}

\references{
Baselga A. (2012) The relationship between species replacement, dissimilarity derived from nestedness, and nestedness.
Global Ecology and Biogeography 21, 1223-1232

Bryant JA, Lamanna C, Morlon H, Kerkhoff AJ, Enquist BJ, et al. (2008) Microbes on mountainsides: Contrasting elevational patterns of bacterial and plant diversity. Proceedings of the National Academy of Sciences of the United States of America 105: 11505-11511.

Faith DP, Lozupone CA, Nipperess D, Knight R (2009) The Cladistic Basis for the Phylogenetic Diversity (PD) Measure Links Evolutionary Features to Environmental Gradients and Supports Broad Applications of Microbial Ecology's "Phylogenetic Beta Diversity" Framework. Int J Mol Sci 10: 4723-4741. doi: 10.3390/ijms10114723.

Leprieur F, Albouy C, De Bortoli J, Cowman PF, Bellwood DR, et al. (2012) Quantifying Phylogenetic Beta Diversity: Distinguishing between "True" Turnover of Lineages and Phylogenetic Diversity Gradients. PLoS ONE 7(8): e42760. doi:10.1371/journal.pone.0042760

Lozupone C, Knight R (2005) UniFrac: a new phylogenetic method for comparing microbial communities. Applied and Environmental Microbiology 71: 8228-8235.
}

\author{
Julien De Bortoli (juldebortoli@yahoo.fr), Fabien Leprieur(fabien.leprieur@univ-montp2.fr), Andrés Baselga and David Orme
}


\seealso{
\code{\link{phylo.beta.pair}}, \code{\link{phylo.beta.multi}}
}
\examples{

# toy tree for 6 species (sp1 to sp6)
require(ape)
toy.tree<-read.tree(text="(((sp1:1,sp2:1):5,(sp3:3,sp4:3):3):2,(sp5:7,sp6:7):1);")
plot(toy.tree)

# toy community table with 6 assemblages (A to F) with 6 species (sp1 to sp6)
toy.comm<-matrix(nrow=6, ncol=6)
rownames(toy.comm)<-c("A","B","C","D","E","F")
colnames(toy.comm)<-c("sp1","sp2","sp3","sp4","sp5","sp6")
toy.comm[1,]<-c(1,1,1,0,0,0)
toy.comm[2,]<-c(0,1,1,1,0,0)
toy.comm[3,]<-c(0,0,1,1,1,0)
toy.comm[4,]<-c(0,0,1,1,1,1)
toy.comm[5,]<-c(0,0,0,1,1,1)
toy.comm[6,]<-c(1,0,0,1,1,1)

toy.phylocore<-phylo.betapart.core(toy.comm, toy.tree)
}
