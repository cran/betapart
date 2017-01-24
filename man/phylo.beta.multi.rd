\encoding{utf8}
\name{phylo.beta.multi}
\alias{phylo.beta.multi}

\title{
Multiple-site phylogenetic dissimilarities
}
\description{
Computes 3 distance values accounting for the multiple-site phylogenetic turnover and nestedness components
of phylogenetic beta diversity, and the sum of both values.
Phylogenetic dissimilarities are based on Faith's phylogenetic diversity.
}
\usage{
phylo.beta.multi(x, tree, index.family="sorensen")
}

\arguments{
\item{x}{ a community matrix or data frame, where rows are sites and columns are species. Alternatively \code{x} can 
be a \code{phylo.betapart} object derived from the \code{phylo.betapart.core} function}
\item{tree}{ a phylogenetic tree of class phylo with tips names identic to species names from the community matrix.}
\item{index.family}{ family of dissimilarity indices, partial match of \code{"sorensen"} or \code{"jaccard"}.}
}

\details{
The Sorensen dissimilarity index allows computing  the PhyloSor index (Bryant et al. 2008) whereas the Jaccard dissimilarity index allows computing the UniFrac index (Lozupone & Knight 2005).
}

\value{
The function returns a list with three phylogenetic dissimilarity values.

For \code{index.family="sorensen"} the three values are:
\item{phylo.beta.sim}{ \code{dist} object, dissimilarity value accounting for phylogenetic turnover, measured as Simpson derived multiple-site phylogenetic dissimilarity}
\item{phylo.beta.sne}{ \code{dist} object, dissimilarity value accounting for nestedness-resultant phylogenetic dissimilarity, measured as the nestedness-fraction of Sorensen derived multiple-site phylogenetic dissimilarity}
\item{phylo.beta.sor}{ \code{dist} object, dissimilarity value accounting for phylogenetic beta diversity, measured as Sorensen derived multiple-site phylogenetic dissimilarity}

For \code{index.family="jaccard"} the three values are:
\item{phylo.beta.jtu}{ \code{dist} object, dissimilarity value accounting for phylogenetic turnover, measured as the turnover-fraction of Jaccard derived multiple-site phylogenetic dissimilarity}
\item{phylo.beta.jne}{ \code{dist} object, dissimilarity value accounting for nestedness-resultant phylogenetic dissimilarity, measured as the nestedness-fraction of Jaccard derived multiple-site phylogenetic dissimilarity}
\item{phylo.beta.jac}{ \code{dist} object, dissimilarity value accounting for phylogenetic beta diversity, measured as Jaccard derived multiple-site phylogenetic dissimilarity}
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
\code{\link{phylo.betapart.core}}, \code{\link{beta.multi}}
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

toy.phylobetamulti<-phylo.beta.multi(toy.comm, toy.tree, index.family="sor")
toy.betamulti<-beta.multi(toy.comm, index.family="sor")

}