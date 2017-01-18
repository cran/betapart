\name{functional.beta.pair}
\alias{functional.beta.pair}
\encoding{utf8}

\title{
Pair-wise functional dissimilarities
}
\description{
Computes 3 distance matrices accounting for the spatial turnover and nestedness components 
of functional beta diversity, and the sum of both values. 
Functional dissimilarities are based on volume of convex hulls intersections in a multidimensional functional space.
}
\usage{
functional.beta.pair(x, traits, index.family="sorensen")
}

\arguments{
\item{x}{ data frame, where rows are sites and columns are species. Alternatively \code{x} can 
be a \code{functional.betapart} object derived from the \code{functional.betapart.core} function}
\item{traits}{ if \code{x} is not a \code{functional.betapart} object, a data frame, where rows are species and columns are functional space dimensions (i.e. quantitative traits or synthetic axes after PCoA). Number of species in each site must be strictly higher than number of dimensions.}
\item{index.family}{ family of dissimilarity indices, partial match of \code{"sorensen"} or \code{"jaccard"}.}
}

\value{
The function returns a list with three functional dissimilarity matrices.

For \code{index.family="sorensen"} the three matrices are:
\item{funct.beta.sim}{ \code{dist} object, dissimilarity matrix accounting for functional turnover, measured as Simpson derived pair-wise functional dissimilarity}
\item{funct.beta.sne}{ \code{dist} object, dissimilarity matrix accounting for nestedness-resultant functional dissimilarity, measured as the nestedness-fraction of Sorensen derived pair-wise functional dissimilarity}
\item{funct.beta.sor}{ \code{dist} object, dissimilarity matrix accounting for functional beta diversity, measured as Sorensen derived pair-wise functional dissimilarity}


For \code{index.family="jaccard"} the three matrices are:
\item{funct.beta.jtu}{ \code{dist} object, dissimilarity matrix accounting for functional turnover, measured as the turnover-fraction of Jaccard derived pair-wise functional dissimilarity}
\item{funct.beta.jne}{ \code{dist} object, dissimilarity matrix accounting for nestedness-resultant functional dissimilarity, measured as the nestedness-fraction of Jaccard derived pair-wise functional dissimilarity}
\item{funct.beta.jac}{ \code{dist} object, dissimilarity matrix accounting for functional beta diversity, measured as Jaccard derived pair-wise functional dissimilarity}
}

\references{
Villéger S., Novack-Gottshal P. & Mouillot D. 2011. The multidimensionality of the niche reveals functional diversity changes in benthic marine biotas across geological time. Ecology Letters. 14, 561–568

Baselga, A. 2012. The relationship between species replacement, dissimilarity derived from nestedness, and nestedness. 
Global Ecology and Biogeography 21, 1223-1232

Villéger, S. Grenouillet, G., Brosse, S. 2013. Decomposing functional beta-diversity reveals that low functional 
beta-diversity is driven by low functional turnover in European fish assemblages. Global Ecology and Biogeography, 22: 671-681
}

\author{
Sébastien Villéger, Andrés Baselga and David Orme
}


\seealso{
\code{\link{functional.beta.multi}}, \code{\link{functional.betapart.core}}, \code{\link{beta.pair}}
}
\examples{
##### 4 communities in a 2D functional space (convex hulls are rectangles)
traits.test<-cbind( c(1,1,1,2,2,3,3,4,4,5,5) , c(1,2,4,1,2,3,5,1,4,3,5) )
	dimnames(traits.test)<-list(paste("sp",1:11,sep="") , c("Trait 1","Trait 2") ) 

comm.test<-matrix(0,4,11,dimnames=list( c("A","B","C","D") , paste("sp",1:11,sep="") ) )
comm.test["A",c(1,2,4,5)]<-1
comm.test["B",c(1,3,8,9)]<-1
comm.test["C",c(6,7,10,11)]<-1
comm.test["D",c(2,4,7,9)]<-1

plot(5,5,xlim=c(0,6), ylim=c(0,6), type="n", xlab="Trait 1",ylab="Trait 2")
points(traits.test[,1],traits.test[,2], pch=21,cex=1.5,bg="black")
rect(1,1,4,4, col="#458B0050", border="#458B00") ; text(2.5,2.5,"B",col="#458B00",cex=1.5)	
polygon(c(2,1,3,4), c(1,2,5,4), col="#DA70D650", border="#DA70D6") ; 
text(2.5,3,"D",col="#DA70D6",cex=1.5)	
rect(1,1,2,2, col="#FF000050" , border="#FF0000") ; text(1.5,1.5,"A",col="#FF0000",cex=1.5)	
rect(3,3,5,5, col="#1E90FF50", border="#1E90FF") ; text(4,4.2,"C",col="#1E90FF",cex=1.5)	

test.pair<-functional.beta.pair(x=comm.test, traits=traits.test, index.family = "jaccard" )
lapply(test.pair,round,2)
}