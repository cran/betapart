\name{functional.beta.multi}
\alias{functional.beta.multi}
\encoding{utf8}

\title{
Multiple-site functional dissimilarities
}
\description{
Computes 3 multiple-site functional dissimilarities accounting for the spatial turnover 
and the nestedness components of functional beta diversity, and the sum of both values.
Functional dissimilarities are based on volume of convex hulls intersections in a multidimensional functional space.
}
\usage{
functional.beta.multi(x, traits, index.family="sorensen", warning.time=TRUE)
}

\arguments{
\item{x}{ data frame, where rows are sites and columns are species. Alternatively \code{x} can 
be a \code{functional.betapart} object derived from the \code{functional.betapart.core} function}
\item{traits}{ if \code{x} is not a \code{functional.betapart} object, a data frame, where rows are species and columns are functional space dimensions (i.e. quantitative traits or synthetic axes after PCoA). Number of species in each site must be strictly higher than number of dimensions. Number of dimensions should not exceed 4 and number of sites should not exceed 10. See \code{Details}.}
\item{index.family}{ family of dissimilarity indices, partial match of \code{"sorensen"} or \code{"jaccard"}.}
\item{warning.time}{ a logical value indicating whether computation of multiple-site dissimilarities would stop if number of dimensions exceeds 4 or if number of sites exceeds 10. If turn to \code{FALSE}, computation process can be tracked in the \code{step.fbc.txt} file, see \code{Details}.}
}

\details{ For multiple-site dissimilarities metrics (N>2 sites), the volume of the union of the N convex hulls is computed using the inclusion-exclusion principle (Villéger et al., 2011). It requires to compute the volume of all the intersections between 2 to N convex hulls. Intersection between k>2 convex hulls is computed as the intersection between the two convex hulls shaping intersections between the corresponding k-1 convex hulls, e.g. V(AnBnC)=V( (AnB)n(BnC) ). For N sites, computing multiple-site dissimilarity metrics thus requires computing 2^N-(N+1) pair-wise intersections between convex hulls in a multidimensional functional space.
Computation time of	the intersection between two convex hulls increases with the number of dimensions (D) of the functional space. Therefore, to prevent from running very long computation process \code{warning.time} is set by default to stop the function if D>4 or N>10.
Computation progress can be tracked in the "step.fbc.txt" file written in the working directory. This table shows proportion of steps completed for computing convex hull volume shaping each site ("FRi") and intersections between them ("intersection_k").
}

\value{
The function returns a list with the three multiple site functional dissimilarity values.

For \code{index.family="sorensen"} the three indices are:
\item{beta.SIM}{ value of the functional turnover component, measured as Simpson derived functional dissimilarity}
\item{beta.SNE}{ value of the functional nestedness component, measured as Nestedness-resultant fraction of Sorensen derived functional dissimilarity}
\item{beta.SOR}{ value of the overall functional beta diversity, measured as Sorensen derived functional dissimilarity}


For \code{index.family="jaccard"} the three indices are:
\item{beta.JTU}{ value of the functional turnover component, measured as turnover fraction of Jaccard derived functional dissimilarity}
\item{beta.JNE}{ value of the functional nestedness component, measured as Nestedness-resultant fraction of Jaccard derived functional dissimilarity}
\item{beta.JAC}{ value of the overall functional beta diversity, measured as Jaccard derived functional dissimilarity}
}

\references{
Villéger S., Novack-Gottshal P. & Mouillot D. 2011. The multidimensionality of the niche reveals functional diversity changes in benthic marine biotas across geological time. Ecology Letters. 14, 561–568

Baselga, A. 2012. The relationship between species replacement, dissimilarity derived from nestedness, and nestedness. 
Global Ecology and Biogeography 21, 1223-1232

Villéger, S. Grenouillet, G., Brosse, S. in press. Decomposing functional beta-diversity reveals that low functional beta-diversity is driven by low functional turnover in European fish assemblages. Global Ecology and Biogeography, DOI: 10.1111/geb.12021
}

\author{
Sébastien Villéger, Andrés Baselga and David Orme
}

\seealso{
\code{\link{functional.beta.pair}}, \code{\link{functional.betapart.core}}, \code{\link{beta.multi}}
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
polygon(c(2,1,3,4), c(1,2,5,4), col="#DA70D650", border="#DA70D6") ; text(2.5,3,"D",col="#DA70D6",cex=1.5)	
rect(1,1,2,2, col="#FF000050" , border="#FF0000") ; text(1.5,1.5,"A",col="#FF0000",cex=1.5)	
rect(3,3,5,5, col="#1E90FF50", border="#1E90FF") ; text(4,4.2,"C",col="#1E90FF",cex=1.5)	

test.multi<-functional.beta.multi(x=comm.test, traits=traits.test, index.family = "jaccard" )
test.multi

test.multi.ABC<-functional.beta.multi(x=comm.test[c("A","B","C"),], traits=traits.test, index.family = "jaccard" )
test.multi.ABC

test.multi.ABD<-functional.beta.multi(x=comm.test[c("A","B","D"),], traits=traits.test, index.family = "jaccard" )
test.multi.ABD
}
