\name{functional.betapart.core}
\alias{functional.betapart.core}
\encoding{utf8}

\title{
Core calculations of functional dissimilarities metrics
}
\description{
Computes the basic quantities needed for computing the multiple-site functional beta diversity measures
and pairwise functional dissimilarity matrices.
}

\usage{
functional.betapart.core(x, traits, multi=TRUE, warning.time=TRUE, return.details=FALSE)
}

\arguments{

\item{x}{ data frame, where rows are sites and columns are species.}
\item{traits}{  data frame, where rows are species and columns are functional space dimensions (i.e. quantitative traits or synthetic axes after PCoA). Number of species in each site must be strictly higher than number of dimensions.}
\item{multi}{ a logical value indicating whether basic quantities for multiple-site functional beta-diversity should be computed. See \code{Details}.}
\item{warning.time}{ a logical value indicating whether computation of multiple-site dissimilarities would stop if number of dimensions exceeds 4 or if number of sites exceeds 10. If turn to \code{FALSE}, computation process can be tracked in the \code{step.fbc.txt} file, see \code{Details}.}
\item{return.details}{ a logical value indicating whether volume and coordinates of vertices of convex hulls shaping each site and their intersections in the functional space should be returned.}
	}

\value{
The function returns an object of class \code{betapart} with the following elements:
\item{sumFRi}{ the sum of the functional richness values of all sites}
\item{FRt}{ the total functional richness in the dataset}
\item{a}{ the multiple-site analog of the shared functional richness term}
\item{shared}{ a matrix containing the functional richness shared between pairs of sites}
\item{not.shared}{ a matrix containing the functional richness not shared between pairs of sites: b, c}
\item{sum.not.shared}{ a matrix containing the total functional richness not shared between pairs of sites: b+c}
\item{max.not.shared}{ a matrix containing the total maximum functional richness not shared between pairs of sites: max(b,c)}
\item{min.not.shared}{ a matrix containing the total minimum functional richness not shared between pairs of sites: min(b,c)}
\item{details}{ if \code{return.details=TRUE} a list of two lists: 
			\code{$CH} a list with a vector (\code{FRi}) of functional richness in each site (i.e. convex hull volume) and \code{coord_vertices} a list of N matrices with the coordinates of species being vertices in the D-dimensions functional space. 
			\code{$intersections} a list of 3 lists: \code{$combinations}, N-1 matrices with all combinations of 2 to N sites (numbers are rank of sites in x) ; \code{$volumes}, N-1 vectors with the volume inside the intersection between each combination of sites ; \code{$coord_vertices}, list of N-1 matrices with the coordinates of the vertices shaping these intersections (NA if no intersection).
			 }
}

\details{ For multiple-site dissimilarities metrics (N>2 sites), the volume of the union of the N convex hulls is computed using the inclusion-exclusion principle (Villéger et al., 2011). It requires to compute the volume of all the intersections between 2 to N convex hulls. Intersection between k>2 convex hulls is computed as the intersection between the two convex hulls shaping intersections between the corresponding k-1 convex hulls, e.g. V(AnBnC)=V( (AnB)n(BnC) ). For N sites, computing multiple-site dissimilarity metrics thus requires computing 2^N-(N+1) pair-wise intersections between convex hulls in a multidimensional functional space.
Computation time of	the intersection between two convex hulls increases with the number of dimensions (D) of the functional space. Therefore, to prevent from running very long computation process \code{warning.time} is set by default to stop the function if D>4 or N>10.
Computation progress can be tracked in the "step.fbc.txt" file written in the working directory. This table shows proportion of steps completed for computing convex hull volume shaping each site ("FRi") and intersections between them ("intersection_k").
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
\code{\link{functional.beta.multi}}, \code{\link{functional.beta.pair}}, \code{\link{betapart.core}}
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


test.core<-functional.betapart.core(x=comm.test, traits=traits.test, multi=TRUE, 
return.details=FALSE)
test.core

# using functional.betapart.core to get details on intersections
# when only pairwise dissimilarity is computed
test.core.pair<-functional.betapart.core(x=comm.test, traits=traits.test, multi=FALSE, 
return.details=TRUE)
test.core.pair

# for multiple dissimilarity
test.core.multi<-functional.betapart.core(x=comm.test, traits=traits.test, multi=TRUE, 
return.details=TRUE)
test.core.multi

# using core outputs to compute pairwise and multiple functional dissimilarities
functional.beta.pair(x=test.core.pair, index.family = "jaccard" )
functional.beta.multi(x=test.core.multi, index.family = "jaccard" )

}