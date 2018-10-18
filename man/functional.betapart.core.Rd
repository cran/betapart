\name{functional.betapart.core}
\alias{functional.betapart.core}
\encoding{utf8}

\title{
Core calculations of functional dissimilarities metrics
}
\description{
Computes the basic quantities needed for computing the multiple-site functional beta diversity measures
and pairwise functional dissimilarity matrices. This version of the function now supports external parallelization for null models. 
}

\usage{
functional.betapart.core(x, traits, multi=TRUE, warning.time=TRUE, return.details=FALSE,  
fbc.step=FALSE, core.ident=NULL) 
}

\arguments{

\item{x}{ data frame, where rows are sites and columns are species.}
\item{traits}{  data frame, where rows are species and columns are functional space dimensions (i.e. quantitative traits or synthetic axes after PCoA). Number of species in each site must be strictly higher than number of dimensions.}
\item{multi}{ a logical value indicating whether basic quantities for multiple-site functional beta-diversity should be computed. See \code{Details}.}
\item{warning.time}{ a logical value indicating whether computation of multiple-site dissimilarities would stop if number of dimensions exceeds 4 or if number of sites exceeds 10. If turn to \code{FALSE}, computation process can be tracked in the \code{step.fbc.txt} file, see \code{Details}.}
\item{return.details}{ a logical value indicating whether volume and coordinates of vertices of convex hulls shaping each site and their intersections in the functional space should be returned.}
\item{fbc.step}{ a logical value indicating whether the computation progress tracking file "step.fbc.txt" should be created; Setting it to \code{FALSE} will speed up the function.}
\item{core.ident}{ computer core identifier to allow for external parallelization of this function (e.g., for null models),see \code{Examples}.}
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

If \code{fbc.step} is set to \code{TRUE}, computation progress can be tracked in the "step.fbc.txt" file written in the working directory. In the case of external parallelization, each core will write a unique computation progress tracking file "step.fbc_core.ident.txt" using the \code{core.ident} argument. This table shows proportion of steps completed for computing convex hull volume shaping each site ("FRi") and intersections between them ("intersection_k"). 
}

\references{
Villéger S., Novack-Gottshal P. & Mouillot D. 2011. The multidimensionality of the niche reveals functional diversity changes in benthic marine biotas across geological time. Ecology Letters. 14, 561-568

Baselga, A. 2012. The relationship between species replacement, dissimilarity derived from nestedness, and nestedness. 
Global Ecology and Biogeography 21, 1223-1232

Villéger, S. Grenouillet, G., Brosse, S. 2012. Decomposing functional beta-diversity reveals that low functional beta-diversity is driven by low functional turnover in European fish assemblages. Global Ecology and Biogeography, in press
}

\author{
Sébastien Villéger, Andrés Baselga, David Orme, Renato Henriques-Silva, Maxime Logez 
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
polygon(c(2,1,3,4), c(1,2,5,4), col="#DA70D650", border="#DA70D6") 
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


# Null model using parallel computing 

# Example 1: pairwise functional beta diversity (functional.beta.pair)
# Note that this is an example with a small number of samples and null model 
# permutations for illustration.
# Real null model analyses should have a much greater number of samples and permutations.

##### 4 communities in a 3D functional space

comm.test<-matrix(0,4,11,dimnames=list( c("A","B","C","D") , paste("sp", 1:11, sep = "")))
comm.test["A",c(1,2,4,5)]<-1
comm.test["B",c(1,3,8,9)]<-1
comm.test["C",c(6,7,10,11)]<-1
comm.test["D",c(2,4,7,9)]<-1

set.seed(1)
traits.test<-matrix(rnorm(11*3, mean=0, sd=1), 11, 3) 
dimnames(traits.test)<-list(paste("sp", 1:11, sep="") , c("Trait 1","Trait 2", "Trait 3"))
     
# Required packages
require(doParallel)
require(picante)
require(fastmatch)
require(foreach)

    
# define number of cores
# Use parallel::detectCores() to determine number of cores available in your machine
nc <- 2 

# create cluster
cl <- parallel::makeCluster(nc)

# register parallel backend
registerDoParallel(cl)
    
# define number of permutations for the null model (the usual is 1000)
# make sure that nperm/nc is a whole number so that all cores have the same number 
# of permutations to work on
nperm <- 10
    
# compute observed values for multiple functional dissimilarities
test.score <- functional.betapart.core(x = comm.test, traits = traits.test, 
multi = FALSE, warning.time = FALSE, return.details = FALSE, 
fbc.step = FALSE, core.ident = NULL)
                                           
obs.pair.func.dis <- functional.beta.pair(x = test.score, index.family = "sorensen")

# number of possible pairs of sites
pairs <- combn(nrow(comm.test), 2)

# create strings for each pair of site
pair_names <- NULL
for (i in 1:ncol(pairs)){
  pair_id <- paste(attr(obs.pair.func.dis$funct.beta.sim, "Labels")[pairs[1,i]],
  attr(obs.pair.func.dis$funct.beta.sim, "Labels")[pairs[2,i]], sep="-")
  pair_names<-c(pair_names, pair_id) 
}

# transform functional.beta.pair results into a matrix
obs.pair.func.dis <- do.call(rbind, lapply(obs.pair.func.dis, function(x) c(x)))

# set names for each pair of site using the set of strings created previsouly
colnames(obs.pair.func.dis) <- pair_names

# export necessary variables and functions to the cluster of cores
parallel::clusterExport(cl = cl, c("nperm", "nc", "comm.test", "traits.test"), 
envir = environment())
                                     
null.pair.func.dis <- foreach::foreach(n = 1:nc, .combine = c, 
.packages=c("picante","betapart","fastmatch")) \%dopar\% {
    
# number of tasks per core (i.e., permutations per core)
nt <- nperm/nc
     
#create a list of lists where the results from "nt" permutations will be stored
  null.pair.temp <- replicate(nt,list())
  
  # for each core "n" perform "nt" permutations
  for (j in 1:nt){ 
    
    # randomize community with chosen null model
    # for this particular example we used the "independent swap algorithm" 
    # but the user can choose other types of permutation, or create it's own null model 
    null.comm.test <- randomizeMatrix(comm.test, null.model="independentswap", 
	iterations=1000)
    
    # execute functional.betapart.core function identifying each "n" core with the 
	  # core.ident argument 
    null.test.score <- try(functional.betapart.core(null.comm.test, 
	  traits = traits.test, multi = FALSE, warning.time = FALSE, return.details = FALSE, 
	  fbc.step = FALSE, core.ident = n), silent=TRUE)
    
    # in this artificial example there are a few combinations of species that the convex hull 
	  # cannot be calculated due to 
    # some odd geometric combination so we need to re-permute the community matrix 
    while(inherits(null.test.score, "try-error")){
      
    null.comm.test <- randomizeMatrix(comm.test, null.model="independentswap", 
	iterations=1000)
    null.test.score <- try(functional.betapart.core(x = null.comm.test, 
	traits = traits.test, multi = FALSE, 
    warning.time = FALSE, return.details = FALSE, fbc.step = FALSE, 
	core.ident = n), silent=TRUE)
    }
    
    # compute the pairwise beta-diversity null values and input them in the temporary 
	# result matrix
    null.pair.temp[[j]] <- functional.beta.pair(x = null.test.score, 
	index.family = "sorensen" )
    
  }
  #retrieve the results from each core
  null.pair.temp
}

# stop cluster
parallel::stopCluster(cl)

# transform the results from each permutation into a matrix
null.pair.func.dis <- lapply(null.pair.func.dis, function(x) do.call(rbind,
lapply(x, function(x) c(x))))

# compute the mean, standard deviation and p-values of dissimilary metrics
# for each pair of site

mean.null.pair.func <-matrix(numeric(),ncol = ncol(obs.pair.func.dis), 
nrow = nrow(obs.pair.func.dis))
sd.null.pair.func <-matrix(numeric(),ncol = ncol(obs.pair.func.dis), 
nrow = nrow(obs.pair.func.dis))
p.pair.func.dis<-matrix(numeric(),ncol = ncol(obs.pair.func.dis), 
nrow = nrow(obs.pair.func.dis))

# for each one of the 3 null dissimilarity metrics (SIN, SNE and SOR) 
for (j in 1:nrow(obs.pair.func.dis)){
  
  #for each pair of sites
  for (i in 1:ncol(obs.pair.func.dis)){
    
# group the null values from each metric of each pair of site into into
#	it's own matrix
null.pair.temp <- do.call(rbind, lapply(null.pair.func.dis, `[`,j,i)) 
    
# compute mean of all null values
mean.null.pair.func[j,i] <- mean(null.pair.temp)
    
# compute standard deviation of all null values
sd.null.pair.func[j,i] <- sd(null.pair.temp)
    
# compute the p-values
p.pair.func.dis[j,i] <-(length(which(obs.pair.func.dis[j,i]<null.pair.temp))+1)/nperm
# the +1 is to take into account that the observed value is one of the possibilities
  }
}

# compute standardized effect sizes
ses.pair.func.dis <- (obs.pair.func.dis - mean.null.pair.func)/sd.null.pair.func

# result matrix
result <- as.data.frame(matrix(numeric(), ncol = 9,nrow = ncol(obs.pair.func.dis),
     dimnames = list(pair_names, c("obs.SIM","ses.SIM","p.SIM","obs.SNE",
     "ses.SNE","p.SNE","obs.SOR","ses.SOR","p.SOR"))))

result$obs.SIM <- obs.pair.func.dis[1,]
result$obs.SNE <- obs.pair.func.dis[2,]
result$obs.SOR <- obs.pair.func.dis[3,]

result$ses.SIM <- ses.pair.func.dis[1,]
result$ses.SNE <- ses.pair.func.dis[2,]
result$ses.SOR <- ses.pair.func.dis[3,]

result$p.SIM <- p.pair.func.dis[1,]
result$p.SNE <- p.pair.func.dis[2,]
result$p.SOR <- p.pair.func.dis[3,]

result

\dontrun{
# Example 2: multiple functional beta diversity (functional.beta.multi)
# Note that this is an example with a small number of samples and null model
# permutations for illustration.
# Real null model analyses should have a much greater number of samples 
# and permutations.

##### 4 communities in a 3D functional space

comm.test<-matrix(0,4,11,dimnames=list( c("A","B","C","D"), 
paste("sp", 1:11, sep = "")))
comm.test["A",c(1,2,4,5)]<-1
comm.test["B",c(1,3,8,9)]<-1
comm.test["C",c(6,7,10,11)]<-1
comm.test["D",c(2,4,7,9)]<-1

set.seed(1)
traits.test<-matrix(rnorm(11*3, mean=0, sd=1), 11, 3) 
dimnames(traits.test)<-list(paste("sp", 1:11, sep=""), 
c("Trait 1","Trait 2", "Trait 3"))

# Required packages
require(doParallel)
require(picante)
require(fastmatch)
require(foreach)


# define number of cores
# Use parallel::detectCores() to determine number of cores available in your machine
nc <- 2 
    
# create cluster
cl <- parallel::makeCluster(nc)
    
# register parallel backend 
registerDoParallel(cl)
    
# define number of permutations for the null model (the usual is 1000)
# make sure that nperm/nc is a whole number so that all cores have the same number 
# of permutations to work on
nperm <- 10
    
# compute observed values for multiple functional dissimilarities
test.score <- functional.betapart.core(x = comm.test, traits = traits.test, 
multi = TRUE, warning.time = FALSE, return.details = FALSE, fbc.step = FALSE, 
core.ident = NULL)
obs.multi.func.dis <- t(as.matrix(unlist(functional.beta.multi(x=test.score, 
index.family = "sorensen"))))
    
# export necessary variables and functions to the cluster of cores
parallel::clusterExport(cl = cl, c("nperm", "nc", "comm.test", "traits.test"), 
envir=environment())
    
null.multi.func.dis <- foreach::foreach(n = 1:nc, .combine = rbind, 
.packages = c("picante","betapart","fastmatch")) \%dopar\% {
  
 # for each core, create temporary matrix to store 3 null multiple functional 
 # dissimilarity indices (SIN, SNE,SOR)
 null.multi.temp <- matrix(numeric(), ncol = 3, nrow = nperm/nc, 
 dimnames = list(NULL, c("funct.beta.SIM", "funct.beta.SNE", "funct.beta.SOR")))

 # number of tasks per core (i.e., permutations per core)
 nt <- nperm/nc

      # for each core "n" perform "nt" permutations
      for (j in 1:nt) { 
  
      # randomize community matrix with chosen null model (for this example 
	  # we chose the "independent swap" algorithm) 
      null.comm.test <- picante::randomizeMatrix(comm.test, 
	  null.model="independentswap", iterations=1000)
    
      # execute functional.betapart.core function identifying each "n" core 
      # with the core.ident argument for external parallelization, 
      null.test.score <- try(functional.betapart.core.modif(x = null.comm.test, 
	  traits = traits.test, 
      multi = TRUE,  warning.time = FALSE, return.details = FALSE, fbc.step = FALSE, 
	  core.ident = n), silent = TRUE)

      # in this artificial example there are a few combinations of species
	  # that the convex hull 
      # cannot be calculated due to some odd geometric combination so we 
	  # need to re-permute the community matrix 
                        
      while(inherits(null.test.score, "try-error")){
      null.comm.test <- randomizeMatrix(comm.test, null.model="independentswap", 
	  iterations=1000)
      null.test.score <- try(functional.betapart.core(x = null.comm.test, 
	  traits = traits.test, 
	  multi = TRUE,  warning.time = FALSE, return.details = FALSE, fbc.step = FALSE, 
	  core.ident = n), silent = TRUE)
                        }    
      # input null values in the temporary result matrix
      null.multi.temp[j,] <- unlist(functional.beta.multi(x = null.test.score, 
	  index.family = "sorensen"))
      }
# recover results from each core
null.multi.temp
}

# close cluster
parallel::stopCluster(cl)

# result matrix
result <- matrix(numeric(), ncol = 3, nrow = 3, dimnames = list(c("obs","ses","p"), 
colnames(obs.multi.func.dis)))

# input observed values for the multiple functional dissimilarity indices (SIN, SNE,SOR)
result[1,] = obs.multi.func.dis

# compute standardized effect sizes (ses) for the multiple functional
# dissimilarity indices (SIN, SNE,SOR)
result[2,] <- (obs.multi.func.dis-colMeans(null.multi.func.dis, na.rm=TRUE))/
apply(null.multi.func.dis,2,sd, na.rm=TRUE)

# compute p-values for the multiple functional dissimilarity indices (SIN, SNE,SOR)
# SIM
result[3,1] <- (length(which(obs.multi.func.dis[1]<null.multi.func.dis[,1])) + 1)/nperm
# SNE
result[3,2] <- (length(which(obs.multi.func.dis[2]<null.multi.func.dis[,2])) + 1)/nperm
# SOR 
result[3,3] <- (length(which(obs.multi.func.dis[3]<null.multi.func.dis[,3])) + 1)/nperm
# the +1 is to take into account that the observed value is one of the possibilities
     
result
###     
}
}