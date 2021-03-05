\name{functional.betapart.core}
\alias{functional.betapart.core}
\encoding{utf8}

\title{
Core calculations of functional dissimilarities metrics
}
\description{
Computes the basic quantities needed for computing the multiple-site functional beta diversity measures
and pairwise functional dissimilarity matrices. This version of the function now supports internal parallelization to fasten the computations and external parallelization for null models. 
}

\usage{
functional.betapart.core(x, traits, multi = TRUE, warning.time = TRUE, 
                         return.details = FALSE, fbc.step = FALSE, 
                         parallel = FALSE, opt.parallel = beta.para.control(),
                         convhull.opt = qhull.opt(),
                         progress = FALSE)
}

\arguments{

\item{x}{ data frame, where rows are sites and columns are species.}
\item{traits}{  data frame, where rows are species and columns are functional space dimensions (i.e. quantitative traits or synthetic axes after PCoA). Number of species in each site must be strictly higher than number of dimensions.}
\item{multi}{ a logical value indicating whether basic quantities for multiple-site functional beta-diversity should be computed. See \code{Details}.}
\item{warning.time}{ a logical value indicating whether computation of multiple-site dissimilarities would stop if number of dimensions exceeds 4 or if number of sites exceeds 10. If turn to \code{FALSE}, computation process can be tracked in the \code{step.fbc.txt} file, see \code{Details}.}
\item{return.details}{ a logical value indicating whether volume and coordinates of vertices of convex hulls shaping each site and their intersections in the functional space should be returned.}
\item{fbc.step}{ a logical value indicating whether the computation progress tracking file "step.fbc.txt" should be created; Setting it to \code{FALSE} will speed up the function. It is automatically turn to \code{FALSE} when parallel is \code{TRUE}.}
\item{parallel}{ a logical value indicating if internal parallelization is used to compute pairwise dissimilarities, see \code{Examples}. If \code{multi} is set to \code{TRUE} parallelization will be turned-off.}
\item{opt.parallel}{ a list of values to replace the default values returned by the function \code{\link{beta.para.control}} to customize the cluster used for parallel computing.}
\item{convhull.opt}{ a list of values to replace the default values returned by the function \code{\link{qhull.opt}}, to pass options to the \code{\link{convhulln}} function. The first element, \code{conv1}, set the options that be used by default with \code{convhulln}, while the second element, \code{conv2}, set the options that will be used if \code{convhulln} returns an error. By default (\code{conv1 = 'QJ'} and \code{conv2 = NULL}). If \code{convhulln} generates an error, a \code{NA} is returned. This prevents the function to stop when encountering an error.}
\item{progress}{ a logical indicating if a progress bar should be displayed (\code{TRUE}) or not (by default).}
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
			\code{$intersections} a list of 3 lists: \code{$combinations}, N-1 matrices with all combinations of 2 to N sites (numbers are rank of sites in x) ; \code{$volumes}, N-1 vectors with the volume inside the intersection between each combination of sites ; \code{$coord_vertices}, list of N-1 matrices with the coordinates of the vertices shaping these intersections (NA if no intersection). See \url{http://www.qhull.org/html/qh-optq.htm} for the possible options.
			 }
}

\details{ For multiple-site dissimilarities metrics (N>2 sites), the volume of the union of the N convex hulls is computed using the inclusion-exclusion principle (Villéger et al., 2011). It requires to compute the volume of all the intersections between 2 to N convex hulls. Intersection between k>2 convex hulls is computed as the intersection between the two convex hulls shaping intersections between the corresponding k-1 convex hulls, e.g. V(AnBnC)=V( (AnB)n(BnC) ). For N sites, computing multiple-site dissimilarity metrics thus requires computing 2^N-(N+1) pair-wise intersections between convex hulls in a multidimensional functional space.
Computation time of	the intersection between two convex hulls increases with the number of dimensions (D) of the functional space. Therefore, to prevent from running very long computation process \code{warning.time} is set by default to stop the function if D>4 or N>10.

If \code{fbc.step} is set to \code{TRUE}, computation progress can be tracked in the "step.fbc.txt" file written in the working directory. This table shows proportion of steps completed for computing convex hull volume shaping each site ("FRi") and intersections between them ("intersection_k"). This is only possible when computations are not performed in parallel, and this whatever the type of parallelization used (external or internal).  

If \code{parallel} is set to \code{TRUE}, computation will be run though the creation of a cluster. This is interesting when beta diversity computation is long. When the number of sites increase and/or when the taxonomic richness is highly variable between sites, parallelization becomes more and more interesting. On small matrices, the running time could inflate due to the creation of the cluster and its management. 
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
traits.test <- cbind(c(1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5), 
                     c(1, 2, 4, 1, 2, 3, 5, 1, 4, 3, 5))
dimnames(traits.test) <- list(paste("sp", 1:11, sep=""), c("Trait 1", "Trait 2")) 

comm.test <- matrix(0, 4, 11, dimnames = list(c("A", "B", "C", "D"), 
                    paste("sp", 1:11, sep="")))
comm.test["A", c(1, 2, 4, 5)] <- 1
comm.test["B", c(1, 3, 8, 9)] <- 1
comm.test["C", c(6, 7, 10, 11)] <- 1
comm.test["D", c(2, 4, 7, 9)] <- 1

plot(5, 5, xlim = c(0, 6), ylim = c(0, 6), type = "n", xlab = "Trait 1", ylab = "Trait 2")
points(traits.test[, 1], traits.test[, 2], pch = 21, cex = 1.5, bg = "black")
rect(1, 1, 4, 4, col = "#458B0050", border = "#458B00")
text(2.5, 2.5, "B" , col = "#458B00", cex = 1.5)	
polygon(c(2, 1, 3, 4), c(1, 2, 5, 4), col = "#DA70D650", border = "#DA70D6") 
text(2.5, 3, "D", col = "#DA70D6", cex = 1.5)	
rect(1, 1, 2, 2, col = "#FF000050", border = "#FF0000")
text(1.5, 1.5, "A", col = "#FF0000", cex = 1.5)	
rect(3, 3, 5, 5, col = "#1E90FF50", border = "#1E90FF")
text(4, 4.2, "C", col = "#1E90FF", cex = 1.5)	

# for multiple dissimilarity, multi = TRUE
test.core <- functional.betapart.core(x = comm.test, traits = traits.test, 
                                      multi = TRUE, return.details = FALSE)
test.core

# for pairwise dissimilarity, multi = FALSE
test.core <- functional.betapart.core(x = comm.test, traits = traits.test, 
                                      multi = FALSE, return.details = FALSE)
test.core

# to use systematilcally the "QJ" options 
test.core <- functional.betapart.core(x = comm.test, traits = traits.test, 
                                      multi = FALSE, return.details = FALSE,
                                      convhull.opt = list(conv1 = "QJ"))

# to use the "QJ" options only if the convhull function generates an error
# instead of returning NA
test.core <- functional.betapart.core(x = comm.test, traits = traits.test, 
                                      multi = FALSE, return.details = FALSE,
                                      convhull.opt = list(conv2 = "QJ"))

# using functional.betapart.core to get details on intersections
# when only pairwise dissimilarity is computed
test.core.pair <- functional.betapart.core(x = comm.test, traits = traits.test, 
                                           multi = FALSE, return.details = TRUE)
test.core.pair

# for multiple dissimilarity
test.core.multi <- functional.betapart.core(x = comm.test, traits = traits.test, 
                                            multi = TRUE, return.details = TRUE)
test.core.multi


# using core outputs to compute pairwise and multiple functional dissimilarities
functional.beta.pair(x = test.core.pair, index.family = "jaccard" )
functional.beta.multi(x = test.core.multi, index.family = "jaccard" )


##### using internal parallelisation to fasten pairiwse dissimilarity
# by default it use serial computation
test.core.pair <- functional.betapart.core(x = comm.test, traits = traits.test, 
                                           multi = FALSE, return.details = FALSE, 
                                           fbc.step = FALSE, parallel = FALSE)
                                          
# by default it uses half of the cores and 1 task per run (this can be customised)
# test.core.pairp <- functional.betapart.core(x = comm.test, traits = traits.test, 
#                                            multi = FALSE, return.details = FALSE, 
#                                            fbc.step = FALSE, parallel = TRUE)

# you can set the number of core to use :
test.core.pairp <- functional.betapart.core(x = comm.test, traits = traits.test, 
                                            multi = FALSE, return.details = FALSE, 
                                            fbc.step = FALSE, parallel = TRUE, 
                                            opt.parallel = beta.para.control(nc = 2))
all.equal(test.core.pair, test.core.pairp)

# library(microbenchmark)
# microbenchmark(serial = 
#                  functional.betapart.core(comm.test, traits.test, multi = FALSE, 
#                                           return.details = FALSE, fbc.step = FALSE,
#                                           parallel = FALSE),
#                nc2 = 
#                  functional.betapart.core(comm.test, traits.test, multi = FALSE, 
#                                           return.details = FALSE, fbc.step = FALSE,
#                                           parallel = TRUE, 
#                                           opt.parallel = beta.para.control(nc = 2)),
#                nc4 =
#                  functional.betapart.core(comm.test, traits.test, multi = FALSE, 
#                                           return.details = FALSE, fbc.step = FALSE,
#                                           parallel = TRUE, 
#                                           opt.parallel = beta.para.control(nc = 4))
# )

\dontrun{
# If the number of species is very different among communities
# load-balancing parallelisation could be very efficient
# especially when the number of community is high
test.core.pairp <- functional.betapart.core(comm.test, traits.test, multi = FALSE, 
                                            return.details = FALSE, fbc.step = FALSE,
                                            parallel = TRUE, 
                                            opt.parallel = 
                                              beta.para.control(nc = 2, LB = TRUE))

# use you can use fork cluster (but not on Windows)
#test.core.pairp <- functional.betapart.core(comm.test, traits.test, multi = FALSE, 
#                                            return.details = FALSE, fbc.step = FALSE,
#                                            parallel = TRUE, 
#                                            opt.parallel = 
#                                              beta.para.control(nc = 2, type = "FORK"))

# finally you can customise the number of task run at each time
test.core.pairp <- functional.betapart.core(comm.test, traits.test, multi = FALSE, 
                                            return.details = FALSE, fbc.step = FALSE,
                                            parallel = TRUE, 
                                            opt.parallel = 
                                              beta.para.control(nc = 2, size = 6))

# using internal parallelisation is not always useful, especially on small data set
# load balancing is very helpful when species richness are highly variable

# Null model using 'external' parallel computing 

# Example 1: pairwise functional beta diversity (functional.beta.pair)
# Note that this is an example with a small number of samples and null model 
# permutations for illustration.
# Real null model analyses should have a much greater number of samples and permutations.

##### 4 communities in a 3D functional space

comm.test <- matrix(0, 4, 11, dimnames = list(c("A", "B", "C", "D"), 
                                              paste("sp", 1:11, sep = "")))
comm.test["A", c(1, 2, 4, 5)] <- 1
comm.test["B", c(1, 3, 8, 9)] <- 1
comm.test["C", c(6, 7, 10, 11)] <- 1
comm.test["D", c( 2, 4, 7, 9)] <- 1

set.seed(1)
traits.test <- matrix(rnorm(11*3, mean = 0, sd = 1), 11, 3) 
dimnames(traits.test) <- list(paste("sp", 1:11, sep = "") , 
                              c("Trait 1", "Trait 2", "Trait 3"))

# Required packages
library(doSNOW)
library(picante)
library(fastmatch)
library(foreach)


# define number of cores
# Use parallel::detectCores() to determine number of cores available in your machine
nc <- 2 

# 4 cores would be better (nc <- 4)

# create cluster
cl <- snow::makeCluster(nc)

# register parallel backend
doSNOW::registerDoSNOW(cl)

# define number of permutations for the null model (the usual is 1000)
# make sure that nperm/nc is a whole number so that all cores have the same number 
# of permutations to work on
nperm <- 100

test.score <- functional.betapart.core(comm.test, traits.test, multi = FALSE, 
                                       warning.time = FALSE, return.details = FALSE, 
                                       fbc.step = FALSE, parallel = FALSE)

obs.pair.func.dis <- functional.beta.pair(x = test.score, index.family = "sorensen")

# transform functional.beta.pair results into a matrix
obs.pair.func.dis <- do.call(rbind, obs.pair.func.dis)

# set names for each pair of site
pair_names <- combn(rownames(comm.test), 2, FUN = paste, collapse = "_")
colnames(obs.pair.func.dis) <- pair_names

# export necessary variables and functions to the cluster of cores
snow::clusterExport(cl = cl, c("comm.test", "traits.test"),
                        envir = environment())

# creation of an iterator to run 1 comparaisons on each core at time
it <- itertools::isplitIndices(nperm, chunkSize = 1)

# parallel computation
null.pair.func.dis <- 
  foreach(n = it, .combine = c, .packages=c("picante","betapart","fastmatch")) \%dopar\% {
                                
          # it enables to adjust the number of permutations (nt) done on each run
          nt <- length(n)
          null.pair.temp <- vector(mode = "list", length = nt)
                                
          # for each core "n" perform "nt" permutations
          for (j in 1:nt){ 
                    
            # randomize community with chosen null model
            # for this particular example we used the "independent swap algorithm" 
            # but the user can choose other types of permutation
            # or create it's own null model 
            null.comm.test <- randomizeMatrix(comm.test, null.model = "independentswap", 
                                              iterations=1000)
                                  
            # execute functional.betapart.core function 
            null.test.score <- 
              try(functional.betapart.core(null.comm.test, traits = traits.test, 
                                           multi = FALSE, warning.time = FALSE, 
                                           return.details = FALSE,  fbc.step = FALSE, 
                                           parallel = FALSE), silent = TRUE)
            # using 'external' parallelisation it is necessary to set parralel to FALSE
                                  
            # in this artificial example there are a few combinations of species that 
            # the convex hull cannot be calculated due to some odd geometric combination 
            # so we need to re-permute the community matrix 
            while(inherits(null.test.score, "try-error")){
                                    
              null.comm.test <- randomizeMatrix(comm.test, null.model = "independentswap", 
                                                iterations = 1000)
              null.test.score <- 
                try(functional.betapart.core(x = null.comm.test, traits = traits.test, 
                                             multi = FALSE, warning.time = FALSE, 
                                             return.details = FALSE, fbc.step = FALSE, 
                                             parallel = FALSE), silent = TRUE)
              }
                                  
            # compute the pairwise beta-diversity null values and input them in the 
            # temporary result matrix
            res <- functional.beta.pair(x = null.test.score, index.family = "sorensen")
            null.pair.temp[[j]] <- do.call(rbind, res)
                                  
        }
        #retrieve the results from each core
        null.pair.temp
    }

# stop cluster
snow::stopCluster(cl)


# compute the mean, standard deviation and p-values of dissimilarity metrics
# for each pair of site

mean.null.pair.func <- matrix(numeric(), ncol = ncol(obs.pair.func.dis), 
                              nrow = nrow(obs.pair.func.dis))
sd.null.pair.func <- matrix(numeric(), ncol = ncol(obs.pair.func.dis), 
                            nrow = nrow(obs.pair.func.dis))
p.pair.func.dis <- matrix(numeric(), ncol = ncol(obs.pair.func.dis), 
                          nrow = nrow(obs.pair.func.dis))

# for each one of the 3 null dissimilarity metrics (SIN, SNE and SOR) 
for (j in 1:nrow(obs.pair.func.dis)){
  matnull <- sapply(null.pair.func.dis, function(x) x[j,])
  mean.null.pair.func[j,] <- rowMeans(matnull)
  sd.null.pair.func[j,] <- sqrt(rowSums((matnull - mean.null.pair.func[j,])^2)/(nperm-1))
  p.pair.func.dis[j,] <- rowSums(matnull >= obs.pair.func.dis[j,]) 
  p.pair.func.dis[j,] <- (pmin(p.pair.func.dis[j,],nperm-p.pair.func.dis[j,])+1)/(nperm+1)
  # the +1 is to take into account that the observed value is one of the possibilities
}

# compute standardized effect sizes
ses.pair.func.dis <- (obs.pair.func.dis - mean.null.pair.func)/sd.null.pair.func


# Example 2: multiple functional beta diversity (functional.beta.multi)
# Note that this is an example with a small number of samples and null model
# permutations for illustration.
# Real null model analyses should have a much greater number of samples 
# and permutations.

##### 4 communities in a 3D functional space

comm.test <- matrix(0, 4, 11,dimnames = list(c("A", "B", "C", "D"), 
                                             paste("sp", 1:11, sep = "")))
comm.test["A", c(1, 2, 4, 5)] <- 1
comm.test["B", c(1, 3, 8, 9)] <- 1
comm.test["C", c(6, 7, 10, 11)] <- 1
comm.test["D", c(2, 4, 7, 9)]<-1

set.seed(1)
traits.test <- matrix(rnorm(11*3, mean=0, sd=1), 11, 3) 
dimnames(traits.test) <- list(paste("sp", 1:11, sep=""), 
                              c("Trait 1", "Trait 2", "Trait 3"))

# Required packages
library(doSNOW)
library(picante)
library(fastmatch)
library(foreach)


# define number of cores
# Use parallel::detectCores() to determine number of cores available in your machine
nc <- 2 

# create cluster
cl <- snow::makeCluster(nc)

# register parallel backend 
doSNOW::registerDoSNOW(cl)

# define number of permutations for the null model (the usual is 1000)
# make sure that nperm/nc is a whole number so that all cores have the same number 
# of permutations to work on
nperm <- 10

# compute observed values for multiple functional dissimilarities
test.score <- functional.betapart.core(comm.test, traits.test, multi = TRUE, 
                                       warning.time = FALSE, return.details = FALSE, 
                                       fbc.step = FALSE, parallel = FALSE)
obs.multi.func.dis <- do.call(cbind, functional.beta.multi(test.score, 
                                                           index.family = "sorensen"))

# export necessary variables and functions to the cluster of cores
snow::clusterExport(cl = cl, c("comm.test", "traits.test"), 
                        envir=environment())

it <- itertools::isplitIndices(nperm, chunkSize = 1)
  
null.multi.func.dis <- 
  foreach(n = it, .combine = rbind, 
          .packages = c("picante","betapart","fastmatch")) \%dopar\% {
                                 
    # for each core, create temporary matrix to store 3 null multiple functional 
    # dissimilarity indices (SIN, SNE,SOR)
    null.multi.temp <- matrix(numeric(), ncol = 3, nrow = length(n), 
                              dimnames = list(NULL, c("funct.beta.SIM", "funct.beta.SNE", 
                                              "funct.beta.SOR")))
                                 
    # number of tasks per core (i.e., permutations per core)
    nt <- length(n)
                                 
    # for each core "n" perform "nt" permutations
    for (j in 1:nt) { 
                                     
      # randomize community matrix with chosen null model (for this example 
      # we chose the "independent swap" algorithm) 
      null.comm.test <- randomizeMatrix(comm.test, null.model="independentswap", 
                                        iterations=1000)
                                   
      # execute functional.betapart.core function identifying each "n" core 
      # with the core.ident argument for external parallelization, 
      null.test.score <- 
        try(functional.betapart.core(null.comm.test, traits.test, multi = TRUE,  
                                     warning.time = FALSE, return.details = FALSE, 
                                     fbc.step = FALSE, parallel = FALSE), 
            silent = TRUE)
                                   
      # in this artificial example there are a few combinations of species
      # that the convex hull 
      # cannot be calculated due to some odd geometric combination so we 
      # need to re-permute the community matrix 
                                   
      while(inherits(null.test.score, "try-error")){
        null.comm.test <- randomizeMatrix(comm.test, null.model="independentswap", 
                                          iterations=1000)
        null.test.score <- 
          try(functional.betapart.core(null.comm.test, traits.test, multi = TRUE,
                                       warning.time = FALSE, return.details = FALSE, 
                                       fbc.step = FALSE, parallel = FALSE), 
              silent = TRUE)
        }    
      # input null values in the temporary result matrix
      null.multi.temp[j,] <- unlist(functional.beta.multi(null.test.score, 
                                                          index.family = "sorensen"))
      }
    # recover results from each core
    null.multi.temp
    }

# close cluster
snow::stopCluster(cl)

# result matrix
result <- matrix(numeric(), ncol = 3, nrow = 3, 
                 dimnames = list(c("obs", "ses", "p"), colnames(obs.multi.func.dis)))

# input observed values for the multiple functional dissimilarity indices (SIN, SNE,SOR)
result[1,] = obs.multi.func.dis

# compute standardized effect sizes (ses) for the multiple functional
# dissimilarity indices (SIN, SNE,SOR)
result[2,] <- (obs.multi.func.dis-colMeans(null.multi.func.dis, na.rm=TRUE))/
  apply(null.multi.func.dis,2, sd, na.rm=TRUE)

# compute p-values for the multiple functional dissimilarity indices (SIN, SNE,SOR)
for (i in 1:3) {
  result[3, i] <- sum(obs.multi.func.dis[i]<=null.multi.func.dis[,i])
  result[3, i] <- (pmin(result[3, i], nperm - result[3, i]) + 1)/(nperm+1)
}
# the +1 is to take into account that the observed value is one of the possibilities

result
###     
}

}