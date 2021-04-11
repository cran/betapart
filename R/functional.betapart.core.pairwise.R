#' functional.betapart.core.pairwise
#' @description Computes the basic quantities needed for computing the pairwise functional dissimilarity matrices. 
#' This function is similar to functional.betapart.core with multi=FALSE but it provides more options for computing convex hulls shaping each assemblage (through option passed to qhull algorithm in 'geometry::convhulln()' as well as their intersections (computed based on library 'geometry' whenever possible, else with 'rcdd' as in functional.betapart.core.

#'
#' @usage functional.betapart.core.pairwise(x, traits, 
#'                                   return.details = TRUE,
#'                                   parallel = FALSE, 
#'                                   opt.parallel = beta.para.control(), 
#'                                   convhull.opt = qhull.opt(),
#'                                   progress = FALSE)
#'  
#' @param x A data frame, where rows are sites and columns are species.
#' @param traits A data frame, where rows are species and columns are functional space 
#' dimensions (i.e. quantitative traits or synthetic axes after PCoA). Number of species
#'  in each site must be strictly higher than number of dimensions.
#' @param return.details A logical value indicating if informations concerning the 
#' computation of the convexhull volumes should be returned (\code{TRUE}) or not (\code{FALSE}). See \code{Details}.
#' @param parallel A logical value indicating if internal parallelization is used to compute pairwise dissimilarities (\code{TRUE}), see Examples..
#' @param opt.parallel A list of four values to modify default values used to define and run the parallel cluster. 
#' See \code{\link{beta.para.control}} and the \code{Details} section.
#' @param convhull.opt A list of two named vectors of character \code{conv1} and \code{conv2} defining the options that 
#' will be passed to the \code{\link{convhulln}} function to compute the convexhull volumes (\url{http://www.qhull.org/html/qh-optq.htm}).
#'  \code{conv1} sets the default option that will be passed tp \code{convhulln} (\code{'QJ'} by default instead of \code{'Qt'}), 
#'  while \code{conv2} set the options if \code{convhulln} return an error with options set by \code{conv1}. 
#'  See \code{Details} and the examples.
#' @param progress A logical indicating if a progress bar should be displayed (\code{TRUE}) or not (by default) (\code{FALSE}).
#' 
#' @details {\itemize{\item{opt.parallel:}{ Among the four options (see \code{\link{beta.para.control}}),
#' the number of cores (\code{nc}) and the number of convexhull volumes computed by each core at each iteration {\code{size}}
#' are the most important with this function. As \code{\link{inter_rcdd}} is very fast, it is necessary to set a large value (>100) for \code{size}.
#' Otherwise the parallelisation would not be so efficient. With a low number of communities using internal parallelisation will 
#' slow down the function.}
#' \item{convhull.opt:}{ Some specific distribution of points could generate errors when computing the convexhull volumes with
#' the default qhull options ('Qt'), that is why \code{'QJ'} was prefered as default values (\code{conv1}). Sometimes, it could be 
#' interesting to use alternative options such as 'Qt' or 'Qs' or to use alternative options to achieve the computation. 
#' These alternative options could be used to compute all the convexhull volumes with \code{conv1} or only when an error occurs
#' using \code{conv2}. It is thus possible to define \code{conv1} as \code{'Qt'} to use the default 'qhull' options but to prevent error 
#' by setting \code{'QJ'} to the \code{conv2} argument.}
#' }
#' }
#' 
#' @return {The function returns an object of class \code{betapart} with the following elements:
#'   \describe{
#'   \item{sumFRi}{ The sum of the functional richness values of all sites}
#'   \item{FRt}{ The total functional richness in the dataset \code{NA}. Kept for compatibility with \code{functional.betapart.core}}
#'   \item{a}{ The multiple-site analog of the shared functional richness term, \code{NA}. Kept for compatibility with \code{functional.betapart.core} }
#'   \item{shared}{ A matrix containing the functional richness shared between pairs of sites}
#'   \item{not.shared}{ A matrix containing the functional richness not shared between pairs of sites: b, c}
#'   \item{sum.not.shared}{ A matrix containing the total functional richness not shared between pairs of sites: b+c}
#'   \item{max.not.shared}{ A matrix containing the total maximum functional richness not shared between pairs of sites: max(b,c)}
#'   \item{min.not.shared}{ A matrix containing the total minimum functional richness not shared between pairs of sites: min(b,c)}
#'   \item{details}{ \code{NA} if \code{return.details = FALSE}. Otherwise a list of two elements: 
#'   \itemize{
#'   \item{\code{$FRi}}{ a data frame with two columns, the \code{FRi} values and the qhull options used to compute them (\code{qhull.opt}).}
#'   \item{\code{$Intersection}}{ a data frame with the pairs of communities (\code{Comms}), 
#'   the function used to compute the volume of their intersections (\code{Inter}) and the qhull options used (\code{qhull.opt})}
#'   }
#'   }
#'   }
#' }
#' 
#' @examples
#' ##### 4 communities in a 2D functional space (convex hulls are rectangles)
#' traits.test <- cbind(c(1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5), 
#'                     c(1, 2, 4, 1, 2, 3, 5, 1, 4, 3, 5))
#' dimnames(traits.test) <- list(paste("sp", 1:11, sep=""), c("Trait 1", "Trait 2")) 
#'
#' comm.test <- matrix(0, 4, 11, dimnames = list(c("A", "B", "C", "D"), 
#'                                              paste("sp", 1:11, sep="")))
#' comm.test["A", c(1, 2, 4, 5)] <- 1
#' comm.test["B", c(1, 3, 8, 9)] <- 1
#' comm.test["C", c(6, 7, 10, 11)] <- 1
#' comm.test["D", c(2, 4, 7, 9)] <- 1
#' 
#' plot(5, 5, xlim = c(0, 6), ylim = c(0, 6), type = "n", xlab = "Trait 1", ylab = "Trait 2")
#' points(traits.test[, 1], traits.test[, 2], pch = 21, cex = 1.5, bg = "black")
#' rect(1, 1, 4, 4, col = "#458B0050", border = "#458B00")
#' text(2.5, 2.5, "B" , col = "#458B00", cex = 1.5)	
#' polygon(c(2, 1, 3, 4), c(1, 2, 5, 4), col = "#DA70D650", border = "#DA70D6") 
#' text(2.5, 3, "D", col = "#DA70D6", cex = 1.5)	
#' rect(1, 1, 2, 2, col = "#FF000050", border = "#FF0000")
#' text(1.5, 1.5, "A", col = "#FF0000", cex = 1.5)	
#' rect(3, 3, 5, 5, col = "#1E90FF50", border = "#1E90FF")
#' text(4, 4.2, "C", col = "#1E90FF", cex = 1.5)	
#' 
#' # for pairwise dissimilarity
#' test.core <- functional.betapart.core.pairwise(x = comm.test, traits = traits.test,
#'                                                return.details = FALSE)
#' test.core
#' # equivalent to 
#' test <- functional.betapart.core(x = comm.test, traits = traits.test,
#'                                  return.details = FALSE,
#'                                  multi = FALSE)
#' all.equal(test.core, test)
#' 
#' # using core outputs to compute pairwise and multiple functional dissimilarities
#' functional.beta.pair(x = test.core, index.family = "jaccard" )
#' 
#' \dontrun{
#' #### using convhulln options
#' # a data set that generates NA (due to errors) with functional.betapart.core
#' data(betatest)
#' # a list of 2 data.frame : comm.test & traits.test
#' comm.test <- betatest$comm.test
#' traits.test <- betatest$traits.test
#' 
#' test <- functional.betapart.core(x = comm.test, traits = traits.test,
#'                                  return.details = FALSE,
#'                                  multi = FALSE)
#' 
#' any(is.na(test$shared))
#' # no NA because the default option was set to QJ
#' # if we use the default option of qhull (Qt) :
#' test <- functional.betapart.core(x = comm.test, traits = traits.test,
#'                                  return.details = FALSE,
#'                                  multi = FALSE,
#'                                  convhull.opt = list(conv1 = "Qt"))
#' any(is.na(test$shared))
#' # some NA arise
#' }
#' # with functional.betapart.core.pairwise
#' test.core <- functional.betapart.core.pairwise(x = comm.test, traits = traits.test,
#'                                                return.details = FALSE,
#'                                                convhull.opt = list(conv1 = "Qt"))
#' any(is.na(test.core$shared))
#' # here no NA were generated because the volumes of the intersections
#' # were computed only with inter_geom
#' # to know which functions were used, set return.details to TRUE
#' test.core <- functional.betapart.core.pairwise(x = comm.test, traits = traits.test,
#'                                                return.details = TRUE,
#'                                                convhull.opt = list(conv1 = "Qt"))
#' test.core$details$Intersection
#' 
#' #### convhull options
#' traits.test <- cbind(c(1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5) , 
#'                      c(1, 2, 4, 1, 2, 3, 5, 1, 4, 3, 5))
#' dimnames(traits.test) <- list(paste("sp", 1:11, sep=""), c("Trait 1", "Trait 2")) 
#' 
#' comm.test <- matrix(0, 4, 11,
#'                     dimnames = list(c("A", "B", "C", "D"), paste("sp", 1:11, sep="")))
#' comm.test["A", c(1, 2, 4, 5)] <- 1
#' comm.test["B", c(1, 3, 8, 9)] <- 1
#' comm.test["C", c(6, 7, 10, 11)] <- 1
#' comm.test["D", c(2, 4, 7, 9)] <- 1 
#' 
#' # simulating a case with species very close to each other
#' # here species 2, 4, 3, 8  close to species 1 
#' # so it is like commA and commB have only 2 species instead of 4 in the 2D space
#' traits.test0<-traits.test
#' traits.test0[c(2,5),]<-traits.test0[1,]+10^-6
#' traits.test0[c(3,8),]<-traits.test0[1,]+10^-6
#' traits.test0
#' 
#' \dontrun{
#' # trying .core function with default qhull option
#' core.test0 <- functional.betapart.core(x = comm.test, traits = traits.test0, multi=FALSE,
#'                                        convhull.opt = list(conv1 = "Qt"))
#' core.test0 # crashing because of coplanarity
#' 
#' # trying new .core.pairwise with default qhull option for convex hull
#' core.pair_test0 <- functional.betapart.core.pairwise(x = comm.test, traits = traits.test0,
#'                                                      convhull.opt = list(conv1 = "Qt"))
#' 
#' # with default qhull options (Qt) it coud be impossible to compute the functional volumes
#' # this is why 'Qj' is now used as the default option for the convexhull volumes
#' 
#' # but other options could be passed to convexhull
#' 
#' # trying new core.pairwise with Qs for convex hull
#' core.pair_test0_Qs <- functional.betapart.core.pairwise(x = comm.test, 
#'                                                         traits = traits.test0,
#'                                                         convhull.opt = list(conv1= "Qs"))
#' # not working
#' }
#' 
#' # trying new .pairwise with QJ (default option) for convex hull
#' core.pair_test0_Qj <- functional.betapart.core.pairwise(x = comm.test, traits = traits.test0,
#'                                                         convhull.opt = list(conv1 = "QJ"),
#'                                                         return.details = TRUE)
#' # OK and QJ applied to each volume computation
#' core.pair_test0_Qj
#' 
#' # equivalent to 
#' core.pair_test0_Qj <- functional.betapart.core.pairwise(x = comm.test, traits = traits.test0,
#'                                                         return.details = TRUE)
#' 
#' # but QJ could be applied only when error are generated with other options :
#' core.pair_test0_Qja <- functional.betapart.core.pairwise(x = comm.test, traits = traits.test0,
#'                                                          convhull.opt = list(conv1 = "Qt",
#'                                                                              conv2 = "QJ"),
#'                                                          return.details = TRUE)
#' # OK and QJ applied only for one volume computation (community 'B')
#' 
#' core.pair_test0_Qja
#' # numerous intersection had to be computed with inter_rcdd
#' 
#' # the results are comparable
#' all.equal(core.pair_test0_Qj[-9], core.pair_test0_Qja[-9]) # -9 to remove details
#' 
#' # the pairwise functional functional betadiversity
#' functional.beta.pair(core.pair_test0_Qj, index.family = "jaccard")
#' 
#' \dontrun{
#' ##### using internal parallelisation to fasten pairiwse dissimilarity
#' # by default (parallel = FALSE) the code is run in serial
#' test.core.pair <- functional.betapart.core.pairwise(x = comm.test, traits = traits.test, 
#'                                                     parallel = FALSE)
#' # when using internal parallelisation and default options 
#' # it uses half of the cores and 1 task per run (this can be customised)
#' # test.core.pairp <- functional.betapart.core(x = comm.test, traits = traits.test, 
#' #                                            multi = FALSE, return.details = FALSE, 
#' #                                            fbc.step = FALSE, parallel = TRUE)
#' # you can set the number of core to use :
#' test.core.pairp <- 
#'   functional.betapart.core.pairwise(x = comm.test, traits = traits.test, 
#'                                     parallel = TRUE, 
#'                                     opt.parallel = beta.para.control(nc = 2))
#' all.equal(test.core.pair, test.core.pairp)
#' 
#' 
#' # as inter_geom is much more faster than inter_rccd it is useful to increase the number
#' # of calculus run per iteration to limit the time consumed by the cluster itself
#' # you can play on size (this would not have sense here as there is only few computation)
#' 
#'   test.core.pairp <- 
#'     functional.betapart.core.pairwise(x = comm.test, traits = traits.test, 
#'                                       parallel = TRUE, 
#'                                       opt.parallel = 
#'                                       beta.para.control(nc = 2,
#'                                                         size = 100))
#'   all.equal(test.core.pair, test.core.pairp)
#'   
#'   
#' library(microbenchmark)
#'   microbenchmark(
#'     serial = functional.betapart.core.pairwise(comm.test, traits.test),
#'     nc2 = functional.betapart.core.pairwise(comm.test, traits.test,
#'                                             parallel = TRUE,
#'                                             opt.parallel = beta.para.control(nc = 2)),
#'     nc4 = functional.betapart.core.pairwise(comm.test, traits.test, multi = FALSE,
#'                                             parallel = TRUE,
#'                                             opt.parallel = beta.para.control(nc = 4))
#'   )
#' 
#' 
#' # If the number of species is very different among communities
#' # load-balancing parallelisation could be more efficient
#' # especially when the number of community is high
#' test.core.pairp <- 
#'   functional.betapart.core.pairwise(comm.test, traits.test, 
#'                                     parallel = TRUE, 
#'                                     opt.parallel = beta.para.control(nc = 2, LB = TRUE))
#' 
#' # use you can use fork cluster (but not on Windows)
#' test.core.pairp <- 
#'   functional.betapart.core.pairwise(comm.test, traits.test,
#'                                       parallel = TRUE,
#'                                       opt.parallel = 
#'                                         beta.para.control(nc = 2, type = "FORK"))
#' 
#'  
#' # a progress bar can be displayed to asses the evolution of the computations
#'   test.core.pairp <- 
#'     functional.betapart.core.pairwise(comm.test, traits.test,
#'                                       parallel = TRUE,
#'                                       opt.parallel = 
#'                                         beta.para.control(nc = 2, LB = TRUE),
#'                                       progress = TRUE)
#'
#' 
#' # using internal parallelisation is not always useful, especially on small data set
#' # load balancing is very helpful when species richness are highly variable
#' 
#' # Null model using 'external' parallel computing 
#' 
#' # Example 1: pairwise functional beta diversity (functional.beta.pair)
#' # Note that this is an example with a small number of samples and null model 
#' # permutations for illustration.
#' # Real null model analyses should have a much greater number of samples and permutations.
#' 
#' ##### 4 communities in a 3D functional space
#' 
#' comm.test <- matrix(0, 4, 11, dimnames = list(c("A", "B", "C", "D"), 
#'                                               paste("sp", 1:11, sep = "")))
#' comm.test["A", c(1, 2, 4, 5)] <- 1
#' comm.test["B", c(1, 3, 8, 9)] <- 1
#' comm.test["C", c(6, 7, 10, 11)] <- 1
#' comm.test["D", c( 2, 4, 7, 9)] <- 1
#' 
#' set.seed(1)
#' traits.test <- matrix(rnorm(11*3, mean = 0, sd = 1), 11, 3) 
#' dimnames(traits.test) <- list(paste("sp", 1:11, sep = "") , 
#'                               c("Trait 1", "Trait 2", "Trait 3"))
#' 
#' # Required packages
#' library(doSNOW)
#' library(picante)
#' library(foreach)
#' library(itertools)
#' 
#' # define number of permutations for the null model (the usual is 1000)
#' # make sure that nperm/nc is a whole number so that all cores have the same number 
#' # of permutations to work on
#' nperm <- 100
#' 
#' test.score <- functional.betapart.core.pairwise(comm.test, traits.test)
#' 
#' obs.pair.func.dis <- functional.beta.pair(x = test.score, index.family = "sorensen")
#' 
#' # transform functional.beta.pair results into a matrix
#' obs.pair.func.dis <- do.call(rbind, obs.pair.func.dis)
#' 
#' # set names for each pair of site
#' pair_names <- combn(rownames(comm.test), 2, FUN = paste, collapse = "_")
#' colnames(obs.pair.func.dis) <- pair_names
#' 
#' # define number of cores
#' # Use parallel::detectCores() to determine number of cores available in your machine
#' nc <- 2 
#' 
#' # 4 cores would be better (nc <- 4)
#' 
#' # create cluster
#' cl <- snow::makeCluster(nc)
#' 
#' # register parallel backend
#' doSNOW:::registerDoSNOW(cl)
#' 
#' # export necessary variables and functions to the cluster of cores
#' snow::clusterExport(cl = cl, c("comm.test", "traits.test"),
#'                     envir = environment())
#' 
#' # creation of an iterator to run 1 comparaisons on each core at time
#' it <- itertools::isplitIndices(nperm, chunkSize = 10)
#' 
#' null.pair.func.dis <- 
#'   foreach(n = it, .combine = c, 
#'           .packages=c("picante","betapart","fastmatch", "rcdd", "geometry")) %dopar% {
#'             
#'             # it enables to adjust the number of permutations (nt) done on each run
#'             nt <- length(n)
#'             null.pair.temp <- vector(mode = "list", length = nt)
#'             
#'             # for each core "n" perform "nt" permutations
#'             for (j in 1:nt){ 
#'               
#'               # randomize community with chosen null model
#'               # for this particular example we used the "independent swap algorithm" 
#'               # but the user can choose other types of permutation
#'               # or create it's own null model 
#'               null.comm.test <- randomizeMatrix(comm.test, null.model = "independentswap", 
#'                                                 iterations=1000)
#'               
#'               # execute functional.betapart.core function 
#'               null.test.score <- 
#'                 functional.betapart.core.pairwise(null.comm.test, traits = traits.test, 
#'                                                   parallel = FALSE)
#'               # using 'external' parallelisation it is necessary to set parralel to FALSE
#'               
#'               # in this artificial example there are a few combinations of species that 
#'               # the convex hull cannot be calculated due to some odd geometric combination 
#'               # so we need to specify to use the 'QJ' options in case of errors
#'               
#'               # compute the pairwise beta-diversity null values and input them in the 
#'               # temporary result matrix
#'               res <- functional.beta.pair(x = null.test.score, index.family = "sorensen")
#'               null.pair.temp[[j]] <- do.call(rbind, res)
#'             }
#'             #retrieve the results from each core
#'             null.pair.temp
#'           }
#' 
#' # stop cluster
#' snow::stopCluster(cl)
#' 
#' #compute the mean, standard deviation and p-values of dissimilarity metrics
#' # for each pair of site
#' 
#' mean.null.pair.func <- matrix(numeric(), ncol = ncol(obs.pair.func.dis), 
#'                               nrow = nrow(obs.pair.func.dis))
#' sd.null.pair.func <- matrix(numeric(), ncol = ncol(obs.pair.func.dis), 
#'                             nrow = nrow(obs.pair.func.dis))
#' p.pair.func.dis <- matrix(numeric(), ncol = ncol(obs.pair.func.dis), 
#'                           nrow = nrow(obs.pair.func.dis))
#' 
#' # for each one of the 3 null dissimilarity metrics (SIN, SNE and SOR) 
#' for (j in 1:nrow(obs.pair.func.dis)){
#'   matnull <- sapply(null.pair.func.dis, function(x) x[j,])
#'   mean.null.pair.func[j,] <- rowMeans(matnull)
#'   sd.null.pair.func[j,] <- sqrt(rowSums((matnull - mean.null.pair.func[j,])^2)/(nperm-1))
#'   p.pair.func.dis[j,] <- rowSums(matnull >= obs.pair.func.dis[j,]) 
#'   p.pair.func.dis[j,] <- (pmin(p.pair.func.dis[j,],nperm-p.pair.func.dis[j,])+1)/(nperm+1)
#'   # the +1 is to take into account that the observed value is one of the possibilities
#' }
#' 
#' # compute standardized effect sizes
#' ses.pair.func.dis <- (obs.pair.func.dis - mean.null.pair.func)/sd.null.pair.func
#' }

functional.betapart.core.pairwise <- function (x, traits, return.details = TRUE,
                                               parallel = FALSE, 
                                               opt.parallel = beta.para.control(), 
                                               convhull.opt = qhull.opt(),
                                               progress = FALSE)
{
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (!is.numeric(x)) 
    stop("The data in 'x' is not numeric.", call. = TRUE)
  xvals <- unique(as.vector(x))
  if (any(!is.element(xvals, c(0, 1)))) 
    stop("The 'x' table contains values other than 0 and 1: data should be presence/absence.", 
         call. = TRUE)
  if (!is.matrix(traits)) {
    traits <- as.matrix(traits)
  }
  if (is.null(row.names(x))) 
    stop("'x' should have row names with site names", call. = TRUE)
  if (is.null(colnames(x))) 
    stop("'x' should have column names with species names", 
         call. = TRUE)
  if (is.null(row.names(traits))) 
    stop("'traits' should have row names with species names", 
         call. = TRUE)
  if (is.null(colnames(traits))) 
    stop("'traits' should have columns names with traits names", 
         call. = TRUE)
  if (any(colnames(x) != row.names(traits))) 
    stop("Species names in 'x' and 'traits' must be identical (including order)", 
         call. = TRUE)
  if (!is.numeric(traits)) 
    stop("The data in 'traits' is not numeric.", call. = TRUE)
  if (any(is.na(traits))) 
    stop("NA are not allowed in 'traits'", call. = TRUE)
  if (ncol(x) != nrow(traits)) 
    stop("Number of species in 'x' and 'traits' must be identical", 
         call. = TRUE)
  # definition of qhull options
  qopts <- qhull.opt()
  if (!missing(convhull.opt)) {
    if (!is.list(convhull.opt))
      stop("convhull.opt must be a list")
    if (length(convhull.opt) > 2)
      stop("convhull.opt must be a list of two elements : conv1 & conv 2")
    if (!all(names(convhull.opt) %in% c("conv1", "conv2")))
      stop("convhull.opt must a list of two elements : conv1 & conv 2")
    if ((!is.null(convhull.opt$conv1) & !is.null(convhull.opt$conv2)) && 
        isTRUE(all.equal(convhull.opt$conv1, convhull.opt$conv2)))
      stop("The two elements of convhull.opt must be different when not null")
    qopts[names(convhull.opt)] <- convhull.opt
  }
  if (!is.null(qopts$conv1)){
    qhull.opt1 <- qopts$conv1
    if (!is.character(qhull.opt1))
      stop("The first element of convhull.opt must be a character vector")
    if (!("FA" %in% qhull.opt1))
      qhull.opt1 <- c("FA", qhull.opt1)
    qhull.opt1 <- paste(qhull.opt1, collapse = " ")
  }else{
    qhull.opt1 <- "FA"
  }
  if (!is.null(qopts$conv2)) {
    qhull.opt2 <- qopts$conv2
    if (!is.character(qhull.opt2))
      stop("The second element of convhull.opt must be a character vector")
    if (!("FA" %in% qhull.opt2))
      qhull.opt2 <- c("FA", qhull.opt2)
    qhull.opt2 <- paste(qhull.opt2, collapse = " ")
  }else{
    qhull.opt2 <- NULL
  }
  if (parallel) {
    control.para <- beta.para.control()
    if (!missing(opt.parallel)) {
      control.para[names(opt.parallel)] <- opt.parallel
    }
    nc <- control.para$nc
    if (!is.numeric(nc)) 
      stop("nc must be numeric (integer)", call. = TRUE)
    nc <- as.integer(nc)
    type <- control.para$type
    if (!type %in% c("SOCK", "PSOCK", "FORK")) 
      stop("type only supoort (P)SOCK or FORK", call. = TRUE)
    if (type == "FORK" && Sys.info()["sysname"] == "Windows") 
      stop("Only SOCK clusters are enabled on Windows", 
           call. = TRUE)
    if (type =="PSOCK") 
      type <- "SOCK"
    LB <- control.para$LB
    if (!is.logical(LB)) 
      stop("LB must be logical", call. = TRUE)
    size <- control.para$size
    if (!is.null(size) && !is.numeric(size)) 
      stop("size must be numeric (integer)", call. = TRUE)
  }
  
  D <- ncol(traits)
  Si <- rowSums(x)
  if (any(Si <= D)) 
    stop(paste("'community ", row.names(x)[which(Si <= D)], 
               " must contain at least ", D + 1, " species", sep = ""))
  N <- nrow(x)
  if (N < 2) 
    stop("Computing dissimilairty requires at least 2 communities", 
         call. = TRUE)
  FRi <- numeric(N)
  names(FRi) <- row.names(x)
  
  if (return.details) {
    details <- vector(mode = "list", length = 2L)
    names(details) <- c("CH", "intersections")
    
    details$CH <- vector(mode = "list", 2L)
    names(details$CH) <- c("FRi", "coord_vertices")
    details$CH$coord_vertices <- vector(mode = "list", N)
    names(details$CH$coord_vertices) <- rownames(x)
    
    details$intersections <- vector("list", 3L)
    names(details$intersections) <- c("combinations", "volumes", "coord_vertices")
    
  }
  
  if (progress) {
    cat("Serial computing of convex hulls shaping assemblages with conv1\n")
    pb <- txtProgressBar(max = N, style = 3)
    progresse <- function(n) setTxtProgressBar(pb, n)
  }
  if (!return.details) {
    for (i in 1:N) { # first try with qhull.opt1
      tr_i <- traits[which(x[i, ] == 1), ]
      FRi[i] <- tryCatch(convhulln(tr_i, options = qhull.opt1)$vol, 
                         error = function(...) NA)
      if (progress)
        progresse(i)
    }
  }else {
    friqopt <- rep(sub("FA ", "", qhull.opt1), N)
    for (i in 1:N) { # first try with qhull.opt1
      tr_i <- traits[which(x[i, ] == 1), ]
      ch <- tryCatch(convhulln(tr_i, options = qhull.opt1), 
                     error = function(...) NA)
      if (!is.na(ch[1])) {
        FRi[i] <- ch$vol
        details$CH$coord_vertices[[i]] <- ch$p[unique(c(ch$hull)),]
      }else {
        FRi[i] <- NA
      }
      if (progress)
        progresse(i)
    }
  }
  if (progress)
    cat("\n")
  if (any(nna <- which(is.na(FRi))) & !is.null(qhull.opt2)) {
    if (progress) {
      cat("Serial computing of convex hulls shaping assemblages with conv2\n")
      pb <- txtProgressBar(max = length(nna), style = 3)
      progresse <- function(n) setTxtProgressBar(pb, n)
    }
    if (!return.details){
      for (i in nna) {
        tr_i <- traits[which(x[i, ] == 1), ]
        FRi[i] <- convhulln(tr_i, options = qhull.opt2)$vol
        if (progress)
          progresse(i)
      }
    }else {
      friqopt[nna] <- sub("FA ", "", qhull.opt2)
      for (i in nna) {
        tr_i <- traits[which(x[i, ] == 1), ]
        ch <- convhulln(tr_i, options = qhull.opt2)
        FRi[i] <- ch$vol
        details$CH$coord_vertices[[i]] <- ch$p[unique(c(ch$hull)),]
        if (progress)
          progresse(i)
      }
    }
    if (progress)
      cat("\n")
  }
  if (any(nafri <- is.na(FRi))) {
    nafri <- which(nafri)
    if (length(nafri) < 5.5){
      ncom <- paste0(paste(row.names(x)[nafri], collapse = ", "), ",")
    }else{
      ncom <- paste0(paste(row.names(x)[nafri[1:5]], collapse = ", "), "...")
    }
    stop(sprintf("For communit%s %s it is not possible to compute %s convex hull volume.
         Please modify the options passed to qhull through the convhull.opt argument.", 
                 ifelse(length(nafri)== 1L, "y", "ies") , ncom,
                 ifelse(length(nafri)==1L, "its", "their")),
         call. = FALSE)
  }
  
  sumFRi <- sum(FRi)
  comb2 <- combn(N, 2)
  vol_inter2_mat <- matrix(0, N, N, dimnames = list(row.names(x), 
                                                    row.names(x)))
  if (return.details) {
    FRid <- data.frame(FRi = FRi, qhull.opt = friqopt)
  }
  
  #### parallel ####
  if (parallel) {
    cl <- snow::makeCluster(nc, type = type)
    doSNOW::registerDoSNOW(cl)
    
    if (progress) {
      n <- choose(nrow(x), 2)
      n <- if (is.null(size)){ 
        round(n/nc)
      }else {
        round(n/size)
      }
    }
    
    if (!return.details) {
      if (progress) {
        cat("Parallel computing of intersections between pairs of assemblages with inter_geom\n")
        pb <- txtProgressBar(max = n, style = 3)
        progressb <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progressb)
      }
      #### no details ####
      if (type %in% c("SOCK", "PSOCK")) 
        snow::clusterExport(cl, c("x", "traits",
                                  "inter_geom", "qhull.opt1",
                                  "inter_rcdd"), envir = environment())

      #### 1. inter_geom 
      iter <- if (is.null(size)) 
        isplitCols(rbind(comb2, seq_len(ncol(comb2))), chunks = nc)
      else isplitCols(rbind(comb2, seq_len(ncol(comb2))), chunkSize = size)
      interp <- foreach(u = iter, .packages = c("rcdd", "geometry"),
                        .inorder = !LB, 
                        .options.snow = if (progress) opts else list()) %dopar% {
                          # u <- as.matrix(u)
                          vol <- numeric(ncol(u))
                          for (k in 1:length(vol)) {
                            i <- u[1, k]
                            j <- u[2, k]
                            seti <- traits[which(x[i, ] == 1), ]
                            setj <- traits[which(x[j, ] == 1), ]
                            interij <- tryCatch(inter_geom(seti, setj, qhull.opt = qhull.opt1),
                                                error = function(...) NA)
                            vol[k] <- interij
                          }
                          res <- cbind(vol, u[3,])
                          res
                        }
      if (progress)
        cat("\n")
      vol_inter2 <- do.call(rbind, interp)
      vol_inter2 <- vol_inter2[order(vol_inter2[,2]), 1]
      
      # now what happens if inter_geom fail
      #### 2. inter_rcdd and qhull.opt1
      nvna <- NULL
      nvna2 <- NULL
      
      if (any(nvna <- which(is.na(vol_inter2)))) {
        if (progress) {
          cat("Parallel computing of intersections between pairs of assemblages with inter_rcdd & conv1\n")
          n <- length(nvna)
          pbr <- txtProgressBar(max = n, style = 3)
          progressbr <- function(n) setTxtProgressBar(pbr, n)
          opts <- list(progress = progressbr)
        }
        combna <- comb2[, nvna, drop = FALSE]
        interna <- isplitCols(combna, chunkSize = 1)
        interp2 <- foreach(u = interna, .packages = c("rcdd", "geometry"),
                           .export = c("inter_rcdd"),
                           .inorder = TRUE,
                           .options.snow = if (progress) opts else list()) %dopar% {
                             u <- as.matrix(u)
                             vol <- numeric(ncol(u))
                             for (k in 1:length(vol)) {
                               i <- u[1, k]
                               j <- u[2, k]
                               seti <- traits[which(x[i, ] == 1), ]
                               setj <- traits[which(x[j, ] == 1), ]
                               interij <- tryCatch(inter_rcdd(seti, setj, qhull.opt = qhull.opt1), 
                                                   error = function(...) NA)
                               vol[k] <- interij
                             }
                             vol
                           }
        if (progress)
          cat("\n")
        vol_interna <- do.call(c, interp2)
        vol_inter2[nvna] <- vol_interna
      }
      #### 3. inter_rcdd and qhull.opt2
      if (any(nvna2 <- which(is.na(vol_inter2))) && !is.null(qhull.opt2)) {
        if (progress) {
          cat("Parallel computing of intersections between pairs of assemblage with inter_rcdd & conv2\n")
          n <- length(nvna2)
          pbr2 <- txtProgressBar(max = n, style = 3)
          progressbr2 <- function(n) setTxtProgressBar(pbr2, n)
          opts <- list(progress = progressbr2)
        }
        combna2 <- comb2[, nvna2, drop = FALSE]
        interna2 <- isplitCols(combna2, chunkSize = 1)
        interp3 <- foreach(u = interna2, .packages = c("rcdd", "geometry"),
                           .export = c("inter_rcdd", "qhull.opt2"),
                           .inorder = TRUE) %dopar% {
                             u <- as.matrix(u)
                             vol <- numeric(ncol(u))
                             for (k in 1:length(vol)) {
                               i <- u[1, k]
                               j <- u[2, k]
                               seti <- traits[which(x[i, ] == 1), ]
                               setj <- traits[which(x[j, ] == 1), ]
                               interij <- tryCatch(inter_rcdd(seti, setj, qhull.opt = qhull.opt2), 
                                                   error = function(...) NA)
                               vol[k] <- interij
                             }
                             vol
                           }
        if (progress)
          cat("\n")
        vol_interna2 <- do.call(c, interp3)
        vol_inter2[nvna2] <- vol_interna2
      }
      
    }else { # if return.details
      #### details ####
      if (progress) {
        cat("Parallel computing of intersections between pairs of assemblages with inter_geom_coord\n")
        pb <- txtProgressBar(max = n, style = 3)
        progressb <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progressb)
      }
      if (type %in% c("SOCK", "PSOCK")) 
        snow::clusterExport(cl, c("x", "traits",
                                  "inter_geom_coord", "qhull.opt1",
                                  "inter_rcdd_coord"), envir = environment())
      #### 1. inter_geom
      iter <- if (is.null(size)) 
        isplitCols(rbind(comb2, seq_len(ncol(comb2))), chunks = nc)
      else isplitCols(rbind(comb2, seq_len(ncol(comb2))), chunkSize = size)
      interp <- foreach(u = iter, .packages = c("rcdd", "geometry"),
                        .inorder = !LB, .combine = c,
                        .options.snow = if (progress) opts else list()) %dopar% {
                          u <- as.matrix(u)
                          n <- ncol(u)
                          vol <- numeric(n)
                          coords <- vector("list", n)
                          ncl <- ncol(traits)
                          for (k in 1:length(vol)) {
                            i <- u[1, k]
                            j <- u[2, k]
                            seti <- traits[which(x[i, ] == 1), ]
                            setj <- traits[which(x[j, ] == 1), ]
                            interij <- tryCatch(inter_geom_coord(seti, setj, qhull.opt = qhull.opt1),
                                                error = function(...) NA)
                            if (!is.na(interij[1])) {
                              vol[k] <- interij$vol
                              coords[[k]] <- interij$coord
                            }else {
                              vol[k] <- NA
                              coords[[k]] <- rep(NA, ncl)
                            }
                          }
                          vol <- cbind(vol, u[3, ])
                          res <- list(coords = coords, vol = vol)
                          res
                        }
      if (progress)
        cat("\n")
      
      coords <- do.call(c, interp[seq(1, length(interp), by = 2)])
      vols <- do.call(rbind, interp[seq(2, length(interp), by = 2)])
      ordo <- order(vols[,2])
      coords <- coords[ordo]
      vol_inter2 <- vols[ordo, 1]
      
      inters <- rep("geom", ncol(comb2))
      qinteropts <- rep(sub("FA ", "", qhull.opt1), ncol(comb2))
      
      #### 2. inter_rcdd qhull.opt1
      nvna <- NULL
      nvna2 <- NULL
      
      if (any(nvna <- which(is.na(vol_inter2)))) {
        if (progress) {
          cat("Parallel computing of intersections between pairs of assemblages with inter_rcdd_coord & conv1\n")
          n <- length(nvna)
          pbr <- txtProgressBar(max = n, style = 3)
          progressbr <- function(n) setTxtProgressBar(pbr, n)
          opts <- list(progress = progressbr)
        }
        combna <- comb2[, nvna, drop = FALSE]
        interna <- isplitCols(rbind(combna, nvna), chunkSize = 1)
        interp2 <- foreach(u = interna, .packages = c("rcdd", "geometry"),
                           .export = c("inter_rcdd_coord"),
                           .inorder = !LB, .combine = c,
                           .options.snow = if (progress) opts else list()) %dopar% {
                             
                             i <- u[1,]
                             j <- u[2,]
                             seti <- traits[which(x[i, ] == 1), ]
                             setj <- traits[which(x[j, ] == 1), ]
                             interij <- tryCatch(inter_rcdd_coord(seti, setj, qhull.opt = qhull.opt1),
                                                 error = function(...) NA)
                             if (!is.na(interij[1])) {
                               vol <- interij$vol
                               coords <- interij$coord
                             }else {
                               vol <- NA
                               coords <- rep(NA, ncol(traits))
                             }
                             vol <- cbind(vol, u[3,])
                             res <- list(coords = coords, vol = vol)
                             res
                           }
        coordsa <- interp2[seq(1, length(interp2), by = 2)]
        volsa <- do.call(rbind, interp2[seq(2, length(interp2), by = 2)])
        ordoa <- volsa[,2]
        coords[ordoa] <- coordsa
        vol_inter2[ordoa] <- volsa[, 1]
        
        if (progress)
          cat("\n")
        inters[nvna] <- "rcdd"
      }
      
      #### 3. inter_rcdd qhull.opt2
      if (any(nvna2 <- which(is.na(vol_inter2))) && !is.null(qhull.opt2)) {
        if (progress) {
          cat("Parallel computing of intersections between pairs of assemblages with inter_rcdd_coord & conv2\n")
          n <- length(nvna2)
          pbr <- txtProgressBar(max = n, style = 3)
          progressbr <- function(n) setTxtProgressBar(pbr, n)
          opts <- list(progress = progressbr)
        }
        combna2 <- comb2[, nvna2, drop = FALSE]
        interna2 <- isplitCols(rbind(combna2, nvna2), chunkSize = 1)
        interp3 <- foreach(u = interna2, .packages = c("rcdd", "geometry"),
                           .export = c("inter_rcdd_coord", "qhull.opt2"),
                           .inorder = !LB, .combine = c,
                           .options.snow = if (progress) opts else list()) %dopar% {
                             
                             i <- u[1,]
                             j <- u[2,]
                             seti <- traits[which(x[i, ] == 1), ]
                             setj <- traits[which(x[j, ] == 1), ]
                             interij <- tryCatch(inter_rcdd_coord(seti, setj, qhull.opt = qhull.opt2),
                                                 error = function(...) NA)
                             if (!is.na(interij[1])) {
                               vol <- interij$vol
                               coords <- interij$coord
                             }else {
                               vol <- NA
                               coords <- rep(NA, ncol(traits))
                             }
                             vol <- cbind(vol, u[3,])
                             res <- list(coords = coords, vol = vol)
                             res
                           }
        coordsa2 <- interp3[seq(1, length(interp3), by = 2)]
        volsa2 <- do.call(rbind, interp3[seq(2, length(interp3), by = 2)])
        ordoa2 <- volsa2[,2]
        coords[ordoa2] <- coordsa2
        vol_inter2[ordoa2] <- volsa2[, 1]
        
        if (progress)
          cat("\n")
        inters[nvna2] <- "rcdd"
        qinteropts[nvna2] <- sub("FA ", "", qhull.opt2)
      }
    }# not return.details
    snow::stopCluster(cl)
  } # end parallel
  else { # in serial
    if (progress) {
      if (return.details){
        cat("Serial computing of intersections between pairs of assemblages with inter_geom_coord\n")
      }else{
        cat("Serial computing of intersections between pairs of assemblages with inter_geom\n")  
        }
      # initialisation of the progressbar for the interaction
      n <- ncol(comb2)
      pbss <- txtProgressBar(max = n, style = 3)
      progressbs <- function(n) setTxtProgressBar(pbss, n)
      # opts <- list(progress = progressb)
    }
    #### serial ####
    vol_inter2 <- numeric(N)
    #### no details ####
    if (!return.details){
      for (k in 1:ncol(comb2)) {
        i <- comb2[1, k]
        j <- comb2[2, k]
        seti <- traits[which(x[i, ] == 1), ]
        setj <- traits[which(x[j, ] == 1), ]
        vol_inter2[k] <-  tryCatch(inter_geom(seti, setj, qhull.opt = qhull.opt1), 
                                   error = function(...) NA)
        if (progress)
          progressbs(k)
      }
      if (progress)
        close(pbss)
    }else {
      nct <- ncol(traits)
      coords <- vector("list", ncol(comb2))
      inters <- rep("geom", ncol(comb2))
      qinteropts <- rep(sub("FA ", "", qhull.opt1), ncol(comb2))
      #### details ####
      for (k in 1:ncol(comb2)) {
        i <- comb2[1, k]
        j <- comb2[2, k]
        seti <- traits[which(x[i, ] == 1), ]
        setj <- traits[which(x[j, ] == 1), ]
        interij <- tryCatch(inter_geom_coord(seti, setj, qhull.opt = qhull.opt1),
                            error = function(...) NA)
        if (!is.na(interij[1])) {
          vol_inter2[k] <- interij$vol
          coords[[k]] <- interij$coord
        }else {
          vol_inter2[k] <- NA
          coords[[k]] <- rep(NA, nct)
        }
        if (progress)
          progressbs(k)
      }
      if (progress)
        close(pbss)
    }
    nvna <- which(is.na(vol_inter2))
    nvn <- setdiff(1:ncol(comb2), nvna)
    if (any(nvna)) {
      nvna <- sort(nvna)
      if (progress) {
        if (return.details) {
          cat("Serial computing of intersections between pairs of assemblages with inter_rcdd_coord & qhull.opt1\n")
        }else{
          cat("Serial computing of intersections between pairs of assemblages with inter_rcdd & qhull.opt1\n")
        }
        pbs <- txtProgressBar(max = length(nvna), style = 3)
        progressbs <- function(n) setTxtProgressBar(pbs, n)
        u <- 1
      }
      if (!return.details) {
        for (k in nvna) {
          i <- comb2[1, k]
          j <- comb2[2, k]
          seti <- traits[which(x[i, ] == 1), ]
          setj <- traits[which(x[j, ] == 1), ]
          vol_inter2[k] <-  inter_rcdd(seti, setj, qhull.opt1)
          if (progress) {
            progressbs(u)
            u <- u+1
          }
        }
        if (progress)
          cat("\n")
      }
      else{
        for (k in nvna) {
          i <- comb2[1, k]
          j <- comb2[2, k]
          seti <- traits[which(x[i, ] == 1), ]
          setj <- traits[which(x[j, ] == 1), ]
          interij <- tryCatch(inter_rcdd_coord(seti, setj, qhull.opt = qhull.opt1),
                              error = function(...) NA)
          if (!is.na(interij[1])) {
            vol_inter2[k] <- interij$vol
            coords[[k]] <- interij$coord
          }else {
            vol_inter2[k] <- NA
            coords[[k]] <-  rep(NA, nct)
          }
          if (progress) {
            progressbs(u)
            u <- u+1
          }
        }
        if (progress)
          cat("\n")
        inters[nvna] <- "rcdd"
      }
      if (any(nvna2 <- which(is.na(vol_inter2))) && !is.null(qhull.opt2)) {
        if (progress) {
          if (return.details){
            cat("Serial computing of intersections between pairs of assemblages with inter_rcdd_coord & qhull.opt2\n")
          }else {
            cat("Serial computing of intersections between pairs of assemblages with inter_rcdd & qhull.opt2\n")
          }
          pbs <- txtProgressBar(max = length(nvna2), style = 3)
          progressbs <- function(n) setTxtProgressBar(pbs, n)
          u <- 1
        }
        if (!return.details){
          for (k in nvna2) {
            i <- comb2[1, k]
            j <- comb2[2, k]
            seti <- traits[which(x[i, ] == 1), ]
            setj <- traits[which(x[j, ] == 1), ]
            vol_inter2[k] <-  inter_rcdd(seti, setj, qhull.opt2)
            if (progress) {
              progressbs(u)
              u+1
            }
          }
          if (progress)
            cat("\n")# end for nvna2
        }else {
          for (k in nvna2) {
            i <- comb2[1, k]
            j <- comb2[2, k]
            seti <- traits[which(x[i, ] == 1), ]
            setj <- traits[which(x[j, ] == 1), ]
            interij <- tryCatch(inter_rcdd_coord(seti, setj, qhull.opt = qhull.opt2),
                                error = function(...) NA)
            if (!is.na(interij[1])) {
              vol_inter2[k] <- interij$vol
              coords[[k]] <- interij$coord
            }else {
              vol_inter2[k] <- NA
              coords[[k]] <-  rep(NA, nct)
            }
            if (progress) {
              progressbs(u)
              u+1
            }
          }
          if (progress)
            cat("\n")
          if(return.details) {
            qinteropts[nvna2] <- sub("FA ", "", qhull.opt2)
          }
        }
      } # end nvna2
    } # end nvna
  } # end serial
  if (return.details) {
    details$CH$FRi <- FRid
    
    couple <- paste0(sprintf(paste0("%0", nchar(N), "d"), comb2[1,]),
                     "_",
                     sprintf(paste0("%0", nchar(N), "d"), comb2[2,]))
    details$intersections$combinations <- vector("list", 1L)
    details$intersections$combinations[[1]] <-  data.frame(Comms = couple, Inter = inters, 
                                                           qhull.opt = qinteropts)
    details$intersections$volumes <- vector("list", 1L)
    details$intersections$volumes[[1]] <- vol_inter2
    
    details$intersections$coord_vertices <- vector("list", 1L)
    details$intersections$coord_vertices[[1]] <- coords

    names(details$intersections$coord_vertices[[1]]) <- 
      details$intersections$combinations[[1]]$Comms
  }else {
    details <- NA
  }
  vol_inter2 <- pmax(vol_inter2, 0)
  comb2 <- t(comb2)
  vol_inter2_mat[comb2[,2:1]] <- vol_inter2
  shared <- vol_inter2_mat
  not.shared <- matrix(0, N, N, dimnames = list(row.names(x), 
                                                row.names(x)))
  not.shared[comb2] <- FRi[comb2[, 1]] - vol_inter2_mat[comb2[, 2:1, drop = FALSE]]
  not.shared[comb2[, 2:1, drop=FALSE]] <- FRi[comb2[, 2]] -
    vol_inter2_mat[comb2[, 2:1, drop = FALSE]]  
  
  if (any(not.shared < 0 & abs(not.shared)>1e-3))
    stop("Negative volume computed")
  not.shared <- pmax(not.shared, 0)
  
  sum.not.shared <- not.shared + t(not.shared)
  max.not.shared <- pmax(not.shared, t(not.shared))
  min.not.shared <- pmin(not.shared, t(not.shared))
  
  functional.computations <- list(sumFRi = sumFRi, FRt = NA, a = NA, 
                                  shared = shared, not.shared = not.shared,
                                  sum.not.shared = sum.not.shared, max.not.shared = max.not.shared,
                                  min.not.shared = min.not.shared,
                                  details = details)
  
  class(functional.computations) <- "functional.betapart"
  return(functional.computations)
}


#' Internal function to compute convexhull volume
#'
#' @description Estimation of the convexhull volume of the intersection of two hypervolumes based on the \code{\link[geometry]{intersectn}} function 
#'
#' @param ps1 A matrix of coordinates.
#' @param ps2 A second matrix of coordinates.
#' @param options Options pass to \code{\link[geometry]{halfspacen}}.
#' @param tol Tolerance, see \code{\link[geometry]{intersectn}}.
#' @param fp Coordinates of feasible point (NULL).
#' @param qhull.opt qhull options.
#' 
#' @return the convex hull volume of the intersection of two hypervolumes
#' 
#' @seealso \code{\link{inter_rcdd}}, \code{\link{intersectn}}
#' 
#' @examples
#' \dontrun{mat1 <- matrix(runif(30), 10)
#' mat2 <- matrix(runif(30), 10)
#' inter_geom(mat1, mat2)}
inter_geom <- function (ps1, ps2, options = "Tv", tol = 0,
                        fp = NULL, qhull.opt = "n FA") 
{
  if (!grepl("^n ", qhull.opt)) {
    qhull.opt <- paste0("n ", qhull.opt)
  }
  distinct <- any(apply(ps1, 2, min) > apply(ps2, 2, max)) || 
    any(apply(ps1, 2, max) < apply(ps2, 2, min))
  if (distinct) {
    vol <- 0
    return(vol)
  }
  ch1 <- convhulln(ps1, "n FA")
  ch2 <- convhulln(ps2, "n FA")
  ch1s <- ch1
  ch2s <- ch2
  if (is.null(fp)) {
    fp <- feasible.point(ch1s, ch2s, tol = tol)
    if (all(is.na(fp))) {
      stop("all fp is.na")
    }
  }
  else {
    if (!is.numeric(fp)) {
      stop("fp should be numeric")
    }
    if (length(fp) != ncol(ps1)) {
      stop("fp should have same dimension as ps1 and ps2")
    }
  }
  ps <- halfspacen(rbind(ch1s$normals, ch2s$normals), fp, options = options)
  if (all(is.na(ps)) || is.null(ps)) {
    stop("ps problem")
  }
  if (tol != 0) {
    ps <- round(ps, ceiling(-log10(abs(tol/2))))
    ps <- ps[!duplicated(ps), ]
  }
  ch <- convhulln(ps, qhull.opt)
  if ((ch$vol > ch1$vol * (1 + 1e-04))) {
    warning("Volume of final intersection hull is bigger than first of the original hulls\n", 
            "ch1 vol = ", ch1$vol, "\n", "ch vol = ", ch$vol, 
            "\n", "Returning ch1")
    ch <- ch1
  }
  if ((ch$vol > ch2$vol * (1 + 1e-04))) {
    warning("Volume of final intersection hull is bigger than first of the original hulls\n", 
            "ch2 vol = ", ch2$vol, "\n", "ch vol = ", ch$vol, 
            "\n", "Returning ch2")
    ch <- ch2
  }
  vol <- ch$vol
  return(vol)
}

#' Internal function to compute convexhull volume
#' @description Estimation of the convexhull volume of the intersection of two hypervolumes based on rcdd functions
#' 
#' @param set1 A matrix of coordinates 
#' @param set2 A matrix of coordinates 
#' @param qhull.opt Qhull options, see \url{http://www.qhull.org/html/qh-optq.htm}
#' @param conv2 A function applyed if the convexhull function crashes
#'
#' @return A volume corresponding to the intersection of the two hypervolumes
#' 
#' @seealso \code{\link{inter_geom}}
#' @examples
#' \dontrun{mat1 <- matrix(runif(30), 10)
#' mat2 <- matrix(runif(30), 10)
#' inter_rcdd(mat1, mat2)}

inter_rcdd <- function(set1, set2, qhull.opt = "FA", conv2 = function(...) NA) {
  set1rep <- d2q(cbind(0, cbind(1, set1)))
  set2rep <- d2q(cbind(0, cbind(1, set2)))
  polytope1 <- redundant(set1rep, representation = "V")$output
  polytope2 <- redundant(set2rep, representation = "V")$output
  H_chset1 <- scdd(polytope1, representation = "V")$output
  H_chset2 <- scdd(polytope2, representation = "V")$output
  H_inter <- rbind(H_chset1, H_chset2)
  V_inter <- scdd(H_inter, representation = "H")$output
  vert_1n2 <- q2d(V_inter[, -c(1, 2)])
  # coord_vert_inter <- rep(NA, ncol(set1))
  vol_inter <- 0 
  if (is.matrix(vert_1n2)) { 
    if (nrow(vert_1n2) > ncol(vert_1n2)) {
      # coord_vert_inter <- vert_1n2
      vol_inter <- tryCatch(convhulln(vert_1n2, qhull.opt)$vol,
                            error = function(...) conv2(vert_1n2))
    }
  }
  return(vol_inter)
}

#' Internal function to compute convexhull volume and vertice coordinates
#' @description Estimation of the convexhull volume and the vertices of the intersection of two hypervolumes based on geometry functions
#' 
#' @param ps1 A matrix of coordinates.
#' @param ps2 A second matrix of coordinates.
#' @param options Options pass to \code{\link[geometry]{halfspacen}}.
#' @param tol Tolerance, see \code{\link[geometry]{intersectn}}.
#' @param fp Coordinates of feasible point (NULL).
#' @param qhull.opt qhull options.
#'
#' @return {A list of 2 elements.
#' \describe{\item{coord}{ the vertice coordinates}
#'   \item{vol}{ a volume corresponding to the intersection of the two hypervolumes}
#'   }
#' } 
#' @seealso \code{\link{inter_rcdd_coord}}
#' @examples
#' \dontrun{mat1 <- matrix(runif(30), 10)
#' mat2 <- matrix(runif(30), 10)
#' inter_geom_coord(mat1, mat2)}

inter_geom_coord <- function (ps1, ps2, options = "Tv", tol = 0,
                              fp = NULL, qhull.opt = "n FA") 
{
  if (!grepl("^n ", qhull.opt)) {
    qhull.opt <- paste0("n ", qhull.opt)
  }
  distinct <- any(apply(ps1, 2, min) > apply(ps2, 2, max)) || 
    any(apply(ps1, 2, max) < apply(ps2, 2, min))
  if (distinct) {
    coord <- rep(NA, ncol(ps1))
    vol <- 0
    res <- list(coord = coord, vol = vol)
    return(res)
  }
  ch1 <- convhulln(ps1, "n FA")
  ch2 <- convhulln(ps2, "n FA")
  ch1s <- ch1
  ch2s <- ch2
  if (is.null(fp)) {
    fp <- feasible.point(ch1s, ch2s, tol = tol)
    if (all(is.na(fp))) {
      stop("all fp is.na")
    }
  }
  else {
    if (!is.numeric(fp)) {
      stop("fp should be numeric")
    }
    if (length(fp) != ncol(ps1)) {
      stop("fp should have same dimension as ps1 and ps2")
    }
  }
  ps <- halfspacen(rbind(ch1s$normals, ch2s$normals), fp, options = options)
  if (all(is.na(ps)) || is.null(ps)) {
    stop("ps problem")
  }
  if (tol != 0) {
    ps <- round(ps, ceiling(-log10(abs(tol/2))))
    ps <- ps[!duplicated(ps), ]
  }
  ch <- convhulln(ps, qhull.opt)
  if ((ch$vol > ch1$vol * (1 + 1e-04))) {
    warning("Volume of final intersection hull is bigger than first of the original hulls\n", 
            "ch1 vol = ", ch1$vol, "\n", "ch vol = ", ch$vol, 
            "\n", "Returning ch1")
    ch <- ch1
  }
  if ((ch$vol > ch2$vol * (1 + 1e-04))) {
    warning("Volume of final intersection hull is bigger than first of the original hulls\n", 
            "ch2 vol = ", ch2$vol, "\n", "ch vol = ", ch$vol, 
            "\n", "Returning ch2")
    ch <- ch2
  }
  coord <- ch$p[unique(c(ch$hull)),]
  vol <- ch$vol
  res <- list(coord = coord, vol = vol)
  return(res)
}

#' Internal function to compute convexhull volume and vertice coordinates
#' @description Estimation of the convexhull volume and the vertices of the intersection of two hypervolumes based on rcdd functions
#' 
#' @param set1 A matrix of coordinates 
#' @param set2 A matrix of coordinates 
#' @param qhull.opt Qhull options, see \url{http://www.qhull.org/html/qh-optq.htm}
#' @param conv2 A function applyed if the convexhull function crashes
#'
#' @return {A list of 2 elements.
#' \describe{\item{coord}{ the vertice coordinates}
#'   \item{vol}{ a volume corresponding to the intersection of the two hypervolumes}
#'   }
#' } 
#' @seealso \code{\link{inter_geom_coord}}
#' @examples
#' \dontrun{mat1 <- matrix(runif(30), 10)
#' mat2 <- matrix(runif(30), 10)
#' inter_rcdd_coord(mat1, mat2)}

inter_rcdd_coord <- function(set1, set2, qhull.opt = "FA", 
                             conv2 = function(...) NA) {
  set1rep <- d2q(cbind(0, cbind(1, set1)))
  set2rep <- d2q(cbind(0, cbind(1, set2)))
  polytope1 <- redundant(set1rep, representation = "V")$output
  polytope2 <- redundant(set2rep, representation = "V")$output
  H_chset1 <- scdd(polytope1, representation = "V")$output
  H_chset2 <- scdd(polytope2, representation = "V")$output
  H_inter <- rbind(H_chset1, H_chset2)
  V_inter <- scdd(H_inter, representation = "H")$output
  vert_1n2 <- q2d(V_inter[, -c(1, 2)])
  coord_vert_inter <- rep(NA, ncol(set1))
  vol_inter <- 0 
  if (is.matrix(vert_1n2)) { 
    if (nrow(vert_1n2) > ncol(vert_1n2)) {
      vol_inter <- tryCatch(convhulln(vert_1n2, qhull.opt)$vol, 
                            error = function(...) conv2(vert_1n2))
      coord_vert_inter <- vert_1n2
    }
  }
  res <- list(coord = coord_vert_inter, vol = vol_inter)
  return(res)
}

#' Specifying control Values for internal parallel cluster
#'
#' @description Set the values used to generate the parallel cluster.
#' 
#' @param nc number of cores to use. Default is half of the available cores.
#' @param conv2 character - the type of cluster to be used, either "SOCK", "PSOCK" or "FORK" (not available on Windows).
#' @param LB logical indicating if load balancing has to be used. Default is \code{TRUE}.
#' @param size number of operation run on each core at each iteration. Default is \codre{1}.
#'
#' @return a list with 4 elements.
#' 
#' @examples 
#' \dontrun{str(beta.para.control(nc = 2, LB = FALSE))}
#' 
beta.para.control <- function(nc = floor(parallel::detectCores()/2), type = "SOCK", 
                              LB = TRUE, size = 1) {
  list(nc = nc, type = type, LB = LB, size = size)
}

#' Specifying control Values for convexhull volume estimation
#'
#' @description Set the default value to use inside the \code{functional.betapart.core} and \code{functional.betapart.pairwise} functions. 
#' It defined the options to use with the \code{convhulln} fonction.
#' 
#' @param conv1 A character vector specifying qhull options. 
#' @param conv2 A character vector specifying qhull options to use if the internal computation of convexhull volumes generates an error. By default (\code{NULL}) \code{NA} is returned.
#'
#' @details \code{conv1} defined options that will be use systematically when calling convhulln to estimate the convexhull volume of the intersection, while
#' \code{conv2} whill be used only if the first call to convhulln would have generated an error. For the complete list of possible options see For the list of options see  \url{http://www.qhull.org/html/qh-optq.htm}. 
#' By default, no option are passed which would generates \code{NA} if internal error. Setting one of the two elements to "QJ" can solve different issues 
#' with very close numerical estimation (difference lower than 1e-4 in our tests).
#' 
#' @return A named list of two elements.
#' 
#' @seealso \code{\link{inter_rcdd}}, \code{\link{inter_geom}}, \code{\link{convhulln}}.
#' 
#' @examples
#' \dontrun{qhull.opt()
#' qhull.opt(conv1 = 'QJ')
#' qhull.opt(conv1 = "Qt", conv2 = 'QJ')}
qhull.opt <- function(conv1 = "QJ", conv2 = NULL) {
  list(conv1 = conv1, conv2 = conv2)
}

