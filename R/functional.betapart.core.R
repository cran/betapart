functional.betapart.core <- function (x, traits, multi = TRUE, warning.time = TRUE, return.details = FALSE, 
                      fbc.step = FALSE, parallel = FALSE, 
                      opt.parallel = beta.para.control())
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
  # also both matrices should have rows and columns names
  if (is.null(row.names(x)))
    stop("'x' should have row names with site names", call. = TRUE)
  if (is.null(colnames(x)))
    stop("'x' should have column names with species names", call. = TRUE)
  if (is.null(row.names(traits)))
    stop("'traits' should have row names with species names", call. = TRUE)
  if (is.null(colnames(traits)))
    stop("'traits' should have columns names with traits names", call. = TRUE)
  # species should be in the same order in columns of 'x' and rows of 'traits' to ensure proper filtering
  if ( any(colnames(x) != row.names(traits) ) )
    stop("Species names in 'x' and 'traits' must be identical (including order)",
         call. = TRUE)
  if (!is.numeric(traits)) 
    stop("The data in 'traits' is not numeric.", call. = TRUE)
  if (any(is.na(traits))) 
    stop("NA are not allowed in 'traits'", call. = TRUE)
  if (ncol(x) != nrow(traits)) 
    stop("Number of species in 'x' and 'traits' must be identical", 
         call. = TRUE)
  if (multi) {
    if (fbc.step) {
      fbc.step <- FALSE
      warnings("As multi = TRUE, fbc.step was set to FALSE")
    }
  }
  
  if (parallel){
    control.para <- beta.para.control()
    if(!missing(opt.parallel)) {
      control.para[names(opt.parallel)] <- opt.parallel
    }
    nc <- control.para$nc
    if (!is.numeric(nc))
      stop("nc must be numeric (integer)", call. = TRUE)
    nc <- as.integer(nc)
    type <- control.para$type
    if (!type %in% c("SOCK", "PSOCK", "FORK"))
      stop("type only supoort (P)SOCK or FORK", call. = TRUE)
    if (type == "FORK" && Sys.info()['sysname'] == "Windows")
      stop("Only SOCK clusters are enabled on Windows", call. = TRUE)
    LB <- control.para$LB
    if (!is.logical(LB))
      stop("LB must be logical", call. = TRUE)
    size <- control.para$size
    if(!is.null(size) && !is.numeric(size))
      stop("size must be numeric (integer)", call. = TRUE)
  }
  
  f1 <- function(z, N, D) {
    comb_z <- combn(N, z, FUN = paste, collapse = "_")
    comb_inter_z2 <- comb_inter[[z - 2]]
    coord_vert_inter_e <- coord_vert_inter[[z - 2]]
    vol_inter_z <- rep(0, length(comb_z))
    coord_vert_inter_z <- list()
    n1 <- sub("_\\d+$", "", comb_z)
    n2 <- sub("^\\d+_", "", comb_z)
    n1 <- fmatch(n1, comb_inter_z2, nomatch = NA)
    n2 <- fmatch(n2, comb_inter_z2, nomatch = NA)
    for (k in 1:length(comb_z)) {
      seti <- coord_vert_inter_e[[n1[k]]]
      setj <- coord_vert_inter_e[[n2[k]]]
      coord_vert_inter_z[[k]] <- rep(NA, D)
      if (!is.na(sum(seti) + sum(setj))) {
        interij <- inter(seti, setj)
        vol_inter_z[k] <- interij$vol_inter
        coord_vert_inter_z[[k]] <- interij$coord_vert_inter
      }
    }
    return(list(comb_z = comb_z, coord_vert_inter_z = coord_vert_inter_z, 
                vol_inter_z = vol_inter_z))
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
  nb.step <- 2
  if (multi) 
    nb.step <- N
  
  if (fbc.step) {
    step.fbc <- as.data.frame(matrix("", nb.step, 1, dimnames = list(c("           FRi", 
                                                                       paste("intersection", 2:nb.step, sep = "_")), c("iteration"))))
    step.fbc[, 1] <- as.character(step.fbc[, 1])
    step.fbc[1, 1] <- paste("0/", N, sep = "")
    for (k in 2:nb.step) step.fbc[k, 1] <- paste("0/", choose(N, 
                                                              k), sep = "")
  }
  FRi <- numeric(N)
  names(FRi) <- row.names(x)
  coord_vert_i <- vector(mode = "list", length = N)
  for (i in 1:N) {
    tr_i <- traits[which(x[i, ] == 1), ]
    ch_i <- geometry::convhulln(tr_i, options="FA") # convex hull (vertices + volume)
    FRi[i] <-ch_i$vol # volume
    coord_vert_i[[i]] <- tr_i[unique(as.integer(ch_i$hull)),]
  }
  names(coord_vert_i) <- row.names(x)
  sumFRi <- sum(FRi)
  inter <- function(set1, set2) {
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
    if (is.matrix(vert_1n2))
      if (nrow(vert_1n2) > ncol(vert_1n2)) {
        coord_vert_inter <- vert_1n2
        vol_inter <- convhulln(vert_1n2, "FA")$vol
      }
    return(list(vol_inter = vol_inter, coord_vert_inter = coord_vert_inter))
  }
  
  comb2 <- combn(N, 2)
  vol_inter2_mat <- matrix(0, N, N, dimnames = list(row.names(x), 
                                                    row.names(x)))
  vol_inter2 <- numeric(ncol(comb2))
  # coord_vert_inter2 <- list()
  coord_vert_inter2 <- vector(mode = "list", length = ncol(comb2))
  if (parallel) {
    combi <- function(x, y) {
      vol <- c(x$vol, y$vol)
      coord <- c(x$coord, y$coord)
      k <- c(x$k, y$k)
      return(list(vol = vol, coord = coord, k = k))
    }
    # paramter to add to the options
    iter <- if (is.null(size)) isplitIndices(ncol(comb2), chunks = nc) else
      isplitIndices(ncol(comb2), chunkSize = size)
    
    # creating the cluster
    cl <- parallel::makeCluster(nc, type = type)
    # register parallel backend
    doParallel::registerDoParallel(cl)
    if (type %in% c("SOCK", "PSOCK"))
      parallel::clusterExport(cl, c("x", "traits", "comb2", "inter"), 
                            envir = environment())
    interp <- foreach(i = iter, .packages = c("rcdd", "geometry"),
                      .combine = combi, .inorder = LB) %dopar% {
                        seqs <- i
                        vol <- numeric(length(seqs))
                        coord <- vector(mode = "list", length = length(seqs))
                        u <- 1
                        for (k in seqs) {
                          i <- comb2[1, k]
                          j <- comb2[2, k]
                          seti <- traits[which(x[i, ] == 1), ]
                          setj <- traits[which(x[j, ] == 1), ]
                          interij <- inter(seti, setj)
                          vol[u] <- interij$vol_inter
                          coord[[u]] <- interij$coord_vert_inter
                          u <- u+1
                        }
                        res <- list(vol = vol, coord = coord, k = seqs)
                        res
                      }
    parallel::stopCluster(cl)
    ordo <- order(interp$k)
    vol_inter2 <- interp$vol
    vol_inter2 <- vol_inter2[ordo]
    coord_vert_inter2 <- interp$coord
    coord_vert_inter2 <- coord_vert_inter2[ordo]
    vol_inter2_mat[t(comb2[2:1,])] <- vol_inter2
    
  }else{ 
    for (k in 1:ncol(comb2)) {
      i <- comb2[1, k]
      j <- comb2[2, k]
      seti <- traits[which(x[i, ] == 1), ]
      setj <- traits[which(x[j, ] == 1), ]
      interij <- inter(seti, setj)
      vol_inter2_mat[j, i] <- interij$vol_inter
      vol_inter2[k] <- interij$vol_inter
      coord_vert_inter2[[k]] <- interij$coord_vert_inter
      if (fbc.step) {
        step.fbc["intersection_2", 1] <- paste(k, "/", ncol(comb2),
                                               sep = "")
        write.table(step.fbc, file = "step.fbc.txt",
                    row.names = T, col.names = F, sep = "\\t")
      }
    }
  }
  shared <- not.shared <- matrix(0, N, N, dimnames = list(row.names(x), row.names(x)))
  for (i in 1:(N - 1)) for (j in (i + 1):N) {
    shared[j, i] <- vol_inter2_mat[j, i]
    not.shared[i, j] <- FRi[i] - vol_inter2_mat[j, i]
    not.shared[j, i] <- FRi[j] - vol_inter2_mat[j, i]
  }
  sum.not.shared <- not.shared + t(not.shared)
  max.not.shared <- pmax(not.shared, t(not.shared))
  min.not.shared <- pmin(not.shared, t(not.shared))
  comb_inter <- list()
  comb_inter[[1]] <- combn(N, 2, paste, collapse = "_")
  coord_vert_inter <- list()
  coord_vert_inter[[1]] <- coord_vert_inter2
  vol_inter <- list()
  vol_inter[[1]] <- vol_inter2
  FRt <- NA
  a <- NA
  if (N > 2 & multi) {
    if (warning.time & N > 10) 
      stop(paste("Computing mulitple functional dissimilarity on more than 10 communities may take a long time. \n    \t\t\t\t\t\t\t\t\tSet 'multi' or 'warning.time' to FALSE"))
    if (warning.time & D > 4) 
      stop(paste("Computing mulitple functional dissimilarity in a", 
                 D, "-dimensions functional space may take a long time. \n    \t\t\t\t\t\t\t\t\tSet 'multi' or 'warning.time' to FALSE"))
    for (z in 3:N) {
      res <- f1(z, N, D)
      comb_inter[[z - 1]] <- res$comb_z
      coord_vert_inter[[z - 1]] <- res$coord_vert_inter_z
      vol_inter[[z - 1]] <- res$vol_inter_z
      if (fbc.step) {
        step.fbc[paste("intersection", z, sep = "_"), 
                 1] <- paste(ncol(res$comb_z), "/", ncol(res$comb_z), 
                             sep = "")
        write.table(step.fbc, file = "step.fbc.txt", 
                    row.names = T, col.names = F, sep = "\\t")
        
      }
    }
    sumvol_sign <- rep(NA, N - 1)
    for (k in 2:N) {
      sumvol_sign[k - 1] <- (-1)^(k - 1) * sum(vol_inter[[k - 
                                                            1]])
    }
    FRt <- sumFRi + sum(sumvol_sign)
    a <- sumFRi - FRt
  }
  details <- NA
  if (return.details) {
    names(coord_vert_i) <- names(FRi)
    CH <- list(FRi = FRi, coord_vertices = coord_vert_i)
    intersections <- list(combinations = comb_inter, volumes = vol_inter, 
                          coord_vertices = coord_vert_inter)
    details <- list(CH = CH, intersections = intersections)
  }
  functional.computations <- list(sumFRi = sumFRi, FRt = FRt, 
                                  a = a, shared = shared, not.shared = not.shared, sum.not.shared = sum.not.shared, 
                                  max.not.shared = max.not.shared, min.not.shared = min.not.shared, 
                                  details = details)
  class(functional.computations) <- "functional.betapart"
  return(functional.computations)
}

# function controlling the parallel framework
# nc : number of cores to use
# type : type of cluster, SOCK or PSOCK (dedied memory) available on all plateforms
# but FORK not available on Windows
# LB : load-balancing (TRUE, by default) or not (FALSE). If species richness or composition 
# is highy varying between sites, then should be prefered.
# size : number of comparisons to do on each core per load. One is the default and recommended when
# species composition/richness is highy variable between sites, otherwise it could be greater.
# size could be NULL, then the comparisons between sites will be performed in one load 
# on each core, with the same number of calculs per cores (approximatively). This is clearly not
# recommended, when the number of sites increases.
beta.para.control <- function(nc = floor(parallel::detectCores()/2), type = "PSOCK",

                              LB = TRUE, size = 1) {

  list(nc = nc, type = type, LB = LB, size = size)
}
