functional.betapart.core <- function (x, traits, multi = TRUE, warning.time = TRUE, 
                                      return.details = FALSE, fbc.step = FALSE, 
                                      parallel = FALSE, opt.parallel = beta.para.control(), 
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
  
  def_conv2 <- function(qhull2) {
    if (is.null(qhull2)) {
      return(function(...) NA)
    }else {
      # return(function(x) convhulln(x, options = qhull2)$vol)
      # remplacer par
      function(x) tryCatch(convhulln(x, options = qhull2)$vol, error = function(...) NA)
    }
  }
  conv2 <- def_conv2(qhull.opt2)

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
  
  interi <- function(set1, set2, qhull.opt = "FA", conv2) {
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
        vol_inter <- tryCatch(convhulln(vert_1n2, qhull.opt)$vol, 
                              error = function(...) conv2(vert_1n2))
      }
    return(list(vol_inter = vol_inter, coord_vert_inter = coord_vert_inter))
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
        interij <- interi(seti, setj, qhull.opt1, conv2)
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
  if (multi) {
    if (fbc.step) {
      fbc.step <- FALSE
      warnings("As multi = TRUE, fbc.step was set to FALSE")
    }
    if (warning.time & N > 10) 
      stop(paste("Computing mulitple functional dissimilarity on more than 10 communities may take a long time. \n    \t\t\t\t\t\t\t\t\tSet 'multi' or 'warning.time' to FALSE"))
    if (warning.time & D > 4) 
      stop(paste("Computing mulitple functional dissimilarity in a", 
                 D, "-dimensions functional space may take a long time. \n    \t\t\t\t\t\t\t\t\tSet 'multi' or 'warning.time' to FALSE"))
  }
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
    ch_i <- tryCatch(convhulln(tr_i, options = qhull.opt1), 
                     error = function(...) NA)
    # ch_i <-  convhulln(tr_i, options = qhull.opt1)
    FRi[i] <- if (!is.na(ch_i)[1]) ch_i$vol else NA
    coord_vert_i[[i]] <- if (!is.na(ch_i)[1]) tr_i[unique(as.integer(ch_i$hull)), ] else NA 
  }
  if (any(nna <- which(is.na(FRi))) & !is.null(qhull.opt2)) {
    for (i in nna) {
      tr_i <- traits[which(x[i, ] == 1), ]
      ch_i <- tryCatch(convhulln(tr_i, options = qhull.opt2), 
                       error = function(...) NA)
      FRi[i] <- if (!is.na(ch_i[1])) ch_i$vol else NA
      coord_vert_i[[i]] <- if (!is.na(ch_i)[1]) tr_i[unique(as.integer(ch_i$hull)), ] else NA 
    }
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
                 ifelse(length(nafri)==1L, "its", "their")))
  }
  names(coord_vert_i) <- row.names(x)
  sumFRi <- sum(FRi)
  comb2 <- combn(N, 2)
  
    
  if (parallel) {
    if (progress) {
      cat("Parallel process\n")
      # initialisation of the progressbar for the interaction
        n <- ncol(comb2)
        n <- if (is.null(size)){ 
          round(n/nc)
        }else {
          round(n/size)
        }
      pb <- txtProgressBar(max = n, style = 3)
      progressb <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progressb)
    }
    
    cl <- snow::makeCluster(nc, type = type)
    doSNOW::registerDoSNOW(cl)
    if (type %in% c("SOCK", "PSOCK")) 
      snow::clusterExport(cl, c("x", "traits",
                                    "interi", "qhull.opt1",
                                    "conv2"), envir = environment())
    if (LB){ #  LB
      iter <- if (is.null(size))
        isplitCols(rbind(comb2, seq_len(ncol(comb2))), chunks = nc)
      else isplitCols(rbind(comb2, seq_len(ncol(comb2))), chunkSize = size)
      interp <- foreach(u = iter, .packages = c("rcdd", "geometry", "betapart"), 
                        .inorder = FALSE, 
                        .options.snow = if (progress) opts else list()) %dopar% {
                          u <- as.matrix(u)
                          vol <- numeric(ncol(u))
                          for (k in 1:length(vol)) {
                            i <- u[1, k]
                            j <- u[2, k]
                            seti <- traits[which(x[i, ] == 1), ]
                            setj <- traits[which(x[j, ] == 1), ]
                            vol[k] <- inter_rcdd(seti, setj, qhull.opt1, conv2)
                          }
                          res <- cbind(vol, u[3,])
                          res
                        }
      if (progress)
        cat("\n")
      snow::stopCluster(cl)
      vol_inter2 <- do.call(rbind, interp)
      vol_inter2 <- vol_inter2[order(vol_inter2[,2]), 1]
    }else{ # no LB
      iter <- if (is.null(size)) 
        isplitCols(comb2, chunks = nc)
      else isplitCols(comb2, chunkSize = size)
      interp <- foreach(u = iter, .packages = c("rcdd", "geometry", "betapart"), 
                        .inorder = TRUE, 
                        .options.snow = if (progress) opts else list()) %dopar% {
                          u <- as.matrix(u)
                          vol <- numeric(ncol(u))
                          for (k in 1:length(vol)) {
                            i <- u[1, k]
                            j <- u[2, k]
                            seti <- traits[which(x[i, ] == 1), ]
                            setj <- traits[which(x[j, ] == 1), ]
                            vol[k] <- inter_rcdd(seti, setj, qhull.opt1, conv2)
                          }
                          vol
                        }
      if (progress)
        cat("\n")
      snow::stopCluster(cl)
      vol_inter2 <- do.call(c, interp)
    } #  end LB
    vol_inter2 <- pmax(vol_inter2, 0)
    comb2 <- t(comb2)
    vol_inter2_mat <- matrix(0, N, N, dimnames = list(row.names(x), 
                                                      row.names(x)))
    vol_inter2_mat[comb2[,2:1]] <- vol_inter2
    shared <- vol_inter2_mat
    not.shared <- matrix(0, N, N, dimnames = list(row.names(x), 
                                                  row.names(x)))
    not.shared[comb2] <- FRi[comb2[,1]]-vol_inter2_mat[comb2[,2:1]]
    not.shared[comb2[,2:1]] <- FRi[comb2[,2]]-vol_inter2_mat[comb2[,2:1]]  
    
    not.shared <- pmax(not.shared, 0)
    coord_vert_inter2 <- NA
  }else {
    if (progress) {
      cat("Serial process\n")
      # initialisation of the progressbar for the interaction
      n <- ncol(comb2)
      pb <- txtProgressBar(max = n, style = 3)
      progressb <- function(n) setTxtProgressBar(pb, n)
      # opts <- list(progress = progressb)
    }
    vol_inter2_mat <- matrix(0, N, N, dimnames = list(row.names(x), 
                                                      row.names(x)))
    not.shared <- matrix(0, N, N, dimnames = list(row.names(x), 
                                                  row.names(x)))
    
    vol_inter2 <- numeric(ncol(comb2))
    coord_vert_inter2 <- vector(mode = "list", length = ncol(comb2))
    for (k in 1:ncol(comb2)) {
      
      i <- comb2[1, k]
      j <- comb2[2, k]
      seti <- traits[which(x[i, ] == 1), ]
      setj <- traits[which(x[j, ] == 1), ]
      interij <- interi(seti, setj, qhull.opt1, conv2)
      vol_inter2_mat[j, i] <- interij$vol_inter
      vol_inter2[k] <- interij$vol_inter
      coord_vert_inter2[[k]] <- interij$coord_vert_inter
      if (fbc.step) {
        step.fbc["intersection_2", 1] <- paste(k, "/", 
                                               ncol(comb2), sep = "")
        write.table(step.fbc, file = "step.fbc.txt", 
                    row.names = T, col.names = F, sep = "\\t")
      }
      if (progress)
        progressb(k)
    }
    comb2 <- t(comb2)
    vol_inter2_mat <- pmax(vol_inter2_mat, 0)
    vol_inter2 <- pmax(vol_inter2)
    shared <- vol_inter2_mat
    not.shared[comb2] <- FRi[comb2[,1]]-vol_inter2_mat[comb2[,2:1]]
    not.shared[comb2[,2:1]] <- FRi[comb2[,2]]-vol_inter2_mat[comb2[,2:1]]
    not.shared <- pmax(not.shared, 0)
  }
  sum.not.shared <- not.shared + t(not.shared)
  max.not.shared <- pmax(not.shared, t(not.shared))
  min.not.shared <- pmin(not.shared, t(not.shared))
  
  FRt <- NA
  a <- NA
  
  comb_inter <- list()
  comb_inter[[1]] <- combn(N, 2, paste, collapse = "_")
  coord_vert_inter <- list()
  coord_vert_inter[[1]] <- coord_vert_inter2
  vol_inter <- list()
  vol_inter[[1]] <- vol_inter2
  
  if (N > 2 & multi) {
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

