functional.betapart.core<- function (x, traits, multi = TRUE, warning.time = TRUE, 
                                                 return.details = FALSE, fbc.step = FALSE, 
                                                 core.ident = NULL) 
  
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
  if (!is.numeric(traits)) 
    stop("The data in 'traits' is not numeric.", call. = TRUE)
  if (any(is.na(traits))) 
    stop("NA are not allowed in 'traits'", call. = TRUE)
  if (ncol(x) != nrow(traits)) 
    stop("Number of species in 'x' and 'traits' must be identical", 
         call. = TRUE)
  
  f1 <- function(z, N, D, core.ident) {
    
    
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
      if (is.na(sum(seti) + sum(setj)) == F) {
        interij <- inter(seti, setj)
        vol_inter_z[k] <- interij$vol_inter
        coord_vert_inter_z[[k]] <- interij$coord_vert_inter
      }
    }
    return(list(comb_z = comb_z, coord_vert_inter_z = coord_vert_inter_z, vol_inter_z = vol_inter_z))
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
  if (multi == T) 
    nb.step <- N
  step.fbc <- as.data.frame(matrix("", nb.step, 1, dimnames = list(c("           FRi", 
                                                                     paste("intersection", 2:nb.step, sep = "_")), c("iteration"))))
  step.fbc[, 1] <- as.character(step.fbc[, 1])
  step.fbc[1, 1] <- paste("0/", N, sep = "")
  if (fbc.step){
    for (k in 2:nb.step) step.fbc[k, 1] <- paste("0/", choose(N, 
                                                              k), sep = "")
  }
  filec <- !is.null(core.ident)
  
  
  
  FRi <- rep(NA, N)
  names(FRi) <- row.names(x)
  coord_vert_i <- list()
  for (i in 1:N) {
    tr_i <- traits[which(x[i, ] == 1), ]
    if(filec){
      
      vert0 <- convhulln(tr_i, sprintf("Fx TO 'vert_%s.txt'", core.ident))
      vert1 <- scan(sprintf("vert_%s.txt", core.ident), quiet = T)
    }else{
      
      vert0 <- convhulln(tr_i, "Fx TO 'vert.txt'")
      vert1 <- scan("vert.txt", quiet = T)
    }
    
    verti <- (vert1 + 1)[-1]
    
    coord_vert_i[[i]] <- tr_i[verti, ]
    FRi[i] <- convhulln(tr_i[verti, ],"FA")$vol
    if (fbc.step){
      step.fbc["           FRi", 1] <- paste(i, "/", N, sep = "")
      step.fbc[, 1] <- as.character(step.fbc[, 1])
      if(filec){
        write.table(step.fbc, file = sprintf("step.fbc_%s.txt", core.ident), row.names = T, 
                    col.names = F, sep = "\t")
      }else{
        write.table(step.fbc, file = "step.fbc.txt", row.names = T, 
                    col.names = F, sep = "\\t")
      }
    }
  }
  
  
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
    res <- list(coord_vert_inter = coord_vert_inter, vol_inter = vol_inter)
    return(res)
  }
  
  
  comb2 <- combn(N, 2)
  
  vol_inter2_mat <- matrix(0, N, N, dimnames = list(row.names(x), row.names(x)))
  vol_inter2 <- rep(0, ncol(comb2))
  coord_vert_inter2 <- list()
  
  for (k in 1:ncol(comb2)) {
    i <- comb2[1, k]
    j <- comb2[2, k]
    seti <- traits[which(x[i, ] == 1), ]
    setj <- traits[which(x[j, ] == 1), ]
    interij <- inter(seti, setj)
    vol_inter2_mat[j, i] <- interij$vol_inter
    vol_inter2[k] <- interij$vol_inter
    coord_vert_inter2[[k]] <- interij$coord_vert_inter
    
    if (fbc.step){
      step.fbc["intersection_2", 1] <- paste(k, "/", ncol(comb2), 
                                             sep = "")
      if(filec){
        write.table(step.fbc, file = sprintf("step.fbc_%s.txt",core.ident), row.names = T, 
                    col.names = F, sep = "\t")
      }else{
        write.table(step.fbc, file = "step.fbc.txt", row.names = T, 
                    col.names = F, sep = "\\t")
      }
    }
  }
  
  
  matNN <- matrix(0, N, N, dimnames = list(row.names(x), row.names(x)))
  shared <- matNN
  not.shared <- matNN
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
  
  if (N > 2 & multi == T) {
    if (warning.time == T & N > 10)
      stop(paste("Computing mulitple functional dissimilarity on more than 10 communities may take a long time. \n    \t\t\t\t\t\t\t\t\tSet 'multi' or 'warning.time' to FALSE"))
    if (warning.time == T & D > 4)
      stop(paste("Computing mulitple functional dissimilarity in a",
                 D, "-dimensions functional space may take a long time. \n    \t\t\t\t\t\t\t\t\tSet 'multi' or 'warning.time' to FALSE"))
    
    for (z in 3:N) {
      res <- f1(z, N, D, core.ident)
      comb_inter[[z - 1]] <- res$comb_z
      coord_vert_inter[[z - 1]] <- res$coord_vert_inter_z
      vol_inter[[z - 1]] <- res$vol_inter_z
      if (fbc.step) {
        step.fbc[paste("intersection", z, sep = "_"), 
                 1] <- paste(ncol(res$comb_z), "/", ncol(res$comb_z), sep = "")
        if(filec){
          write.table(step.fbc, file = sprintf("step.fbc_%s.txt",core.ident), row.names = T, 
                      col.names = F, sep = "\t")
        }else{
          write.table(step.fbc, file = "step.fbc.txt", row.names = T, 
                      col.names = F, sep = "\\t")
        }
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
  if (return.details == T) {
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


