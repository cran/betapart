beta.sample.abund<-function (x, index.family = "bray", sites = nrow(x), samples = 1) 
{
    index.family <- match.arg(index.family, c("bray", "ruzicka"))
    
    if (sites > nrow(x)) 
        stop("More sites requested for sample than are in the dataset")
    pb <- txtProgressBar(min = 0, max = samples, style = 3)
    results.n <- as.data.frame(matrix(nrow = samples, ncol = 3))
    
    for (i in 1:samples) {
        position <- as.vector(1:nrow(x))
        sample.position <- sample(position, sites)
        x.sample <- x[sample.position,]
        x.beta <- beta.multi.abund(x.sample, index.family)
        results.n[i, ] <- unlist(x.beta)
        setTxtProgressBar(pb, i)
    }
    close(pb)
    names(results.n) <- names(x.beta)
    result <- list(sampled.values = results.n, mean.values = sapply(results.n, 
        mean), sd.values = sapply(results.n, sd))
    return(result)
}