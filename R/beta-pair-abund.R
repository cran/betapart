
beta.pair.abund<-function (x, index.family = "bray") 
{
    index.family <- match.arg(index.family, c("bray", "ruzicka"))
    if (!inherits(x, "betapart.abund")) {
        x <- betapart.core.abund(x)
    }
    switch(index.family, bray = {
        beta.bray.bal <- x$pair.min.not.shared.abund/(x$pair.min.not.shared.abund + x$pair.shared.abund)
        beta.bray.gra <- ((x$pair.max.not.shared.abund - x$pair.min.not.shared.abund)/((2 * 
            x$pair.shared.abund) + x$pair.sum.not.shared.abund)) * (x$pair.shared.abund/(x$pair.min.not.shared.abund + 
            x$pair.shared.abund))
        beta.bray <- x$pair.sum.not.shared.abund/(2 * x$pair.shared.abund + x$pair.sum.not.shared.abund)
        pairwise <- list(beta.bray.bal = as.dist(beta.bray.bal), beta.bray.gra = as.dist(beta.bray.gra), 
            beta.bray = as.dist(beta.bray))
    }, ruzicka = {
        beta.ruz.bal <- (2 * x$pair.min.not.shared.abund)/((2 * x$pair.min.not.shared.abund) + 
            x$pair.shared.abund)
        beta.ruz.gra <- ((x$pair.max.not.shared.abund - x$pair.min.not.shared.abund)/(x$pair.shared.abund + 
            x$pair.sum.not.shared.abund)) * (x$pair.shared.abund/((2 * x$pair.min.not.shared.abund) + 
            x$pair.shared.abund))
        beta.ruz <- x$pair.sum.not.shared.abund/(x$pair.shared.abund + x$pair.sum.not.shared.abund)
        pairwise <- list(beta.ruz.bal = as.dist(beta.ruz.bal), beta.ruz.gra = as.dist(beta.ruz.gra), 
            beta.ruz = as.dist(beta.ruz))
    })
    return(pairwise)
}
