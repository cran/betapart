beta.multi.abund<-function (x, index.family = "bray") 
{
    index.family <- match.arg(index.family, c("bray", "ruzicka"))
    if (!inherits(x, "betapart.abund")) {
        x <- betapart.core.abund(x)
    }
    maxbibj <- sum(x$pair.max.not.shared.abund)
    minbibj <- sum(x$pair.min.not.shared.abund)
    switch(index.family, bray = {
        beta.bray.bal <- minbibj/(minbibj + x$multiple.shared.abund)
        beta.bray.gra <- (x$multiple.shared.abund/(minbibj + x$multiple.shared.abund)) * ((maxbibj - minbibj)/((2 * 
            x$multiple.shared.abund) + maxbibj + minbibj))
        beta.bray <- (minbibj + maxbibj)/(minbibj + maxbibj + 
            (2 * x$multiple.shared.abund))
        multi <- list(beta.BRAY.BAL = beta.bray.bal, beta.BRAY.GRA = beta.bray.gra, 
            beta.BRAY = beta.bray)
    }, ruzicka = {
        beta.ruz.bal <- (2 * minbibj)/((2 * minbibj) + x$multiple.shared.abund)
        beta.ruz.gra <- (x$multiple.shared.abund/((2 * minbibj) + x$multiple.shared.abund)) * ((maxbibj - 
            minbibj)/((x$multiple.shared.abund) + maxbibj + minbibj))
        beta.ruz <- (minbibj + maxbibj)/(minbibj + maxbibj + 
            x$multiple.shared.abund)
        multi <- list(beta.RUZ.BAL = beta.ruz.bal, beta.RUZ.GRA = beta.ruz.gra, 
            beta.RUZ = beta.ruz)
    })
    return(multi)
}
