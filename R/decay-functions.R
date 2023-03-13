########################################################################################
##### A function for fitting distance-decay models


decay.model<-function (y, x, model.type = "exponential", y.type = "similarities", 
    perm = 100, st.val = c(1, 0)) 
{
    model.type <- match.arg(model.type, c("exponential", "power", "gompertz"))
    y.type <- match.arg(y.type, c("similarities", "dissimilarities"))
    switch(y.type, similarities = {
        ydis <- y
    }, dissimilarities = {
        ydis <- (1 - y)
    })
    
    nsites <- dim(as.matrix(ydis))[1]
    
    switch(model.type, exponential = {        
    	model.formula <- ydis~a*exp(-b*x)
     st.val <- list(a=st.val[1],b=st.val[2])  

    	}, power = {
    	model.formula <- ydis~a*(x^-b)
     st.val <- list(a=st.val[1],b=st.val[2])  

		}, gompertz ={
		model.formula <- ydis~1-(exp(-b*exp(-c*x)))
     st.val <- list(b=st.val[1],c=st.val[2])  

})
    
    null.x<-rep(1, length(y))
	null.nlm<-nlsLM(ydis~a*(null.x),start=list(a=1))
	null.deviance<-deviance(null.nlm)
	
	model<-nlsLM(model.formula, start=st.val)
    pseudo.r.squared <- 1 - deviance(model)/null.deviance

    perm.dev.R2.distr <- numeric(perm)
    ydismat <- as.matrix(ydis)
   
    for (i in 1:perm) {
        perm.sites <- sample(1:nsites)  
        ydis <- ydismat[perm.sites, perm.sites][lower.tri(ydismat)]
        mod.perm <- nlsLM(model.formula, start=st.val)
        perm.dev.R2.distr[i] <- 1 - deviance(mod.perm)/null.deviance
    }
    
    p.value <- mean(perm.dev.R2.distr > pseudo.r.squared)
    
    switch(y.type, similarities = {
     	parameters <- list(data.y = y, data.x = x, model = model, 
     		model.type = model.type, y.type = y.type, 
     		first.parameter=coefficients(model)[1], 
     		second.parameter=-coefficients(model)[2], 
     		aic=AIC(model), pseudo.r.squared = pseudo.r.squared,
     		p.value = ifelse(p.value == 0, 1/perm, p.value))
    }, dissimilarities = {
     	parameters <- list(data.y = y, data.x = x, model = model, 
     		model.type = model.type, y.type = y.type, 
     		first.parameter=1-coefficients(model)[1], 
     		second.parameter=coefficients(model)[2], 
     		aic=AIC(model), pseudo.r.squared = pseudo.r.squared,
     		p.value = ifelse(p.value == 0, 1/perm, p.value))
    })
        
    class(parameters) <- "decay"
    return(parameters)
}






################################################################################




################################################################################
##### A function for plotting distance-decay models

plot.decay<-function (x, xlim = c(0, max(x$data.x)), ylim = c(0, 1), add = FALSE, 
    remove.dots = FALSE, col = "black", pch = 1, lty = 1, lwd = 5, 
    cex = 1, ...) 
{
    if (!inherits(x, "decay")) {
        stop("The input is not a distance-decay model fitted with decay.model().", 
            call. = TRUE)
    }
    data <- data.frame(as.vector(x$data.x), as.vector(x$data.y))
    if (!remove.dots) {
        pch = pch
    }
    else {
        pch = ""
    }
      
    if (!add) {
	plot(x$data.x, x$data.y, xlim = xlim, ylim = ylim, xlab = "Distance", 
     ylab = ifelse(x$y.type=="similarities", "Similarity", "Disimilarity"),
     col = col, pch = pch, cex = cex, ...)
	}
	
	if (add) {
     points(x$data.x, x$data.y, xlim = xlim, ylim = ylim, xlab = "Distance", 
     ylab = ifelse(x$y.type=="similarities", "Similarity", "Disimilarity"),
     col = col, pch = pch, cex = cex, ...)
	}   
	
	switch(x$y.type, similarities = {
	lines(fitted(x$model)[order(x$data.x)] ~ sort(x$data.x), col = col, lty = lty, lwd = lwd, ...)
	}, dissimilarities = {
	lines(1-fitted(x$model)[order(x$data.x)] ~ sort(x$data.x), col = col, lty = lty, lwd = lwd, ...)	
	})
   	
}



################################################################################




################################################################################
##### A function for bootstrapping the parameters of distance-decay models


boot.coefs.decay<-function (m1, resamples, st.val=c(1,0)) 
{
    if (!inherits(m1, "decay")) {
        stop("The input is not a distance-decay model fitted with decay.model().", 
            call. = TRUE)
    }
    ptm <- proc.time()
    original.coefs <- c(m1$first.parameter, m1$second.parameter)
    names(original.coefs) <- c("first.parameter", "second.parameter")
    boot.coefs <- matrix(nrow = resamples, ncol = 2)
    colnames(boot.coefs) <- c("first.parameter", "second.parameter")
                            
m1.nsites<-dim(as.matrix(m1$data.y))[1]

            
switch(m1$y.type, dissimilarities = {
        m1$data.y <- (1 - m1$data.y)
     })  

## convert dist objects used to build the models into block matrices

m1.y.mat<-as.matrix(m1$data.y)
m1.y.blockmat<-matrix(m1.y.mat[!row(m1.y.mat)==col(m1.y.mat)], nrow=nrow(m1.y.mat)-1, ncol=ncol(m1.y.mat))

m1.x.mat<-as.matrix(m1$data.x)
m1.x.blockmat<-matrix(m1.x.mat[!row(m1.x.mat)==col(m1.x.mat)], nrow=nrow(m1.x.mat)-1, ncol=ncol(m1.x.mat))

# Block bootstrap to compute the variance of parameters 

# Vectors for storage the resampled parameters in the site-block resampling
	m1.par1<-vector(length=resamples)
   	m1.par2<-m1.par1

for (i in 1:resamples){
	
# Index vectors for resampling with replacement
    m1.index<-sample(x = 1:m1.nsites, size=m1.nsites, replace=T)             

# Resampled matrices
    m1.x.resample <- m1.x.blockmat[, m1.index]                    
    m1.y.resample <- m1.y.blockmat[, m1.index]  
 
# Index vector for downsizing resampled matrices
    m1.down.index <- sample(1:length(m1.y.resample), size = m1.nsites*(m1.nsites - 1)/2)  
  

# Downsize resampled matrices

    m1.boot.xdis <- m1.x.resample[m1.down.index]                                         
    m1.boot.ydis <- m1.y.resample[m1.down.index]   
                                                                               

# compute the block bootstrap samples of model parameters
   	
   	switch(m1$model.type, exponential = {        
    	model.formula <- ydis~a*exp(-b*x)
     st.val <- list(a=st.val[[1]],b=st.val[[2]])  

    	}, power = {
    	model.formula <- ydis~a*(x^-b)
     st.val <- list(a=st.val[[1]],b=st.val[[2]])  

		}, gompertz ={
		model.formula <- ydis~1-(exp(-b*exp(-c*x)))
     st.val <- list(b=st.val[[1]],c=st.val[[2]])  

})
   
   
  	x<-m1.boot.xdis
   	ydis<-m1.boot.ydis
	m1.boot<-nlsLM(model.formula,start=st.val)
	
	boot.coefs[i, 1]<-ifelse(m1$y.type=="similarities", coefficients(m1.boot)[1], 1-coefficients(m1.boot)[1])
	boot.coefs[i, 2]<-ifelse(m1$y.type=="similarities", -coefficients(m1.boot)[2], coefficients(m1.boot)[2])
	
		
	}

    mean.boot <- c(mean(boot.coefs[, 1]), mean(boot.coefs[, 2]))
    names(mean.boot) <- c("first.parameter", "second.parameter")
    sd.boot = c(sd(boot.coefs[, 1]), sd(boot.coefs[, 2]))
    names(sd.boot) <- c("first.parameter", "second.parameter")
    result <- list(model.type = m1$model.type, y.type = m1$y.type, 
        boot.coefs = boot.coefs, original.coefs = original.coefs, 
        mean.boot = mean.boot, sd.boot = sd.boot, time = proc.time() - 
            ptm)
    result
}


################################################################################




################################################################################
##### A function to assess differences in parameters between two distance-decay models


zdep<-function(m1, m2, resamples, st.val=c(1,0)){

if (m1$y.type!=m2$y.type) 
        stop("Both models should have the same type of response variable, either similarities or dissimilarities.", 
            call. = TRUE)
if (m1$model.type!=m2$model.type) 
        stop("Both models should fit the same function, either negative exponential, power law, or gompertz.", 
            call. = TRUE)
            
switch(m1$y.type, dissimilarities = {
        m1$data.y <- (1 - m1$data.y)
        m2$data.y <- (1 - m2$data.y)
    })         
                    
                        
m1.nsites<-dim(as.matrix(m1$data.y))[1]
m2.nsites<-dim(as.matrix(m2$data.y))[1]

## convert dist objects used to build the models into block matrices

m1.y.mat<-as.matrix(m1$data.y)
m1.y.blockmat<-matrix(m1.y.mat[!row(m1.y.mat)==col(m1.y.mat)], nrow=nrow(m1.y.mat)-1, ncol=ncol(m1.y.mat))

m1.x.mat<-as.matrix(m1$data.x)
m1.x.blockmat<-matrix(m1.x.mat[!row(m1.x.mat)==col(m1.x.mat)], nrow=nrow(m1.x.mat)-1, ncol=ncol(m1.x.mat))


m2.y.mat<-as.matrix(m2$data.y)
m2.y.blockmat<-matrix(m2.y.mat[!row(m2.y.mat)==col(m2.y.mat)], nrow=nrow(m2.y.mat)-1, ncol=ncol(m2.y.mat))

m2.x.mat<-as.matrix(m2$data.x)
m2.x.blockmat<-matrix(m2.x.mat[!row(m2.x.mat)==col(m2.x.mat)], nrow=nrow(m2.x.mat)-1, ncol=ncol(m2.x.mat))
   

# Block bootstrap to compute the variance of parameters 

# Vectors for storage the resampled parameters in the site-block resampling
	m1.par1<-vector(length=resamples)
   	m1.par2<-m1.par1
   	m2.par1<-m1.par1
   	m2.par2<-m1.par1
 
for (i in 1:resamples){
	
# Index vectors for resampling with replacement
    m1.index<-sample(x = 1:m1.nsites, size=m1.nsites, replace=T)             
    m2.index<-sample(x = 1:m2.nsites, size=m2.nsites, replace=T)

# Resampled matrices
    m1.x.resample <- m1.x.blockmat[, m1.index]                    
    m1.y.resample <- m1.y.blockmat[, m1.index]  
          
    m2.x.resample <- m2.x.blockmat[, m2.index]                    
    m2.y.resample <- m2.y.blockmat[, m2.index] 
        
# Index vectors for downsizing resampled matrices
    m1.down.index <- sample(1:length(m1.y.resample), size = m1.nsites*(m1.nsites - 1)/2)  
    m2.down.index <- sample(1:length(m2.y.resample), size = m2.nsites*(m2.nsites - 1)/2)  

# Downsize resampled matrices

    m1.boot.xdis <- m1.x.resample[m1.down.index]                                         
    m1.boot.ydis <- m1.y.resample[m1.down.index]   
                                          
    m2.boot.xdis <- m2.x.resample[m2.down.index]                                         
    m2.boot.ydis <- m2.y.resample[m2.down.index]                                         

# compute the block bootstrap samples of model parameters
   	
   	switch(m1$model.type, exponential = {        
    	model.formula <- ydis~a*exp(-b*x)
     st.val <- list(a=st.val[[1]],b=st.val[[2]])  

    	}, power = {
    	model.formula <- ydis~a*(x^-b)
     st.val <- list(a=st.val[[1]],b=st.val[[2]]) 

		}, gompertz ={
		model.formula <- ydis~1-(exp(-b*exp(-c*x)))
     st.val <- list(b=st.val[[1]],c=st.val[[2]])  

})
   
   
  	x<-m1.boot.xdis
   	ydis<-m1.boot.ydis
	m1.boot<-nlsLM(model.formula,start=st.val)
	
	x<-m2.boot.xdis
   	ydis<-m2.boot.ydis
	m2.boot<-nlsLM(model.formula,start=st.val)
	
	m1.par1[i]<-coefficients(m1.boot)[1]
	m1.par2[i]<-coefficients(m1.boot)[2]
	m2.par1[i]<-coefficients(m2.boot)[1]
	m2.par2[i]<-coefficients(m2.boot)[2]
		
	}
	
# Z.dep statistic

z.dep.par1<-(m1$first.parameter-m2$first.parameter)/sqrt((var(m1.par1)+var(m2.par1))-2*cov(m1.par1,m2.par1))
z.dep.par2<-(m1$second.parameter-m2$second.parameter)/sqrt((var(m1.par2)+var(m2.par2))-2*cov(m1.par2,m2.par2))


# p values, computed as 2*(P|N|>|Z.dep|)

p.par1<-ifelse(z.dep.par1<=0, 2*pnorm(z.dep.par1,lower.tail = T), 2*pnorm(z.dep.par1,lower.tail = F))
  
p.par2<-ifelse(z.dep.par2<=0, 2*pnorm(z.dep.par2,lower.tail = T), 2*pnorm(z.dep.par2,lower.tail = F))

results<-matrix(ncol=2, nrow=2)

results<-matrix(c(z.dep.par1, z.dep.par2,p.par1, p.par2), ncol=2, nrow=2)
rownames(results)<-c("First parameter", "Second parameter")
colnames(results)<-c("Z.dep", "p value")

return(results)

}	


################################################################################