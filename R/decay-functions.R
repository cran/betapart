########################################################################################
##### A function for fitting distance-decay models


decay.model<-function(y, x, model.type="exponential", y.type="similarities", perm=100, st.val = c(0.005, -0.00001)){

model.type <- match.arg(model.type, c("exponential", "power"))
y.type <- match.arg(y.type, c("similarities", "dissimilarities"))

switch(y.type, similarities = {ydis<-y
}, dissimilarities = {
ydis<-(1-y)
})

switch(model.type, exponential = {xdis<-x
}, power = {
xdis<-log(x)
})

nsites <- dim(as.matrix(xdis))[1]  
  
perm.dev.R2.distr <- numeric(perm) 
                  
    
for (i in 1:perm){ 
perm.sites <- sample(1:nsites)
ydismat <- as.matrix(ydis)                                                         
perm.ydis <- ydismat[perm.sites, perm.sites][lower.tri(ydismat)]
mod.perm <- glm(perm.ydis ~ xdis, family = gaussian(link = "log"), start = st.val)
perm.dev.R2.distr[i] <- 1-mod.perm$deviance/mod.perm$null.deviance
}

log.glm <- glm(ydis ~ xdis, family = gaussian(link = "log"), start = st.val)

pseudo.r.squared <- 1-log.glm$deviance/log.glm$null.deviance        
    
p.value <- mean(perm.dev.R2.distr > pseudo.r.squared )   


switch(y.type, similarities = {
switch(model.type, exponential = {
	parameters <- list(data.y = y, data.x = x, model = log.glm, model.type = model.type, y.type = y.type, 
	pseudo.r.squared = 1-log.glm$deviance/log.glm$null.deviance, 
	a.intercept = exp(log.glm$coefficients[1]), b.slope = log.glm$coefficients[2], 
	p.value = ifelse(p.value==0,1/perm,p.value))
}, 

power = {
parameters <- list(data.y = y, data.x = x, model = log.glm, model.type = model.type, y.type = y.type, 
	pseudo.r.squared = 1-log.glm$deviance/log.glm$null.deviance, 
	a.intercept = exp(log.glm$coefficients[1]), b.slope = log.glm$coefficients[2], 
	p.value = ifelse(p.value==0,1/perm,p.value))


})	
	
}, 

dissimilarities = {

switch(model.type, exponential = {
parameters <- list(data.y = y, data.x = x, model = log.glm, model.type = model.type, y.type = y.type, 
	pseudo.r.squared = 1-log.glm$deviance/log.glm$null.deviance, 
	a.intercept = 1-exp(log.glm$coefficients[1]), b.slope = -log.glm$coefficients[2], 
	p.value = ifelse(p.value==0,1/perm,p.value))
}, 

power = {
parameters <- list(data.y = y, data.x = x, model = log.glm, model.type = model.type, y.type = y.type, 
	pseudo.r.squared = 1-log.glm$deviance/log.glm$null.deviance, 
	a.intercept = 1-exp(log.glm$coefficients[1]), b.slope = -log.glm$coefficients[2], 
	p.value = ifelse(p.value==0,1/perm,p.value))
})
})

class(parameters)<-"decay"	
return(parameters)

}



################################################################################




################################################################################
##### A function for plotting distance-decay models

plot.decay<-function(x, xlim=c(0,max(x$data.x)), ylim=c(0,1), add=FALSE, remove.dots=FALSE, 
col="black", pch=1, lty=1, lwd=5, cex=1, ...){
if (!inherits(x, "decay")){	
		stop("The input is not a distance-decay model fitted with decay.model().",call.=TRUE)
	}
data<-data.frame(as.vector(x$data.x), as.vector(x$data.y))
if(!remove.dots){pch=pch}
else{pch=""}
dista<-sort(unique(data[,1]))
model.type <- match.arg(x$model.type, c("exponential", "power"))
y.type <- match.arg(x$y.type, c("similarities", "dissimilarities"))

switch(model.type, exponential = {
if(!add){
switch(y.type, similarities = {
plot(data, xlim=xlim, ylim=ylim, xlab="Distance", ylab="Similarity", col=col, pch=pch, cex=cex, ...)
lines(dista, x$a.intercept*exp(x$b.slope*dista), col=col, lty=lty, lwd=lwd, ...)
}, dissimilarities = {
plot(data, xlim=xlim, ylim=ylim, xlab="Distance", ylab="Dissimilarity", col=col, pch=pch, cex=cex, ...)
lines(dista, 1-(1-x$a.intercept)*exp(-x$b.slope*dista), col=col, lty=lty, lwd=lwd, ...)
})
}
if(add){
switch(y.type, similarities = {
points(data, xlim=xlim, ylim=ylim, xlab="Distance", ylab="Similarity", col=col, pch=pch, cex=cex, ...)
lines(dista, x$a.intercept*exp(x$b.slope*dista), col=col, lty=lty, lwd=lwd, ...)
}, dissimilarities = {
points(data, xlim=xlim, ylim=ylim, xlab="Distance", ylab="Dissimilarity", col=col, pch=pch, cex=cex, ...)
lines(dista, 1-(1-x$a.intercept)*exp(-x$b.slope*dista), col=col, lty=lty, lwd=lwd, ...)
})
}
}
, power = {
if(!add){
switch(y.type, similarities = {
plot(data, xlim=xlim, ylim=ylim, xlab="Distance", ylab="Similarity", col=col, pch=pch, cex=cex, ...)
lines(dista, x$a.intercept*dista^x$b.slope, col=col, lty=lty, lwd=lwd, ...)
}, dissimilarities = {
plot(data, xlim=xlim, ylim=ylim, xlab="Distance", ylab="Dissimilarity", col=col, pch=pch, cex=cex, ...)
lines(dista, 1-(1-x$a.intercept)*dista^-x$b.slope, col=col, lty=lty, lwd=lwd, ...)
})
}
if(add){
switch(y.type, similarities = {
points(data, xlim=xlim, ylim=ylim, xlab="Distance", ylab="Similarity", col=col, pch=pch, cex=cex, ...)
lines(dista, x$a.intercept*dista^x$b.slope, col=col, lty=lty, lwd=lwd, ...)
}, dissimilarities = {
points(data, xlim=xlim, ylim=ylim, xlab="Distance", ylab="Dissimilarity", col=col, pch=pch, cex=cex, ...)
lines(dista, 1-(1-x$a.intercept)*dista^-x$b.slope, col=col, lty=lty, lwd=lwd, ...)
})
}
})
}

################################################################################




################################################################################
##### A function for bootstrapping the parameters of distance-decay models


boot.coefs.decay<-function(x, resamples, st.val = c(0.005, -0.00001)){
if (!inherits(x, "decay")){	
		stop("The input is not a distance-decay model fitted with decay.model().",call.=TRUE)
	}
ptm <- proc.time()
data<-data.frame(as.vector(x$data.x), as.vector(x$data.y))
original.coefs<-c(x$a.intercept, x$b.slope)
names(original.coefs)<-c("a.intercept", "b.slope")
boot.coefs<-matrix(nrow=resamples, ncol=2)
colnames(boot.coefs)<-c("a.intercept", "b.slope")

for (i in 1:resamples){
data.sample<-data[sample(1:nrow(data), replace=TRUE), ]

switch(x$y.type, similarities = {ydis<-data.sample[,2]
}, 
dissimilarities = {
ydis<-(1-data.sample[,2])
})

switch(x$model.type, exponential = {xdis<-data.sample[,1]
}, power = {
xdis<-log(data.sample[,1])
})


boot.mod <- glm(ydis ~ xdis, family = gaussian(link = "log"), start = st.val)


switch(x$y.type, similarities = {
boot.coefs[i,]<-c(exp(boot.mod$coefficients[1]), boot.mod$coefficients[2])},

x$y.type,dissimilarities = {
boot.coefs[i,]<-c(1-exp(boot.mod$coefficients[1]), -boot.mod$coefficients[2])})

}


mean.boot<-c(mean(boot.coefs[,1]), mean(boot.coefs[,2]))
names(mean.boot)<-c("a.intercept", "b.slope")
sd.boot=c(sd(boot.coefs[,1]), sd(boot.coefs[,2]))
names(sd.boot)<-c("a.intercept", "b.slope")

result<-list( model.type=x$model.type,y.type=x$y.type, boot.coefs=boot.coefs, original.coefs=original.coefs, mean.boot=mean.boot, sd.boot=sd.boot,time=proc.time() - ptm)
result
}

################################################################################
