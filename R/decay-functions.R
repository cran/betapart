########################################################################################
##### A function for fitting a distance-decay models


decay.model<-function(y, x, model.type="exponential", y.type="similarities", perm=100){
model.type <- match.arg(model.type, c("exponential", "power"))

switch(model.type, exponential = {

y.type <- match.arg(y.type, c("similarities", "dissimilarities"))

switch(y.type, similarities = {
y<-as.vector(y)
x<-as.vector(x)
log.glm<-glm(y~x, family=gaussian(link = "log"))

null.dev<-vector(mode="numeric", length=perm)
for (i in 1:perm){
null.dev[i]<-glm(y~sample(x), family=gaussian(link = "log"))$deviance
}
p.value<-mean(null.dev<log.glm$deviance)

parameters <- list(data = data.frame(x,y), model = log.glm, model.type = model.type, y.type = y.type, 
	pseudo.r.squared = 1-log.glm$deviance/log.glm$null.deviance, 
	a.intercept = exp(log.glm$coefficients[1]), b.slope = log.glm$coefficients[2], 
	p.value = ifelse(p.value==0,1/perm,p.value))
class(parameters)<-"decay"	
return(parameters)

}, dissimilarities = {
y<-as.vector(1-y)
x<-as.vector(x)
log.glm<-glm(y~x, family=gaussian(link = "log"))

null.dev<-vector(mode="numeric", length=perm)
for (i in 1:perm){
null.dev[i]<-glm(y~sample(x), family=gaussian(link = "log"))$deviance
}
p.value<-mean(null.dev<log.glm$deviance)


parameters <- list(data = data.frame(x,1-y), model = log.glm, model.type = model.type, y.type = y.type, 
	pseudo.r.squared = 1-log.glm$deviance/log.glm$null.deviance, 
	a.intercept = 1-exp(log.glm$coefficients[1]), b.slope = -log.glm$coefficients[2], 
	p.value = ifelse(p.value==0,1/perm,p.value))
class(parameters)<-"decay"	
return(parameters)
})
}

, power = {

y.type <- match.arg(y.type, c("similarities", "dissimilarities"))


switch(y.type, similarities = {
y<-as.vector(y)
x<-as.vector(log(x))
log.glm<-glm(y~x, family=gaussian(link = "log"))

null.dev<-vector(mode="numeric", length=perm)
for (i in 1:perm){
null.dev[i]<-glm(y~sample(x), family=gaussian(link = "log"))$deviance
}
p.value<-mean(null.dev<log.glm$deviance)

parameters <- list(data = data.frame(exp(x),y), model = log.glm, model.type = model.type, y.type = y.type, 
	pseudo.r.squared = 1-log.glm$deviance/log.glm$null.deviance, 
	a.intercept = exp(log.glm$coefficients[1]), b.slope = log.glm$coefficients[2], 
	p.value = ifelse(p.value==0,1/perm,p.value))
class(parameters)<-"decay"	
return(parameters)

}, dissimilarities = {
y<-as.vector(1-y)
x<-as.vector(log(x))
log.glm<-glm(y~x, family=gaussian(link = "log"))

null.dev<-vector(mode="numeric", length=perm)
for (i in 1:perm){
null.dev[i]<-glm(y~sample(x), family=gaussian(link = "log"))$deviance
}
p.value<-mean(null.dev<log.glm$deviance)


parameters <- list(data = data.frame(exp(x),1-y), model = log.glm, model.type = model.type, y.type = y.type, 
	pseudo.r.squared = 1-log.glm$deviance/log.glm$null.deviance, 
	a.intercept = 1-exp(log.glm$coefficients[1]), b.slope = -log.glm$coefficients[2], 
	p.value = ifelse(p.value==0,1/perm,p.value))
class(parameters)<-"decay"	
return(parameters)

})
})
}


################################################################################




################################################################################
##### A function for plotting the decay models

plot.decay<-function(x, xlim=c(0,max(x$data[,1])), ylim=c(0,1), add=FALSE, remove.dots=FALSE, 
col="black", pch=1, lty=1, lwd=5, cex=1, ...){
if (!inherits(x, "decay")){	
		stop("The input is not a distance-decay model fitted with decay.model().",call.=TRUE)
	}

if(!remove.dots){pch=pch}
else{pch=""}
dista<-sort(unique(x$data[,1]))
model.type <- match.arg(x$model.type, c("exponential", "power"))
y.type <- match.arg(x$y.type, c("similarities", "dissimilarities"))

switch(model.type, exponential = {
if(!add){
switch(y.type, similarities = {
plot(x$data, xlim=xlim, ylim=ylim, xlab="Distance", ylab="Similarity", col=col, pch=pch, cex=cex, ...)
lines(dista, x$a.intercept*exp(x$b.slope*dista), col=col, lty=lty, lwd=lwd, ...)
}, dissimilarities = {
plot(x$data, xlim=xlim, ylim=ylim, xlab="Distance", ylab="Dissimilarity", col=col, pch=pch, cex=cex, ...)
lines(dista, 1-(1-x$a.intercept)*exp(-x$b.slope*dista), col=col, lty=lty, lwd=lwd, ...)
})
}
if(add){
switch(y.type, similarities = {
points(x$data, xlim=xlim, ylim=ylim, xlab="Distance", ylab="Similarity", col=col, pch=pch, cex=cex, ...)
lines(dista, x$a.intercept*exp(x$b.slope*dista), col=col, lty=lty, lwd=lwd, ...)
}, dissimilarities = {
points(x$data, xlim=xlim, ylim=ylim, xlab="Distance", ylab="Dissimilarity", col=col, pch=pch, cex=cex, ...)
lines(dista, 1-(1-x$a.intercept)*exp(-x$b.slope*dista), col=col, lty=lty, lwd=lwd, ...)
})
}
}
, power = {
if(!add){
switch(y.type, similarities = {
plot(x$data, xlim=xlim, ylim=ylim, xlab="Distance", ylab="Similarity", col=col, pch=pch, cex=cex, ...)
lines(dista, x$a.intercept*dista^x$b.slope, col=col, lty=lty, lwd=lwd, ...)
}, dissimilarities = {
plot(x$data, xlim=xlim, ylim=ylim, xlab="Distance", ylab="Dissimilarity", col=col, pch=pch, cex=cex, ...)
lines(dista, 1-(1-x$a.intercept)*dista^-x$b.slope, col=col, lty=lty, lwd=lwd, ...)
})
}
if(add){
switch(y.type, similarities = {
points(x$data, xlim=xlim, ylim=ylim, xlab="Distance", ylab="Similarity", col=col, pch=pch, cex=cex, ...)
lines(dista, x$a.intercept*dista^x$b.slope, col=col, lty=lty, lwd=lwd, ...)
}, dissimilarities = {
points(x$data, xlim=xlim, ylim=ylim, xlab="Distance", ylab="Dissimilarity", col=col, pch=pch, cex=cex, ...)
lines(dista, 1-(1-x$a.intercept)*dista^-x$b.slope, col=col, lty=lty, lwd=lwd, ...)
})
}
})
}

################################################################################




################################################################################
##### A function to bootstrap the parameters
##### It takes a decay model (x) and bootstrap the parameters R times

boot.coefs.decay<-function(x, R){
if (!inherits(x, "decay")){	
		stop("The input is not a distance-decay model fitted with decay.model().",call.=TRUE)
	}
ptm <- proc.time()
original.coefs<-c(x$a.intercept, x$b.slope)
names(original.coefs)<-c("a.intercept", "b.slope")
boot.coefs<-matrix(nrow=R, ncol=2)
colnames(boot.coefs)<-c("a.intercept", "b.slope")
for (i in 1:R){
data.sample<-x$data[sample(1:nrow(x$data), replace=TRUE), ]
boot.mod<-decay.model(data.sample[,2], data.sample[,1], model.type=x$model.type, y.type=x$y.type, perm=1)
boot.coefs[i,]<-c(boot.mod$a.intercept, boot.mod$b.slope)
}
mean.boot<-c(mean(boot.coefs[,1]), mean(boot.coefs[,2]))
names(mean.boot)<-c("a.intercept", "b.slope")
sd.boot=c(sd(boot.coefs[,1]), sd(boot.coefs[,2]))
names(sd.boot)<-c("a.intercept", "b.slope")

result<-list( model.type=x$model.type,y.type=x$y.type, boot.coefs=boot.coefs, original.coefs=original.coefs, mean.boot=mean.boot, sd.boot=sd.boot,time=proc.time() - ptm)
result
}

################################################################################
