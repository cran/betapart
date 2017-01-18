betapart.core.abund <- function(x){
	
	# test for a numeric matrix or data.frame
	if(! is.matrix(x)){
		x<-as.matrix(x)
  	}
	
	if(! is.numeric(x))
    	stop("The data in x is not numeric.",call.=TRUE)
	
	# simple test for positive data
	xvals <-  unique(as.vector(x))
	if (any(xvals<0)) 
        stop("The table contains negative values: data should be abundances.", call. = TRUE)

	pair.shared.abund<-matrix(nrow=nrow(x),ncol=nrow(x))
	rownames(pair.shared.abund)<-rownames(x)
	colnames(pair.shared.abund)<-rownames(x)

	pair.not.shared.abund<-matrix(nrow=nrow(x),ncol=nrow(x))
	rownames(pair.not.shared.abund)<-rownames(x)
	colnames(pair.not.shared.abund)<-rownames(x)

	for(i in 1:nrow(x)) {
    	for(j in i:nrow(x)) {
		pair.shared.abund[j,i]<-sum(pmin(x[i,],x[j,]))
		pair.not.shared.abund[i,j]<-sum(x[i,])-sum(pmin(x[i,],x[j,]))
		pair.not.shared.abund[j,i]<-sum(x[j,])-sum(pmin(x[i,],x[j,]))
		}
		}

	pair.shared.abund<-as.dist(pair.shared.abund)
	pair.max.not.shared.abund<-pmax(as.dist(t(upper.tri(pair.not.shared.abund)*pair.not.shared.abund)), as.dist(pair.not.shared.abund))
	pair.min.not.shared.abund<-pmin(as.dist(t(upper.tri(pair.not.shared.abund)*pair.not.shared.abund)), as.dist(pair.not.shared.abund))

	multiple.shared.abund<-sum(x)-sum(apply(x,2,max))



	computations<-list(data=x, multiple.shared.abund=multiple.shared.abund, pair.shared.abund=pair.shared.abund, pair.min.not.shared.abund=pair.min.not.shared.abund, 
		pair.max.not.shared.abund=pair.max.not.shared.abund, pair.sum.not.shared.abund=pair.min.not.shared.abund+pair.max.not.shared.abund)
    class(computations)<-"betapart.abund"

	return(computations)
} 
