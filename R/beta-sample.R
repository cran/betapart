beta.sample <- function(x, index.family="sorensen", sites=10, samples=100){
		
	# test for a valid index
	index.family <- match.arg(index.family, c('jaccard','sorensen'))

	# test for pre-existing betapart objects - validates data in x if not
	if (! inherits(x, "betapart")){	
		x <- betapart.core(x)
	}
	
	# check sensible inputs
	if(sites > nrow(x$data))
		stop('More sites requested for sample than are in the dataset')
	
	pb <- txtProgressBar(min = 0, max = samples, style = 3)

	# Create a matrix to save the results.
	results.n<-as.data.frame(matrix(nrow=samples, ncol=3))

	# Loop on the selected number of samples
	for(i in 1:samples){

  		# Take a sample of the dataset with the specified number of sites 
  		position <- as.vector(1:nrow(x$data))
  		sample.position <- sample(position, sites)
  		x.sample <- x$data[sample.position, ]
  
  		# Compute the three indices for this sample and save in the results matrix
		x.beta <- beta.multi(x.sample, index.family)	
  		results.n[i,] <- unlist(x.beta)
		
   		# update progress bar
   		setTxtProgressBar(pb, i)         
	}
	
	close(pb)
		
	# give the results the right column names from the relevant family index
	names(results.n) <- names(x.beta)

	result <- list(sampled.values=results.n, mean.values=sapply(results.n, mean),
	               sd.values=sapply(results.n, sd))
	return(result)
}
