# betapart 1.5.4

### Updated functional.betapart.core.pairwise() to get vertice coordinates in the output 'details'

# betapart 1.5.3

* the libraries doParallel/parallel were replaced by doSNOW/snow for parallel computations.

## new features
### functional.betapart.core() was updated. Options can be passed to qhull to prevent some crashes and a progress bar can be displayed.
### When setting multi=TRUE, the function stop earlier if the number of communities is too important.

### New fonction to control options passed to qhull for convexhull estimation:
* `qhull.opt()`

### New function to compute rapidly pair-wise dissimilarty matrices:
* `functional.betapart.core.pairwise()` 

### functional.beta.pair was updated to integrate functional.betapart.core.pairwise


# betapart 1.5.2

## New features

### Updated functional.betapart.core() to allow internal parallel computing
### New function to customize parameters for the internal parallel computing :
* `beta.para.control()`

# betapart 1.5.1

## New features

### Updated functional.betapart.core() to allow parallel computing


# betapart 1.5.0

## New features

### New functions to fit, plot and bootstrap distance-decay patterns

`betapart 1.5.0` includes three new functions:

* `decay.model()` fits a negative-exponential or mower law function describing the decay of assemblage similarity with sptatial distance.

* `plot.decay()` allows plotting the curves fitted with `decay.model()`.

* `boot.coefs.decay()` bootstraps the parameters of the functions fitted with `decay.model()`.

