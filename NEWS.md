# betapart 1.6

## New set of functions for distance-decay models:

* `decay.model()` now allows fitting the Gompertz function (as described in Martín-Devasa et al. 2022, Diversity and Distributions), in addition to negative exponential and power law functions. To do this, nonlinear models are now fitted via the `nls.lm()` function in the `minpack.lm package` (which uses the Levenberg-Marquardt Nonlinear Least-Squares Algorithm), instead of `glm()` as in the previous versions of `betapart`

* `boot.coefs.decay()` now implements the site-block resampling method introduced in Martínez-Santalla et al. 2022 (J. Biogeogr.)

* A new function, `zdep()` allows assessing the significance of differences between parameters of two distance-decay models, as introduced in Martín-Devasa et al. 2022 (Ecol. Informatics)

* `plot.decay()` was modified to handle the new `decay.model()` function


# betapart 1.5.6

* Modified function `decay.model()` to implement the block permutation described in Martínez-Santalla et al 2022 (J. Biogeogr.) for assessing the significance of the distance decay model


# betapart 1.5.5

* A bug of `functional.betapart.core()` was fixed to run it in parallel with multi = TRUE


# betapart 1.5.4

* Updated `functional.betapart.core.pairwise()` to get vertice coordinates in the output `details`


# betapart 1.5.3

## The libraries doParallel/parallel were replaced by doSNOW/snow for parallel computations.

* `functional.betapart.core()` was updated. Options can be passed to qhull to prevent some crashes and a progress bar can be displayed. When setting multi=TRUE, the function stop earlier if the number of communities is too important.

* New function to control options passed to qhull for convexhull estimation: `qhull.opt()`

* New function to compute rapidly pair-wise dissimilarty matrices: `functional.betapart.core.pairwise()` 

* `functional.beta.pair()` was updated to integrate functional.betapart.core.pairwise


# betapart 1.5.2

## New features:

* Updated `functional.betapart.core()` to allow internal parallel computing
* New function to customize parameters for the internal parallel computing : `beta.para.control()`


# betapart 1.5.1

## New features:

* Updated `functional.betapart.core()` to allow parallel computing


# betapart 1.5.0

### New functions to fit, plot and bootstrap distance-decay patterns

* `decay.model()` fits a negative-exponential or mower law function describing the decay of assemblage similarity with spatial distance.

* `plot.decay()` allows plotting the curves fitted with `decay.model()`.

* `boot.coefs.decay()` bootstraps the parameters of the functions fitted with `decay.model()`.

