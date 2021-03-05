\name{beta.para.control}
\alias{beta.para.control}
\encoding{utf8}

\title{
Specifying Control Values for Internal Parallel Cluster
}
\description{
The values supplied in the \code{beta.para.control()} call replace the defaults, and a list with all settings (i.e., values for all possible arguments) is returned. The returned list is used to define the internal parallel cluster of the \code{functional.betapart.core} function.

}
\usage{
beta.para.control(nc = floor(parallel::detectCores()/2), type = "SOCK", 
                              LB = TRUE, size = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{nc}{ number of cores to use. Default is half of the available cores.}
\item{type}{ character - the type of cluster to be used, either \code{"SOCK"}, \code{"PSOCK"} or \code{"FORK"} (not on Windows).}
\item{LB}{ logical indicating if load balancing has to be used. Default is \code{TRUE}}
\item{size}{ number of operation run on each core at each iteration. Default is \code{1}.}
}

\value{
a list with components for each of the possible arguments.
}
\author{
Maxime Logez
}
\examples{
str(beta.para.control(nc = 2, LB = FALSE))
}
