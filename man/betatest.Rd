\encoding{utf8}

\name{betatest}
\alias{betatest}
\docType{data}
\title{A data set of 4 communities, 107 species and a 4D functional space}
\description{
A data set to test the functions to dissimilarities matricies.
}
\usage{data("betatest")}
\format{
  \code{betatest} is a list of two elements.
  \describe{
    \item{\code{comm.test}}{ a dataframe with the presence/absence of 107 species among 4 sites}
    \item{\code{traits.test}}{ a dataframe of the traits (4 axes) of the 107 species}
  }
}
\examples{
data(betatest)
str(betatest$comm.test)
str(betatest$traits.test)
}
\keyword{datasets}
