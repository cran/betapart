\name{bbsData}
\alias{bbs1980}
\alias{bbs2000}
\docType{data}
\title{BBS data by state for two timeslices}
\description{
The data consists of binary presence/absence matrices for 569 bird species across 49 US states for two time slices (1980 - 1985 and 2000 - 2005). Only species (identified by AOU number) recorded during both time periods are included. The data are taken from the North American Breeding Bird Survey dataset and from a version of the database downloaded in May 2009.
}
\usage{data(bbsData)}
\format{
  Two matrices (bbs1980 and bbs2000) of identical structure showing the presence/absence of the species as binary data.
  \describe{
    \item{\code{statenum}}{US states by number.}
    \item{\code{aou}}{Species idenity by AOU ID number.}
  }
}

\source{
http://www.pwrc.usgs.gov/BBS/
}
\examples{
data(bbsData)
str(bbs1980)
str(bbs2000)
}
\keyword{datasets}
