\name{full_matching}
\alias{full_matching}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Matching case and control }
\description{
  The function matches case and control based on distances as measured 
by the significant 
eigenvectors. Internal function.
}
\usage{
full_matching(ext)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ext}{the file extension of the distance matrix}
}
\details{
  The function calls the full matching program in optmatch library.
}
\value{
  \item{ma}{A matrix that consists of sampleId, matched strata, and 
case/control status. }
}
\keyword{models}
