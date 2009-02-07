\name{InFile}
\alias{InFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate a text file}
\description{
The function generates a text file that contains the parameters needed in 
the main fortran program and saves it in the current directory.
}
\usage{
InFile(identifier="smal",stage="1",directory="./",
                    MMfile="MMprime.txt",excludefile="exclude.txt",
                    idlength=8,mincluster=10,logtype=0,
                    outfile="matching_input.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{identifier}{must be 4-letters long}
  \item{stage}{must be 1-letter long}
  \item{directory}{the directory where the files will be generated}
  \item{MMfile}{the MMprime matrix file name}
  \item{excludefile}{the name of the file where the outliers have been
    or will be written to}
  \item{idlength}{the number of letters the longest sample name has}
  \item{mincluster}{the smallest cluster required when do clustering}
  \item{logtype}{the type of log information output to the log file:
    0=limited, 2=lots}
  \item{outfile}{The name of the output file.}
}
\details{
 The function generates a text file with one column which contains all the 
parameters needed. The parameters need to be in the right order in order 
for the main fortran program to use.  
}
\value{
 A text file is saved in the current directory.
}
\references{http://wpicr.wpic.pitt.edu/WPICCompGen/Spectral-GEM/GEM+.htm }
\author{Ann Lee, Diana Luca, Bert Klei, Bernie Devlin, and Kathryn
  Roeder
  
        Maintainer: Jing Wu jwu@stat.cmu.edu}
\seealso{http://wpicr.wpic.pitt.edu/WPICCompGen/Spectral-GEM/matching\_input.txt}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ file }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
