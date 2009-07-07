\name{MMfile}
\alias{MMfile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Produce a .txt file that contains the MM prime matrix.}
\description{
 The function produces a .txt file of the proper format that contains the 
input data and the adjaceny matrix. Because a common choice for the latter 
is MM', this is called the MM prime matrix by default. This 
file is 
loaded to the main fortran program to do eigenvalue decomposition. 
Alternatively, the 
output file can be created from an input file containing M by the fortran 
program located at
http://wpicr.wpic.pitt.edu/WPICCompGen/MMp/MMp\_page.htm.
}
\usage{
MMfile(H = H,sampleInfo = id.info,n=dim(H)[1],ntag=ntag,outfile = "MMprime.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{H}{A square symmetirc matrix.}
  \item{sampleInfo}{A 4-column matrix containing sample\_id, sex,
    case\_control status, and group\_id}
  \item{n}{Number of individual in the H matrix}
  \item{ntag}{Number of tag SNPs used to generate the H matrix.}
  \item{outfile}{The text file to write to. The default is set as
MMprime.txt.}
}
\details{
  The first two lines in the output file contains the information about 
the  square matrix. The first line is the number of individuals in the 
matrix. The second line is the number of tag SNPs.
  The rest of the file is a table. Each line in the table contains the 
sample ID, gender, case/control status, group ID, and the square matrix.
}
\value{
 A text file is produced in the current directory.
}
\references{http://wpicr.wpic.pitt.edu/WPICCompGen/Spectral-GEM/directions.pdf}
\author{Ann Lee, Diana Luca, Bert Klei, Bernie Devlin, and Kathryn
  Roeder
  
        Maintainer: Jing Wu jwu@stat.cmu.edu}
\seealso{http://wpicr.wpic.pitt.edu/WPICCompGen/MMp/MMp\_page.htm}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{file}

