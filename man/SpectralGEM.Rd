\name{SpectralGEM}
\alias{SpectralGEM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Software for Matching}
\description{
  SpectralGEM is designed to find the ancestry vectors and match cases and
controls for association analysis.
}
\usage{
SpectralGEM(InputFile = "matching_input.txt",CM="CM",outlier=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{InputFile}{The name of the file that contains the input parameters.}
  \item{CM}{options are C, M, and CM. C: for clustering only with 
removing outliers. M: for matching only without examing the outliers. CM: 
for clustering, removing outliers, and matching.}
  \item{outlier}{An option to remove outliers by checking the distributions of
    the distances between cases and contrls. Only applicable when CM="C"}
}
\details{
  %~~ If necessary, more details than the description above ~~
}
\value{
  \item{cl}{a two column matrix: the first column is the sample ID, the 
    second column is the cluster id}
  \item{U}{a matrix, the first column is the sample ID, the second
    columns is group id, the third column is the trivial eigenvector U0, 
and rest are the 
significant eigenvectors}
  \item{lambda}{eigenvalues corresponding to the eigenvectors}
  \item{d}{the distance between case and control}
  The program performs clustering and matching or matching only. The cl 
values are generated at the clustering stage. The significant 
eigenvectors are generated at the matching stage.

  A series of files are produced in the current directory.   
}
\references{ http://wpicr.wpic.pitt.edu/WPICCompGen/Spectral-GEM/GEM+.htm}
\author{Ann Lee, Diana Luca, Bert Klei, Bernie Devlin, and Kathryn
  Roeder
  
       Maintainer: Jing Wu jwu@stat.cmu.edu}
\note{ The function depends on the local fortran executable. The function 
asks whether the user would like to download the executable before it  
automatically downloads the executable from the reference website.} 
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
data(H)
data(id.info)
MMfile(H=H,sampleInfo=id.info,ntag=1000,outfile="MMprime.txt")
InFile(identifier="smal",stage=1,directory="./",
                    MMfile="MMprime.txt",excludefile="exclude.txt",
                    idlength=8,mincluster=10,logtype=0,
                    outfile="matching_input.txt")

#not run#
#out=SpectralGEM() #first do clustering and remove outliers
                                 #then do matching
#out=SpectralGEM(CM="C") # do clustering and remove 
                                         #outliers
#out=SpectralGEM(CM="M") # do matching without removing 
                                         #outliers

# For continuous response, create new id.info for the H matrix 
data(H)
n=345
y=sample(c(rnorm(mean=0,172),rnorm(mean=1,173)))
cc=rep(1,345)
cc[y>median(y)]=2  #create case control
newid.info=cbind(c(1:n),rep(1,n),cc,cc)
MMfile(H=H,sampleInfo=newid.info,ntag=1000,outfile="MMprime1.txt")
InFile(identifier="smal",stage=1,directory="./",
                    MMfile="MMprime1.txt",excludefile="exclude.txt",
                    idlength=8,mincluster=10,logtype=0,
                    outfile="input.txt")
#not run#
#out=SpectralGEM(InputFile="input.txt",CM="C") 

# buildin ancestry plots
#Current version: ext= smal1 ;
#not run
#pc_graphs_GEMpClusters("smal1")
#pc_graphs_GEMp("smal1")

#plot from the SpectralGEM output, 
#significant eigenvectors start from out$U[,4] with
#significant eigenvalues start from from out$lambda[2]
#not run
#plot(sqrt(out$lambda[2])*out$U[,4],sqrt(out$lambda[3])*out$U[,5],
#col=out$U[,2],xlab="PC 1",ylab="PC 2") #for PC plots
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}
\keyword{models}
\keyword{regression}
% __ONLY ONE__ keyword per line
