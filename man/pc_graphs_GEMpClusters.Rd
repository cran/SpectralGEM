\name{pc_graphs_GEMpClusters}
\alias{pc_graphs_GEMpClusters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot ancestral structures with cluster id}
\description{
  The function generates ancestrol plots in .pdf and .ps form.
}
\usage{
pc_graphs_GEMpClusters(ext)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ext}{the identifier and stage concatenated}
}
\details{
The function plots pairs of principle componets corresponding to the
significant eigenvectors. Use only after the clustering step.
}
\value{
%  ~Describe the value returned
%  If it is a LIST, use
%  \item{comp1 }{Description of 'comp1'}
%  \item{comp2 }{Description of 'comp2'}
%  ...
  The plots are put in the current folder.
}
\references{ 
http://wpicr.wpic.pitt.edu/WPICCompGen/Spectral-GEM/GEM+.htm}
\author{Ann Lee, Diana Luca, Bert Klei, Bernie Devlin, and 
Kathryn  Roeder

       Maintainer: Jing Wu jwu@stat.cmu.edu}
\examples{
#   pc_graphs_GEMpClusters("smal1")
}
% R documentation directory.
\keyword{dplot}

