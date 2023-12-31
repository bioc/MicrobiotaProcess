% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/method-mp_cal_clust.R
\name{get_clust}
\alias{get_clust}
\alias{get_clust.dist}
\alias{get_clust.data.frame}
\alias{get_clust.phyloseq}
\title{Hierarchical cluster analysis for the samples}
\usage{
get_clust(obj, ...)

\method{get_clust}{dist}(obj, distmethod, sampleda = NULL, hclustmethod = "average", ...)

\method{get_clust}{data.frame}(
  obj,
  distmethod = "euclidean",
  taxa_are_rows = FALSE,
  sampleda = NULL,
  tree = NULL,
  method = "hellinger",
  hclustmethod = "average",
  ...
)

\method{get_clust}{phyloseq}(
  obj,
  distmethod = "euclidean",
  method = "hellinger",
  hclustmethod = "average",
  ...
)
}
\arguments{
\item{obj}{phyloseq, phyloseq class or dist class, or 
data.frame, data.frame, default is nrow samples * ncol features.}

\item{..., }{additional parameters.}

\item{distmethod}{character, the method of dist, when the 
obj is data.frame or phyloseq default is "euclidean". see also 
\code{\link[MicrobiotaProcess]{get_dist}}.}

\item{sampleda}{data.frame, nrow sample * ncol factor. default is NULL.}

\item{hclustmethod}{character, the method of hierarchical cluster, 
default is average.}

\item{taxa_are_rows}{logical, if the features of data.frame(obj) 
is in column, it should set FALSE.}

\item{tree}{phylo, the phylo class, see also \code{\link[ape]{as.phylo}}.}

\item{method}{character, the standardization methods for community 
ecologists, see also \code{\link[vegan]{decostand}}}
}
\value{
treedata object.
}
\description{
Hierarchical cluster analysis for the samples
}
\examples{
\dontrun{
library(phyloseq)
data(GlobalPatterns)
subGlobal <- subset_samples(GlobalPatterns, 
         SampleType \%in\% c("Feces", "Mock", "Ocean", "Skin"))
hcsample <- get_clust(subGlobal, distmethod="jaccard",
                  method="hellinger", hclustmethod="average")
}
}
\author{
Shuangbin Xu
}
