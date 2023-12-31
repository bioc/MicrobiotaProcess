% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/method-io.R
\name{mp_import_biom}
\alias{mp_import_biom}
\title{building MPSE object from biom-format file.}
\usage{
mp_import_biom(
  biomfilename,
  mapfilename = NULL,
  otutree = NULL,
  refseq = NULL,
  ...
)
}
\arguments{
\item{biomfilename}{character the biom-format file path.}

\item{mapfilename}{character, the file contained sample information,
the tsv format, default is NULL.}

\item{otutree}{treedata, phylo or character, the file contained reference sequences, or
treedata object, which is the result parsed by functions of treeio, default is NULL.}

\item{refseq}{XStringSet or character, the file contained the representation sequence file or 
XStringSet class to store the representation sequence, default is NULL.}

\item{...}{additional parameter, which is meaningless now.}
}
\value{
MPSE-class
}
\description{
building MPSE object from biom-format file.
}
