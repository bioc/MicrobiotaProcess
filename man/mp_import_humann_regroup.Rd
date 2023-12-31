% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/method-io.R
\name{mp_import_humann_regroup}
\alias{mp_import_humann_regroup}
\title{Import function to load the output of human_regroup_table in HUMAnN.}
\usage{
mp_import_humann_regroup(
  profile,
  mapfilename = NULL,
  rm.unknown = TRUE,
  keep.contribute.abundance = FALSE,
  ...
)
}
\arguments{
\item{profile}{the output file (text format) of human_regroup_table in HUMAnN.}

\item{mapfilename}{the sample information file or data.frame,}

\item{rm.unknown}{logical whether remove the unmapped and ungrouped features.}

\item{keep.contribute.abundance}{logical whether keep the abundance of contributed taxa,
default is FALSE, it will consume more memory if it set to TRUE.}

\item{...}{additional parameters, meaningless now.}
}
\description{
Import function to load the output of human_regroup_table in HUMAnN.
}
\author{
Shuangbin Xu
}
