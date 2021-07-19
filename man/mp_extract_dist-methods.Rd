% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/method-accessor.R
\docType{methods}
\name{mp_extract_dist}
\alias{mp_extract_dist}
\alias{mp_extract_dist,MPSE-method}
\alias{mp_extract_dist,MPSE}
\alias{mp_extract_dist,tbl_mpse-method}
\alias{mp_extract_dist,tbl_mpse}
\alias{mp_extract_dist,grouped_df_mpse-method}
\alias{mp_extract_dist,grouped_df_mpse}
\title{extract the dist object from MPSE or tbl_mpse object}
\usage{
mp_extract_dist(x, distmethod, env = FALSE)

\S4method{mp_extract_dist}{MPSE}(x, distmethod, env = FALSE)

\S4method{mp_extract_dist}{tbl_mpse}(x, distmethod, env = FALSE)

\S4method{mp_extract_dist}{grouped_df_mpse}(x, distmethod, env = FALSE)
}
\arguments{
\item{x}{MPSE object or tbl_mpse object}

\item{distmethod}{character the method of calculated distance.}

\item{env}{logical whether extract the distance of samples calculated 
based on continuous environment factors, default is FALSE.}
}
\value{
dist object.
}
\description{
extract the dist object from MPSE or tbl_mpse object
}