% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/method-accessor.R
\name{mp_extract_rarecurve}
\alias{mp_extract_rarecurve}
\alias{mp_extract_rarecurve,MPSE-method}
\alias{mp_extract_rarecurve,MPSE}
\alias{mp_extract_rarecurve,tbl_mpse-method}
\alias{mp_extract_rarecurve,tbl_mpse}
\alias{mp_extract_rarecurve,grouped_df_mpse-method}
\alias{mp_extract_rarecurve,grouped_df_mpse}
\title{Extract the result of mp_cal_rarecurve with action="add" from MPSE or tbl_mpse object}
\usage{
mp_extract_rarecurve(x, .rarecurve, ...)

\S4method{mp_extract_rarecurve}{MPSE}(x, .rarecurve, ...)

\S4method{mp_extract_rarecurve}{tbl_mpse}(x, .rarecurve, ...)

\S4method{mp_extract_rarecurve}{grouped_df_mpse}(x, .rarecurve, ...)
}
\arguments{
\item{x}{MPSE object or tbl_mpse object}

\item{.rarecurve}{the column name of rarecurve after run mp_cal_rarecurve with action="add".}

\item{...}{additional parameter}
}
\value{
rarecurve object that be be visualized by ggrarecurve
}
\description{
Extract the result of mp_cal_rarecurve with action="add" from MPSE or tbl_mpse object
}
