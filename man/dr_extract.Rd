% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/method-dr-tibble.R
\name{dr_extract}
\alias{dr_extract}
\title{Extracting the internal tbl_df attribute of tibble.}
\usage{
dr_extract(name, .f = NULL)
}
\arguments{
\item{name}{character the name of internal tbl_df attribute.}

\item{.f}{a function (if any, default is NULL) that pre-operate the
data}
}
\value{
tbl_df object
}
\description{
Extracting the internal tbl_df attribute of tibble.
}
\examples{
\dontrun{
library(vegan)
data(varespec, varechem)
mpse <- MPSE(assays=list(Abundance=t(varespec)), colData=varechem)
tbl <- 
mpse \%>\%
  mp_cal_nmds(.abundance=Abundance, action="add") \%>\%
  mp_envfit(.ord=NMDS, .env=colnames(varechem), action="only") 
tbl 
tbl \%>\% attributes \%>\% names
# This function is useful to extract the data to display with ggplot2
# you can also refer to the examples of mp_envfit.
dr_extract(name=NMDS_ENVFIT_tb)(tbl)
# add .f function 
dr_extract(name=NMDS_ENVFIT_tb, 
           .f=td_filter(pvals<=0.05 & label!="Humdepth"))(tbl)
}
}
\author{
Shuangbin Xu
}
