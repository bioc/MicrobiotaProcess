% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-methods.R
\name{mp_plot_venn}
\alias{mp_plot_venn}
\alias{mp_plot_venn,MPSE-method}
\alias{mp_plot_venn,MPSE}
\alias{mp_plot_venn,tbl_mpse-method}
\alias{mp_plot_venn,tbl_mpse}
\alias{mp_plot_venn,grouped_df_mpse-method}
\alias{mp_plot_venn,grouped_df_mpse}
\title{Plotting the different number of OTU between groups with Venn Diagram.}
\usage{
mp_plot_venn(.data, .group, .venn = NULL, ...)

\S4method{mp_plot_venn}{MPSE}(.data, .group, .venn = NULL, ...)

\S4method{mp_plot_venn}{tbl_mpse}(.data, .group, .venn = NULL, ...)

\S4method{mp_plot_venn}{grouped_df_mpse}(.data, .group, .venn = NULL, ...)
}
\arguments{
\item{.data}{MPSE object or tbl_mpse object}

\item{.group}{the column names of group to be visualized}

\item{.venn}{the column names of result after run \code{mp_cal_venn}.}

\item{...}{additional parameters, such as 'size', 'label_size', 'edge_size' etc, 
see also 'ggVennDiagram'.}
}
\description{
Plotting the different number of OTU between groups with Venn Diagram.
}
\examples{
\dontrun{
data(mouse.time.mpse)
mpse <- mouse.time.mpse \%>\%
        mp_rrarefy() \%>\%
        mp_cal_venn(.abundance=RareAbundance, .group=time, action="add")
mpse
p <- mpse \%>\% mp_plot_venn(.group=time, .venn=vennOftime) 
p
}
}
\author{
Shuangbin Xu
}
