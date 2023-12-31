% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-methods.R
\name{mp_plot_dist}
\alias{mp_plot_dist}
\alias{mp_plot_dist,MPSE-method}
\alias{mp_plot_dist,MPSE}
\alias{mp_plot_dist,tbl_mpse-method}
\alias{mp_plot_dist,tbl_mpse}
\alias{mp_plot_dist,grouped_df_mpse-method}
\alias{mp_plot_dist,grouped_df_mpse}
\title{Plotting the distance between the samples with heatmap or boxplot.}
\usage{
mp_plot_dist(
  .data,
  .distmethod,
  .group = NULL,
  group.test = FALSE,
  hclustmethod = "average",
  test = "wilcox.test",
  comparisons = NULL,
  step_increase = 0.1,
  ...
)

\S4method{mp_plot_dist}{MPSE}(
  .data,
  .distmethod,
  .group = NULL,
  group.test = FALSE,
  hclustmethod = "average",
  test = "wilcox.test",
  comparisons = NULL,
  step_increase = 0.1,
  ...
)

\S4method{mp_plot_dist}{tbl_mpse}(
  .data,
  .distmethod,
  .group = NULL,
  group.test = FALSE,
  hclustmethod = "average",
  test = "wilcox.test",
  comparisons = NULL,
  step_increase = 0.1,
  ...
)

\S4method{mp_plot_dist}{grouped_df_mpse}(
  .data,
  .distmethod,
  .group = NULL,
  group.test = FALSE,
  hclustmethod = "average",
  test = "wilcox.test",
  comparisons = NULL,
  step_increase = 0.1,
  ...
)
}
\arguments{
\item{.data}{the MPSE or tbl_mpse object after [mp_cal_dist()] is performed with action="add"}

\item{.distmethod}{the column names of distance of samples, it will generate after
[mp_cal_dist()] is performed.}

\item{.group}{the column names of group, default is NULL, when it is not provided 
the heatmap of distance between samples will be returned. If it is provided and
\code{group.test} is TURE, the comparisons boxplot of distance between the group
will be returned, but when \code{group.test} is FALSE, the heatmap of distance between
samples with group information will be returned.}

\item{group.test}{logical default is FALSE, see the \code{.group} argument.}

\item{hclustmethod}{character the method of \code{\link[stats]{hclust}}, default is
'average' (= UPGMA).}

\item{test}{the name of the statistical test, default is 'wilcox.test'}

\item{comparisons}{A list of length-2 vectors. The entries in the vector are
either the names of 2 values on the x-axis or the 2 integers that
correspond to the index of the columns of interest, default is NULL, meaning
it will be calculated automatically with the names in the .group.}

\item{step_increase}{numeric vector with the increase in fraction of total
height for every additional comparison to minimize overlap, default is 0.1.}

\item{...}{additional parameters, see also \code{\link[ggsignif]{geom_signif}}}
}
\description{
Plotting the distance between the samples with heatmap or boxplot.
}
\examples{
\dontrun{
data(mouse.time.mpse)
mouse.time.mpse \%<>\% mp_decostand(.abundance=Abundance)
mouse.time.mpse
mouse.time.mpse \%<>\% 
  mp_cal_dist(.abundance=hellinger, distmethod="bray")
mouse.time.mpse
p1 <- mouse.time.mpse \%>\% 
        mp_plot_dist(.distmethod=bray)
p2 <- mouse.time.mpse \%>\% 
        mp_plot_dist(.distmethod=bray, .group=time, group.test=TRUE)
p3 <- mouse.time.mpse \%>\% 
        mp_plot_dist(.distmethod=bray, .group=time)
}
}
\seealso{
[mp_cal_dist()] and [mp_extract_dist()]
}
\author{
Shuangbin Xu
}
