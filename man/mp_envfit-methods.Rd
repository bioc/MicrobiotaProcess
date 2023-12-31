% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/method-mp_envfit.R
\name{mp_envfit}
\alias{mp_envfit}
\alias{mp_envfit,MPSE-method}
\alias{mp_envfit,MPSE}
\alias{mp_envfit,tbl_mpse-method}
\alias{mp_envfit,tbl_mpse}
\alias{mp_envfit,grouped_df_mpse-method}
\alias{mp_envfit,grouped_df_mpse}
\title{Fits an Environmental Vector or Factor onto an Ordination With MPSE or tbl_mpse Object}
\usage{
mp_envfit(
  .data,
  .ord,
  .env,
  .dim = 3,
  action = "only",
  permutations = 999,
  seed = 123,
  ...
)

\S4method{mp_envfit}{MPSE}(
  .data,
  .ord,
  .env,
  .dim = 3,
  action = "only",
  permutations = 999,
  seed = 123,
  ...
)

\S4method{mp_envfit}{tbl_mpse}(
  .data,
  .ord,
  .env,
  .dim = 3,
  action = "only",
  permutations = 999,
  seed = 123,
  ...
)

\S4method{mp_envfit}{grouped_df_mpse}(
  .data,
  .ord,
  .env,
  .dim = 3,
  action = "only",
  permutations = 999,
  seed = 123,
  ...
)
}
\arguments{
\item{.data}{MPSE or tbl_mpse object}

\item{.ord}{a name of ordination, option it is DCA, NMDS, RDA, CCA.}

\item{.env}{the names of columns of sample group or environment information.}

\item{.dim}{integer The number of dimensions to be returned, default is 3.}

\item{action}{character "add" joins the envfit result to internal attributes of the object, 
"only" return a non-redundant tibble with the envfit result. "get" return 'envfit' object can
be analyzed using the related vegan funtion.}

\item{permutations}{the number of permutations required, default is 999.}

\item{seed}{a random seed to make the analysis reproducible, default is 123.}

\item{...}{additional parameters see also 'vegan::envfit'}
}
\value{
update object according action
}
\description{
Fits an Environmental Vector or Factor onto an Ordination With MPSE or tbl_mpse Object
}
\examples{
library(vegan)
data(varespec, varechem)
mpse <- MPSE(assays=list(Abundance=t(varespec)), colData=varechem)
envformula <- paste("~", paste(colnames(varechem), collapse="+")) \%>\% as.formula
mpse \%<>\% 
       mp_cal_cca(.abundance=Abundance, .formula=envformula, action="add")
mpse2 <- mpse \%>\%
         mp_envfit(.ord=cca, 
                   .env=colnames(varechem), 
                   permutations=9999, 
                   action="add")
mpse2 \%>\% mp_plot_ord(.ord=cca, .group=Al, .size=Mn, show.shample=TRUE, show.envfit=TRUE)
\dontrun{
tbl <- mpse \%>\%
       mp_envfit(.ord=CCA, 
                 .env=colnames(varechem), 
                 permutations=9999, 
                 action="only")
tbl
library(ggplot2)
library(ggrepel)
x <- names(tbl)[grepl("^CCA1 ", names(tbl))] \%>\% as.symbol()
y <- names(tbl)[grepl("^CCA2 ", names(tbl))] \%>\% as.symbol()
p <- tbl \%>\%
     ggplot(aes(x=!!x, y=!!y)) + 
     geom_point(aes(color=Al, size=Mn)) + 
     geom_segment(data=dr_extract(
                            name="CCA_ENVFIT_tb", 
                            .f=td_filter(pvals<=0.05 & label!="Humdepth")
                       ), 
                  aes(x=0, y=0, xend=CCA1, yend=CCA2), 
                  arrow=arrow(length = unit(0.02, "npc"))
     ) + 
     geom_text_repel(data=dr_extract(
                              name="CCA_ENVFIT_tb", 
                              .f=td_filter(pvals<=0.05 & label!="Humdepth")
                          ), 
                  aes(x=CCA1, y=CCA2, label=label)
     ) +
     geom_vline(xintercept=0, color="grey20", linetype=2) +
     geom_hline(yintercept=0, color="grey20", linetype=2) +
     theme_bw() +
     theme(panel.grid=element_blank())
p
}
}
\author{
Shuangbin Xu
}
