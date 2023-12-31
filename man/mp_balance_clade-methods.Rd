% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/method-binary_offspring.R
\name{mp_balance_clade}
\alias{mp_balance_clade}
\alias{mp_balance_clade,MPSE-method}
\alias{mp_balance_clade,MPSE}
\alias{mp_balance_clade,tbl_mpse-method}
\alias{mp_balance_clade,tbl_mpse}
\alias{mp_balance_clade,grouped_df_mpse-method}
\alias{mp_balance_clade,grouped_df_mpse}
\title{Calculating the balance score of internal nodes (clade) according to the geometric.mean/mean/median abundance of their binary children tips.}
\usage{
mp_balance_clade(
  .data,
  .abundance = NULL,
  force = FALSE,
  relative = TRUE,
  balance_fun = c("geometric.mean", "mean", "median"),
  pseudonum = 0.001,
  action = "get",
  ...
)

\S4method{mp_balance_clade}{MPSE}(
  .data,
  .abundance = NULL,
  force = FALSE,
  relative = TRUE,
  balance_fun = c("geometric.mean", "mean", "median"),
  pseudonum = 0.001,
  action = "get",
  ...
)

\S4method{mp_balance_clade}{tbl_mpse}(
  .data,
  .abundance = NULL,
  force = FALSE,
  relative = TRUE,
  balance_fun = c("geometric.mean", "mean", "median"),
  pseudonum = 0.001,
  action = "get",
  ...
)

\S4method{mp_balance_clade}{grouped_df_mpse}(
  .data,
  .abundance = NULL,
  force = FALSE,
  relative = TRUE,
  balance_fun = c("geometric.mean", "mean", "median"),
  pseudonum = 0.001,
  action = "get",
  ...
)
}
\arguments{
\item{.data}{MPSE object which must contain otutree slot, required}

\item{.abundance}{the column names of abundance.}

\item{force}{logical whether calculate the (relative) abundance forcibly when the abundance
is not be rarefied, default is FALSE.}

\item{relative}{logical whether calculate the relative abundance.}

\item{balance_fun}{function the method to calculate the (relative) abundance of internal nodes
according to their children tips, default is 'geometric.mean', other options are 'mean' and 'median'.}

\item{pseudonum}{numeric add a pseudo numeric to avoid the error of division in calculation, default 
is 0.001 .}

\item{action}{character, "add" joins the new information to the otutree slot if it exists (default).
In addition, "only" return a non-redundant tibble with the just new information. "get" return a new 'MPSE' 
object, and the 'OTU' column is the internal nodes and 'Abundance' column is the balance scores.}

\item{...}{additional parameters, meaningless now.}
}
\value{
a object according to 'action' argument.
}
\description{
Calculating the balance score of internal nodes (clade) according to the geometric.mean/mean/median abundance of their binary children tips.
}
\examples{
\dontrun{
  suppressPackageStartupMessages(library(curatedMetagenomicData))
  xx <- curatedMetagenomicData('ZellerG_2014.relative_abundance', dryrun=F)
  xx[[1]] \%>\% as.mpse -> mpse
  mpse.balance.clade <- mpse \%>\%
    mp_balance_clade(
      .abundance = Abundance,
      force = TRUE,
      relative = FALSE,
      action = 'get',
      pseudonum = .01
    )
  mpse.balance.clade 

  # Performing the Euclidean distance or PCA.

  mpse.balance.clade \%>\%
    mp_cal_dist(.abundance = Abundance, distmethod = 'euclidean') \%>\%
    mp_plot_dist(.distmethod = 'euclidean', .group = disease, group.test = T)

  mpse.balance.clade \%>\%
    mp_adonis(.abundance = Abundance, .formula=~disease, distmethod = 'euclidean', permutation = 9999)

  mpse.balance.clade \%>\%
    mp_cal_pca(.abundance = Abundance) \%>\% 
    mp_plot_ord(.group = disease)

  # Detecting the signal balance nodes.
  mpse.balance.clade \%>\% mp_diff_analysis(
      .abundance = Abundance,
      force = TRUE,
      relative = FALSE,
      .group = disease,
      fc.method = 'compare_mean'
  )
}
}
\references{
Morton JT, Sanders J, Quinn RA, McDonald D, Gonzalez A, Vázquez-Baeza Y, 
Navas-Molina JA, Song SJ, Metcalf JL, Hyde ER, Lladser M, Dorrestein PC, 
Knight R. 2017. Balance trees reveal microbial niche differentiation. 
mSystems 2:e00162-16. https://doi.org/10.1128/mSystems.00162-16.

Justin D Silverman, Alex D Washburne, Sayan Mukherjee, Lawrence A David. 
A phylogenetic transform enhances analysis of compositional microbiota data. 
eLife 2017;6:e21887. https://doi.org/10.7554/eLife.21887.001.
}
