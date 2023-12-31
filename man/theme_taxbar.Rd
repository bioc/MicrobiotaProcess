% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-ggbartaxa.R
\name{theme_taxbar}
\alias{theme_taxbar}
\title{theme_taxbar}
\usage{
theme_taxbar(
  axis.text.x = element_text(angle = -45, hjust = 0, size = 8),
  legend.position = "bottom",
  legend.box = "horizontal",
  legend.text = element_text(size = 8),
  legend.title = element_blank(),
  strip.text.x = element_text(size = 12, face = "bold"),
  strip.background = element_rect(colour = "white", fill = "grey"),
  ...
)
}
\arguments{
\item{axis.text.x}{element_text, x axis tick labels.}

\item{legend.position}{character, default is "bottom".}

\item{legend.box}{character, arrangement of legends, default is "horizontal".}

\item{legend.text}{element_text, legend labels text.}

\item{legend.title}{element_text, legend title text}

\item{strip.text.x}{element_text, strip text of x}

\item{strip.background}{element_rect, the background of x}

\item{...}{additional parameters}
}
\value{
updated ggplot object with new theme
}
\description{
theme_taxbar
}
\examples{
\dontrun{
    library(ggplot2)
    data(test_otu_data)
    test_otu_data \%<>\% as.phyloseq()
    otubar <- ggbartax(test_otu_data, settheme=FALSE) + 
        xlab(NULL) + ylab("relative abundance(\%)") + 
        theme_taxbar()
}
}
\seealso{
\code{\link[ggplot2]{theme}}
}
