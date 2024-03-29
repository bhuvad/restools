% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exploreVariability.R
\name{explorePCVariability}
\alias{explorePCVariability}
\alias{explorePCVariability,DGEList-method}
\alias{explorePCVariability,ExpressionSet-method}
\alias{explorePCVariability,SummarizedExperiment-method}
\title{Compute and plot the results of a PCA analysis on gene expression data}
\usage{
explorePCVariability(edata, annots, ...)

\S4method{explorePCVariability}{DGEList}(edata, annots, maxdim = 25, precomputed = NULL, ...)

\S4method{explorePCVariability}{ExpressionSet}(edata, annots, maxdim = 25, precomputed = NULL, ...)

\S4method{explorePCVariability}{SummarizedExperiment}(
  edata,
  annots,
  assay = 1,
  maxdim = 25,
  precomputed = NULL,
  ...
)
}
\arguments{
\item{edata}{a DGEList, SummarizedExperiment or ExpressionSet object
containing gene expression data.}

\item{annots}{variables to assess for contributions in each PC. ANOVA is used
for categorical variables and linear models for continuous variables.}

\item{...}{aesthetic mappings to pass to `ggplot2::aes_string()`.}

\item{maxdim}{a numeric, specifying the maximum number of dimensions to plot.}

\item{precomputed}{a dimensional reduction results from `stats::prcomp`.}

\item{assay}{a numeric or character, specifying the assay to use (for
SummarizedExperiment).}
}
\value{
a `Heatmap-class` object.
}
\description{
Compute and plot the results of a PCA analysis on gene expression data
}
\examples{
se = emtdata::cursons2015_se()
dge = emtdata::asDGEList(se)
explorePCVariability(dge, c(Subline, Treatment, lib.size))

}
