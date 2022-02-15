#' Plot spatial transcriptomic annotations per spot
#'
#' @param spe a SpatialExperiment or SpatialExperimentList object.
#' @inheritParams plotPCA
#' @return a ggplot2 object
#' @export
#'
#' @examples
plotSpotAnnotation <- function(spe, rl = 1, ...) {
  isSPEList = is(spe, 'SpatialExperimentList')
  #annotate samples
  if (!isSPEList) {
    plotdf = cbind(SummarizedExperiment::colData(spe), SpatialExperiment::spatialCoords(spe))
  } else {
    plotdf = cbind(ExperimentList::colWithExperimentData(spe), SpatialExperiment::spatialCoords(spe))
    #facet annotation for SPEList
    plotdf$SPEListFacet = spe@experimentIndex
    if (!is.null(ExperimentList::experimentNames(spe)))
      plotdf$SPEListFacet = ExperimentList::experimentNames(spe)[plotdf$SPEListFacet]
  }
  plotdf = as.data.frame(plotdf)
  plotSpots(plotdf, isSPEList, rl, ...)
}

#' Plot spatial transcriptomic gene expression per spot
#'
#' @param spe a SpatialExperiment or SpatialExperimentList object.
#' @inheritParams plotPCA
#' @return a ggplot2 object
#' @export
#'
#' @examples
plotSpotExpression <- function(spe, assay, rl = 1, ...) {
  isSPEList = is(spe, 'SpatialExperimentList')
  #annotate samples
  plotdf = cbind(t(SummarizedExperiment::assay(spe, i = assay)), SpatialExperiment::spatialCoords(spe))
  plotdf = as.data.frame(plotdf)
  #facet annotation for SPEList
  if (isSPEList) {
    plotdf$SPEListFacet = spe@experimentIndex
    if (!is.null(ExperimentList::experimentNames(spe)))
      plotdf$SPEListFacet = ExperimentList::experimentNames(spe)[plotdf$SPEListFacet]
  }
  plotSpots(plotdf, isSPEList, rl, ...)
}

#' Plot spatial transcriptomic reduced dimension per spot
#'
#' @param spe a SpatialExperiment or SpatialExperimentList object.
#' @param dimred a string or integer scalar indicating the reduced dimension
#'   result in reducedDims(object) to plot.
#' @inheritParams plotPCA
#' @return a ggplot2 object
#' @export
#'
#' @examples
plotSpotReducedDim <- function(spe, dimred = SingleCellExperiment::reducedDimNames(spe), rl = 1, ...) {
  isSPEList = is(spe, 'SpatialExperimentList')
  #annotate samples
  plotdf = cbind(SingleCellExperiment::reducedDim(spe, type = dimred), SpatialExperiment::spatialCoords(spe))
  plotdf = as.data.frame(plotdf)
  #facet annotation for SPEList
  if (isSPEList) {
    plotdf$SPEListFacet = spe@experimentIndex
    if (!is.null(ExperimentList::experimentNames(spe)))
      plotdf$SPEListFacet = ExperimentList::experimentNames(spe)[plotdf$SPEListFacet]
  }
  plotSpots(plotdf, isSPEList, rl, ...)
}

plotSpots <- function(plotdf, isSPEList, rl = 1, ...) {
  #extract aes
  aesmap = rlang::enquos(...)
  #compute plot
  aesmap = aesmap[!names(aesmap) %in% c('x', 'y')] #remove x,y mappings if present

  #split aes params into those that are not aes i.e. static parametrisation
  if (length(aesmap) > 0) {
    is_aes = sapply(aesmap, rlang::quo_is_symbolic)
    defaultmap = lapply(aesmap[!is_aes], rlang::eval_tidy)
    aesmap = aesmap[is_aes]
  } else {
    defaultmap = list()
  }

  #add size if not present
  if (is.null(defaultmap$size))
    defaultmap$size = 1

  # tidystyle recommends no explicit return statements at end of functions
  p1 = ggplot2::ggplot(plotdf, ggplot2::aes(pxl_col_in_fullres, pxl_row_in_fullres, !!!aesmap)) +
    ggplot2::geom_point(pch = 19) +
    ggplot2::update_geom_defaults('point', defaultmap) +
    vissE::bhuvad_theme(rl) +
    theme(
      axis.text = element_blank()
    )

  #SPEList?
  if (isSPEList) {
    p1 = p1 +
      ggplot2::facet_wrap(~ SPEListFacet)
  }

  return(p1)
}
