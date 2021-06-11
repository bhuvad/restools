#' @import ggplot2
NULL

#' Compute and plot the results of a PCA analysis on gene expression data
#'
#' @param edata a DGEList, SummarizedExperiment or ExpressionSet object
#'   containing gene expression data.
#' @param dims a numeric, containing 2 values specifying the dimensions to plot.
#' @param precomputed a dimensional reduction results from `stats::prcomp`.
#' @param rl a numeric, specifying the relative scale factor to apply to text on
#'   the plot.
#' @param ... aesthetic mappings to pass to `ggplot2::aes_string()`.
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' se = emtdata::cursons2018_se()
#' dge = emtdata::asDGEList(se)
#' plotPCA(dge, colour = Subline)
#'
setGeneric("plotPCA",
           function(edata,
                    dims = c(1, 2),
                    precomputed = NULL,
                    rl = 1,
                    ...) standardGeneric("plotPCA"))

#' @rdname plotPCA
setMethod("plotPCA",
          signature('DGEList','ANY', 'ANY', 'ANY'),
          function(edata, dims, precomputed, rl, ...){
            #compute PCA
            if (is.null(precomputed)) {
              pcdata = stats::prcomp(t(edgeR::cpm(edata, log = TRUE)))
            } else {
              pcdata = precomputed
            }

            #extract sample data
            sdata = edata$samples
            #create data structure
            drdf = pdataPC_intl(pcdata, dims)
            p1 = plotDR_intl(drdf, sdata, rl, ...)

            return(p1)
          })

#' @rdname plotPCA
setMethod("plotPCA",
          signature('ExpressionSet','ANY', 'ANY', 'ANY'),
          function(edata, dims, precomputed, rl, ...){
            #compute PCA
            if (is.null(precomputed)) {
              pcdata = stats::prcomp(t(Biobase::exprs(edata)))
            } else {
              pcdata = precomputed
            }

            #extract sample data
            sdata = Biobase::pData(edata)
            #create data structure
            drdf = pdataPC_intl(pcdata, dims)
            p1 = plotDR_intl(drdf, sdata, rl, ...)

            return(p1)
          })

#' @rdname plotPCA
setMethod("plotPCA",
          signature('SummarizedExperiment','ANY', 'ANY', 'ANY'),
          function(edata, dims, precomputed, rl, ...){
            #compute PCA
            if (is.null(precomputed)) {
              pcdata = stats::prcomp(t(SummarizedExperiment::assay(edata)))
            } else {
              pcdata = precomputed
            }

            #extract sample data
            sdata = BiocGenerics::as.data.frame(SummarizedExperiment::colData(edata))
            #create data structure
            drdf = pdataPC_intl(pcdata, dims)
            p1 = plotDR_intl(drdf, sdata, rl, ...)

            return(p1)
          })

#' Compute and plot the results of a MDS analysis on gene expression data
#'
#' @param precomputed a dimensional reduction results from either
#'   `limma::plotMDS`.
#'
#' @inheritParams plotPCA
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' se = emtdata::cursons2018_se()
#' dge = emtdata::asDGEList(se)
#' plotMDS(dge, colour = Subline)
#'
setGeneric("plotMDS",
           function(edata,
                    dims = c(1, 2),
                    precomputed = NULL,
                    rl = 1,
                    ...) standardGeneric("plotMDS"))

#' @rdname plotMDS
setMethod("plotMDS",
          signature('DGEList','ANY', 'ANY', 'ANY'),
          function(edata, dims, precomputed, rl, ...){
            #compute PCA
            if (is.null(precomputed)) {
              mdsdata = limma::plotMDS(edata, plot = FALSE)
            } else {
              mdsdata = precomputed
            }

            #extract sample data
            sdata = edata$samples
            #create data structure
            drdf = pdataMDS_intl(mdsdata, dims)
            p1 = plotDR_intl(drdf, sdata, rl, ...)

            return(p1)
          })

#' @rdname plotMDS
setMethod("plotMDS",
          signature('ExpressionSet','ANY', 'ANY', 'ANY'),
          function(edata, dims, precomputed, rl, ...){
            #compute PCA
            if (is.null(precomputed)) {
              mdsdata = limma::plotMDS(edata, plot = FALSE)
            } else {
              mdsdata = precomputed
            }

            #extract sample data
            sdata = Biobase::pData(edata)
            #create data structure
            drdf = pdataMDS_intl(mdsdata, dims)
            p1 = plotDR_intl(drdf, sdata, rl, ...)

            return(p1)
          })

#' @rdname plotMDS
setMethod("plotMDS",
          signature('SummarizedExperiment','ANY', 'ANY', 'ANY'),
          function(edata, dims, precomputed, rl, ...){
            #compute PCA
            if (is.null(precomputed)) {
              mdsdata = limma::plotMDS(SummarizedExperiment::assay(edata), plot = FALSE)
            } else {
              mdsdata = precomputed
            }

            #extract sample data
            sdata = BiocGenerics::as.data.frame(SummarizedExperiment::colData(edata))
            #create data structure
            drdf = pdataMDS_intl(mdsdata, dims)
            p1 = plotDR_intl(drdf, sdata, rl, ...)

            return(p1)
          })

#plot data preparation using PCA results
pdataPC_intl <- function(pcdata, dims) {
  stopifnot(length(dims) == 2)
  stopifnot(is(pcdata, 'prcomp'))

  plotdf = data.frame(
    'RestoolsMtchID' = rownames(pcdata$x),
    x = pcdata$x[, dims[1]],
    y = pcdata$x[, dims[2]]
  )
  pca_prop = round(pcdata$sdev ^ 2 / sum(pcdata$sdev ^ 2) * 100, digits = 2)
  pca_labs = paste0('PC', 1:length(pca_prop), ' (', pca_prop, '%)')
  colnames(plotdf)[-1] = pca_labs[dims]

  return(plotdf)
}

#plot data preparation using MDS results
pdataMDS_intl <- function(mdsdata, dims) {
  stopifnot(length(dims) == 2)
  stopifnot(is(mdsdata, 'MDS'))

  lambda = pmax(mdsdata$eigen.values, 0)
  plotdf = data.frame(
    'RestoolsMtchID' = rownames(mdsdata$distance.matrix.squared),
    x = mdsdata$eigen.vectors[, dims[1]] * lambda[dims[1]],
    y = mdsdata$eigen.vectors[, dims[2]] * lambda[dims[2]]
  )
  pca_labs = paste0('Leading logFC dim ', 1:ncol(mdsdata$eigen.vectors))
  pca_labs = paste0(pca_labs, ' (', round(mdsdata$var.explained * 100, digits = 2) , '%)')
  colnames(plotdf)[-1] = pca_labs[dims]

  return(plotdf)
}

#add colour annotation
addSampleAnnot <- function(plotdf, sdata) {
  plotdf = cbind(plotdf, sdata[plotdf$RestoolsMtchID, ])
  return(plotdf)
}

plotDR_intl <- function(drdf, sdata, rl, ...) {

  #annotate samples
  plotdf = addSampleAnnot(drdf, sdata)

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


  # aes requires x & y to be explicit: https://github.com/tidyverse/ggplot2/issues/3176
  x = rlang::sym(colnames(plotdf)[[2]])
  y = rlang::sym(colnames(plotdf)[[3]])

  #add size if not present
  if (is.null(defaultmap$size))
    defaultmap$size = 2

  # tidystyle recommends no explicit return statements at end of functions
  ggplot2::ggplot(plotdf, ggplot2::aes(!!x, !!y, !!!aesmap)) +
    ggplot2::geom_point() +
    ggplot2::update_geom_defaults('point', defaultmap) +
    vissE::bhuvad_theme(rl)
}
