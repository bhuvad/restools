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
#' plotPCA(dge, colour = "Subline")
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
#' plotMDS(dge, colour = "Subline")
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
  #extract aes
  aesmap = list(...)

  #annotate samples
  plotdf = addSampleAnnot(drdf, sdata)

  #compute plot
  aesmap = aesmap[setdiff(names(aesmap), c('x', 'y'))] #remove x,y mappings if present
  aesmap$x = colnames(plotdf)[2] #has to be positionally retained
  aesmap$y = colnames(plotdf)[3] #has to be positionally retained
  aesmap = lapply(aesmap, function(x) paste0("`", x, "`"))
  aesmap = do.call(ggplot2::aes_string, aesmap)

  #split aes params into those that are not aes i.e. static parametrisation
  is_aes = unlist(sapply(aesmap, function(x) x[[2]])) %in% colnames(plotdf)
  defaultmap = sapply(aesmap[!is_aes], function(x) x[[2]])
  aesmap = aesmap[is_aes]

  #add size if not present
  if (is.null(defaultmap$size))
    defaultmap$size = 2

  p1 = ggplot2::ggplot(plotdf, aesmap) +
    ggplot2::geom_point() +
    ggplot2::update_geom_defaults('point', defaultmap) +
    vissE::bhuvad_theme(rl)

  return(p1)
}
