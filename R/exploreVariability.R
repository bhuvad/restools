#' Compute and plot the results of a PCA analysis on gene expression data
#'
#' @param annots variables to assess for contributions in each PC. ANOVA is used
#'   for categorical variables and linear models for continuous variables.
#' @param maxdim a numeric, specifying the maximum number of dimensions to plot.
#'
#' @inheritParams plotPCA
#' @return a `Heatmap-class` object.
#' @export
#'
#' @examples
#' se = emtdata::cursons2015_se()
#' dge = emtdata::asDGEList(se)
#' explorePCVariability(dge, c(Subline, Treatment, lib.size))
#'
setGeneric("explorePCVariability",
           function(edata,
                    annots,
                    maxdim = 25,
                    precomputed = NULL) standardGeneric("explorePCVariability"))

#' @rdname explorePCVariability
setMethod("explorePCVariability",
          signature('DGEList','ANY', 'ANY', 'ANY'),
          function(edata, annots, maxdim, precomputed){
            #compute PCA
            if (is.null(precomputed)) {
              pcdata = stats::prcomp(t(edgeR::cpm(edata, log = TRUE)))
            } else {
              pcdata = checkPrecomputedPCA(edata, precomputed)
            }

            #extract sample data
            sdata = edata$samples
            sdata = extractAnnots(sdata, rlang::enquo(annots))
            #create data structure
            drmat = pmatPC_intl(pcdata, maxdim)
            exploreDRVariability_intl(drmat, sdata)
          })

#' @rdname explorePCVariability
setMethod("explorePCVariability",
          signature('ExpressionSet','ANY', 'ANY', 'ANY'),
          function(edata, annots, maxdim, precomputed){
            #compute PCA
            if (is.null(precomputed)) {
              pcdata = stats::prcomp(t(Biobase::exprs(edata)))
            } else {
              pcdata = checkPrecomputedPCA(edata, precomputed)
            }

            #extract sample data
            sdata = Biobase::pData(edata)
            sdata = extractAnnots(sdata, rlang::enquo(annots))
            #create data structure
            drmat = pmatPC_intl(pcdata, maxdim)
            exploreDRVariability_intl(drmat, sdata)
          })

#' @rdname explorePCVariability
setMethod("explorePCVariability",
          signature('SummarizedExperiment','ANY', 'ANY', 'ANY'),
          function(edata, annots, maxdim, precomputed){
            #compute PCA
            if (is.null(precomputed)) {
              pcdata = stats::prcomp(t(SummarizedExperiment::assay(edata)))
            } else {
              pcdata = checkPrecomputedPCA(edata, precomputed)
            }

            #extract sample data
            sdata = BiocGenerics::as.data.frame(SummarizedExperiment::colData(edata))
            sdata = extractAnnots(sdata, rlang::enquo(annots))
            #create data structure
            drmat = pmatPC_intl(pcdata, maxdim)
            exploreDRVariability_intl(drmat, sdata)
          })


extractAnnots <- function(sdata, annots) {
  sdata = dplyr::select(sdata, !!annots)

  #check that atleast one annotation is present
  if (ncol(sdata) == 0) {
    stop('No annotations present')
  }

  return(sdata)
}

pmatPC_intl <- function(pcdata, maxdim) {
  maxdim = min(maxdim, length(pcdata$sdev))

  #generate labels
  pca_prop = round(pcdata$sdev ^ 2 / sum(pcdata$sdev ^ 2) * 100, digits = 2)
  pca_labs = paste0('PC', 1:length(pca_prop), ' (', pca_prop, '%)')
  #extract projection data
  pcmat = pcdata$x
  colnames(pcmat) = pca_labs
  pcmat = pcmat[, 1:maxdim]

  return(pcmat)
}

exploreDRVariability_intl <- function(drmat, sdata) {
  #convert to factors where needed
  sdata = as.data.frame(lapply(sdata, function(x) {
    if(!is(x, 'numeric'))
      factor(x)
    else
      x
  }))

  #specify model
  fit = lm(drmat ~ ., data = sdata)
  plotmat = lapply(summary(stats::aov(fit)), function(x) x[['Pr(>F)']])
  plotmat = do.call(rbind, plotmat)
  plotmat = plotmat[, 1:ncol(sdata), drop = FALSE]
  colnames(plotmat) = colnames(sdata)
  rownames(plotmat) = colnames(drmat)

  #adjust p.values
  plotmat = apply(plotmat, 2, stats::p.adjust)

  #Scree plot as a row annotation
  #extract proportions from names
  pca_prop = as.numeric(gsub('.*\\((.*)%\\)', '\\1', rownames(plotmat)))
  ha = ComplexHeatmap::HeatmapAnnotation(
    'Scree\nplot' = ComplexHeatmap::anno_lines(pca_prop, add_points = TRUE),
    which = 'row',
    width = unit(2, 'cm')
  )

  #heatmap colours
  colsat = as.numeric(stats::quantile(-log10(plotmat[plotmat < 0.01]), 0.75))
  colsat = ifelse(is.na(colsat), 2, colsat)
  cols = c('#FFFFFF', scico::scico(30, palette = 'acton', direction = -1))
  hm_colmap = circlize::colorRamp2(seq(0, colsat, length.out = 31), cols)

  #heatmap
  ComplexHeatmap::Heatmap(
    -log10(plotmat),
    name = '-log10(p)',
    col = hm_colmap,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = 'left',
    right_annotation = ha,
    cell_fun = function(j, i, x, y, width, height, fill) {
      text = ''
      if (plotmat[i, j] < 0.001) {
        text = '***'
      } else if (plotmat[i, j] < 0.01) {
        text = '**'
      } else if (plotmat[i, j] < 0.1) {
        text = '*'
      }
      grid::grid.text(
        text,
        x,
        y,
        gp = grid::gpar(
          fontsize = 15,
          col = '#F7FCF0',
          fontface = 'bold'
        ),
        hjust = 0.5,
        just = 0
      )
    }
  )
}
