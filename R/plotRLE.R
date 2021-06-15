#' Compute and plot relative log expression (RLE) values of gene expression data
#'
#' @param ordannots a character, specifying the annotations to use to order
#'   samples. Samples are ordered successively using each annotation. If the
#'   columns to use are "A" and "B" (from the annotation table), samples are
#'   first ordered using "A" and then "B" within each category of "A".
#'
#' @inheritParams plotPCA
#' @return a ggplot2 object, containing the RLE plot.
#' @export
#'
#' @examples
#' se = emtdata::cursons2018_se()
#' dge = emtdata::asDGEList(se)
#' plotRLE(dge, colour = Subline, lty = Treatment, lwd = 2, ordannots =
#' c(Subline, Treatment))
#'
setGeneric("plotRLE",
           function(edata,
                    ordannots = c(),
                    rl = 1,
                    ...) standardGeneric("plotRLE"))

#' @rdname plotRLE
setMethod("plotRLE",
          signature('DGEList','ANY', 'ANY'),
          function(edata, ordannots, rl, ...){
            #extract sample data
            sdata = edata$samples
            #extract expression data (and transform)
            edata = edgeR::cpm(edata, log = TRUE)
            #create data structure
            samporder = orderSamples(sdata, rlang::enquo(ordannots))
            rledf = pdataRLE_intl(edata, samporder)
            p1 = plotRLE_intl(rledf, sdata, rl, ...)

            return(p1)
          })

#' @rdname plotRLE
setMethod("plotRLE",
          signature('ExpressionSet','ANY', 'ANY'),
          function(edata, ordannots, rl, ...){
            #extract sample data
            sdata = Biobase::pData(edata)
            #extract expression data (and transform)
            edata = Biobase::exprs(edata)
            #create data structure
            samporder = orderSamples(sdata, rlang::enquo(ordannots))
            rledf = pdataRLE_intl(edata, samporder)
            p1 = plotRLE_intl(rledf, sdata, rl, ...)

            return(p1)
          })

#' @rdname plotRLE
setMethod("plotRLE",
          signature('SummarizedExperiment','ANY', 'ANY'),
          function(edata, ordannots, rl, ...){
            #extract sample data
            sdata = SummarizedExperiment::colData(edata)
            #extract expression data (and transform)
            edata = SummarizedExperiment::assay(edata)
            #create data structure
            samporder = orderSamples(sdata, rlang::enquo(ordannots))
            rledf = pdataRLE_intl(edata, samporder)
            p1 = plotRLE_intl(rledf, sdata, rl, ...)

            return(p1)
          })

#plot data preparation using MDS results
pdataRLE_intl <- function(emat, sampord) {
  #compute RLE
  rle = emat - Biobase::rowMedians(emat)
  #order samples
  rle = rle[, sampord]

  #compute boxplot
  rledf = t(apply(rle, 2, function(x) grDevices::boxplot.stats(x)$stats))
  rledf = as.data.frame(rledf)
  colnames(rledf) = c('ymin', 'lower', 'middle', 'upper', 'ymax')
  rledf$x = 1:nrow(rledf)
  rledf$RestoolsMtchID = rownames(rledf)

  return(rledf)
}

orderSamples <- function(sdata, ordannots) {
  # tidyeval for lists is tricky. See here: https://stackoverflow.com/questions/48847456/
  # TODO: testing it works as expected
  quo_list <- as.list(rlang::quo_squash(ordannots))[-1]
  #if no ordering provided, use default order
  if (length(quo_list) == 0){
    ord = TRUE
  } else {
    ord = rlang::eval_tidy(rlang::expr(order(!!!quo_list)), data = sdata)
  }

  return(rownames(sdata)[ord])
}

plotRLE_intl <- function(plotdf, sdata, rl, ...) {
  #constant - sample size at which standard plot becomes dense
  dense_thresh = 50

  #extract aes
  aesmap = rlang::enquos(...)

  #annotate samples
  plotdf = addSampleAnnot(plotdf, sdata)

  #compute plot
  aesmap = aesmap[!names(aesmap) %in% c('x', 'ymin', 'ymin', 'upper', 'middle', 'lower')] #remove fixed mappings if present

  #split aes params into those that are not aes i.e. static parametrisation
  if (length(aesmap) > 0) {
    is_aes = sapply(aesmap, rlang::quo_is_symbolic)
    defaultmap = lapply(aesmap[!is_aes], rlang::eval_tidy)
    aesmap = aesmap[is_aes]
  } else {
    defaultmap = list()
  }

  #build plot
  p1 = ggplot2::ggplot(plotdf, aes(x = x, y=middle, group = x, !!!aesmap)) +
    ggplot2::geom_boxplot(
      aes(ymin = ymin, ymax = ymax, upper = upper, middle = middle, lower = lower),
      stat = 'identity') +
    ggplot2::geom_hline(yintercept = 0, colour = 2, lty = 2) +
    ggplot2::ylab('Relative log expression') +
    ggplot2::update_geom_defaults('boxplot', defaultmap) +
    vissE::bhuvad_theme(rl) +
    ggplot2::theme(axis.text.x = element_blank())

  #update plot if too many samples are plot
  if (nrow(plotdf) > dense_thresh) {
    ## geom_point will inherit relevant aesthetics from top `aes`, include y=middle
    p1 = p1 +  ggplot2::geom_point()
  }

  return(p1)
}

plotRLEtm <- function(dge, clrannot, ordannots = NA, rl = 1) {
  stopifnot(clrannot %in% colnames(dge$samples))
  if (is.na(ordannots)) ordannots = clrannot
  stopifnot(ordannots %in% colnames(dge$samples))

  #compute RLE
  rle = edgeR::cpm(dge, log = TRUE)
  rle = rle - rowMedians(rle)
  #order samples
  rle = rle[, order(apply(dge$samples[, ordannots], 1, paste, collapse = '_'))]

  #compute boxplot
  rledf = t(apply(rle, 2, function(x) boxplot.stats(x)$stats))
  rledf = as.data.frame(rledf)
  colnames(rledf) = c('ymin', 'lower', 'middle', 'upper', 'ymax')
  rledf = cbind(rledf, dge$samples[rownames(rledf), clrannot, drop = FALSE])
  rledf$x = 1:nrow(rledf)
  p1 = ggplot(rledf, aes(x = x, group = x)) +
    geom_boxplot(
      aes_string(
        ymin = 'ymin',
        lower = 'lower',
        middle = 'middle',
        upper = 'upper',
        ymax = 'ymax',
        colour = clrannot,
      ),
      alpha = 0.5,
      lwd = 0.25,
      stat = 'identity'
    ) +
    geom_point(aes_string(y = 'middle', colour = clrannot)) +
    geom_hline(yintercept = 0, colour = 2, lty = 2) +
    xlab('Tissue') +
    ylab('Relative log expression') +
    # vissE::bhuvad_theme(rl) +
    theme(axis.text.x = element_blank())

  return(p1)
}
