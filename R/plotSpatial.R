#' Plot spatial transcriptomic annotations per spot
#'
#' @param spe a SpatialExperiment object.
#' @param what a character specifying what aspect should be plot, "annotation",
#'   "expression" or reduced dimension ("reduceddim").
#' @param ...
#' @param img a logical indicating whether the tissue image (if present) should
#' @param imgAlpha
#' @inheritParams plotPCA
#' @return a ggplot2 object
#' @export
#'
#' @examples
plotSpots <-
  function(spe,
           what = c('annotation', 'expression', 'reduceddim'),
           ...,
           assay = SummarizedExperiment::assayNames(spe),
           dimred = SingleCellExperiment::reducedDimNames(spe),
           img = TRUE,
           crop = FALSE,
           imgAlpha = 1,
           rl = 1,
           circles = FALSE) {
    what = match.arg(what)

    #----extract aes----
    aesmap = rlang::enquos(...)
    #compute plot
    aesmap = aesmap[!names(aesmap) %in% c('x0', 'y0', 'x', 'y', 'sf')] #remove x,y mappings if present
    #split aes params into those that are not aes i.e. static parametrisation
    if (length(aesmap) > 0) {
      is_aes = sapply(aesmap, rlang::quo_is_symbolic)
      defaultmap = lapply(aesmap[!is_aes], rlang::eval_tidy)
      aesmap = aesmap[is_aes]
    } else {
      defaultmap = list()
    }

    #----prepare spot data----
    plotdf = switch (what,
      annotation = extractAnnotation(spe),
      expression = extractExpression(spe, assay),
      reduceddim = extractReducedDim(spe, dimred)
    )

    #----image data----
    imgdf = NULL
    sf = 1
    if (img) {
      imgdf = extractImage(spe)
      sf = SpatialExperiment::scaleFactors(spe)[1]
    }

    #----add spatial coordinates and scales----
    sf = rep(sf, ncol(spe))
    spatdf = SpatialExperiment::spatialCoords(spe)
    colnames(spatdf) = c('x', 'y')
    plotdf = cbind(plotdf, spatdf)
    plotdf = as.data.frame(plotdf)
    plotdf$sf = sf

    #crop
    if (crop & img) {
      #compute lims
      xlim = range((plotdf$x * plotdf$sf))
      ylim = range((plotdf$y * plotdf$sf))

      #filter image data - remove out-of-bounds pixels
      imgdf = imgdf[imgdf$x >= xlim[1] &
                  imgdf$x <= xlim[2] &
                  imgdf$y >= ylim[1] &
                  imgdf$y <= ylim[2], ]
      imgdf = do.call(rbind, imgdf)
    }

    #----plot----
    #initialise plot
    p1 = ggplot2::ggplot()

    #image
    if (img) {
      requirePkg('ggnewscale')
      p1 = p1 +
        ggplot2::geom_raster(ggplot2::aes(x, y, fill = colour),
                             alpha = imgAlpha,
                             data = imgdf) +
        ggplot2::scale_fill_identity() +
        ggnewscale::new_scale_fill()
    }

    #plot spots/circles
    if (circles) {
      requirePkg('ggforce')
      #add radius if not present
      if (!'r' %in% names(aesmap)) {
        aesmap = c(aesmap, r = rlang::quo_set_env(quo(r), rlang::quo_get_env(aesmap[[1]])))
        plotdf$r = ifelse(is.null(defaultmap$r), 100, defaultmap$r)
        defaultmap = defaultmap[setdiff(names(defaultmap), 'r')]
      }

      #add scale factor
      scquo = paste0(rlang::quo_text(aesmap$r), ' * sf')
      aesmap$r = rlang::set_expr(aesmap$r, rlang::parse_expr(scquo))

      p1 = p1 +
        ggforce::geom_circle(
          ggplot2::aes(, , x0 = x * sf, y0 = y * sf,!!!aesmap),
          data = plotdf
        ) +
        ggplot2::update_geom_defaults(ggforce::GeomCircle, defaultmap)
    } else {
      p1 = p1 +
        ggplot2::geom_point(
          ggplot2::aes(x = x * sf, y = y * sf,!!!aesmap),
          data = plotdf
        ) +
        ggplot2::update_geom_defaults('point', defaultmap)
    }

    #----theme----
    p1 = p1 +
      bhuvad_theme(rl) +
      theme(axis.text = element_blank())

    return(p1)
  }

extractAnnotation <- function(spe) {
  SummarizedExperiment::colData(spe)
}

extractExpression <- function(spe, assay) {
  SummarizedExperiment::assay(spe, assay) |>
    as.matrix() |>
    t()
}

extractReducedDim <- function(spe, dimred) {
  SingleCellExperiment::reducedDim(spe, type = dimred)
}

extractImage <- function(spe) {
  #determine number of images to expect
  nimg = 1 #use the first only for SPEs

  #get img data
  imgdf = SpatialExperiment::imgData(spe)[seq(nimg), 'data']
  names(imgdf) = seq(nimg)
  imgdf = lapply(imgdf, function(x) {
    x = as.matrix(SpatialExperiment::imgRaster(x))
    nc = ncol(x)
    nr = nrow(x)
    data.frame(
      x = rep(seq(nc), each = nr),
      y = rep(seq(nr), nc),
      colour = as.character(x)
    )
  })

  return(imgdf)
}
