#' Plot spatial transcriptomic annotations per spot
#'
#' @param spe a SpatialExperiment or SpatialExperimentList object.
#'   be plot in the background.
#' @param what
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
    isSPEList = is(spe, 'ExperimentList')

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
    if (img) {
      imgdf = extractImage(spe)
    }

    #----add spatial coordinates and scales----
    sf = 1
    if (isSPEList) {
      sf = SpatialExperiment::scaleFactors(spe)[seq(ExperimentList::nexp(spe))]
      sf = rep(sf, times = ExperimentList::elapply(spe, ncol))
    } else {
      sf = SpatialExperiment::scaleFactors(spe)[1]
      sf = rep(sf, ncol(spe))
    }
    plotdf = cbind(plotdf, SpatialExperiment::spatialCoords(spe))
    plotdf = as.data.frame(plotdf)
    plotdf$sf = sf

    #----add SPEList facets if necessary----
    if (isSPEList) {
      #add experiment indices to plotdf
      plotdf$SPEListFacet = spe@experimentIndex
      if (!is.null(ExperimentList::experimentNames(spe))) {
        #add names to plotdf
        plotdf$SPEListFacet = ExperimentList::experimentNames(spe)[plotdf$SPEListFacet]
        if (img)
          #add names to imgdf
          imgdf$SPEListFacet = ExperimentList::experimentNames(spe)[imgdf$SPEListFacet]
      }
    } else {
      plotdf$SPEListFacet = 1
    }

    #crop
    if (crop & img) {
      #crop image data
      nspe = unique(plotdf$SPEListFacet)
      imgdf = lapply(nspe, function(i) {
        ftr = plotdf$SPEListFacet %in% i

        #compute lims
        xlim = range((plotdf$pxl_row_in_fullres * plotdf$sf)[ftr])
        ylim = range((plotdf$pxl_col_in_fullres * plotdf$sf)[ftr])

        #filter image data - remove out-of-bounds pixels
        x = imgdf[imgdf$SPEListFacet %in% i &
                    imgdf$x >= xlim[1] &
                    imgdf$x <= xlim[2] &
                    imgdf$y >= ylim[1] &
                    imgdf$y <= ylim[2], ]
        return(x)
      })
      imgdf = do.call(rbind, imgdf)
    }

    #----plot----
    #initialise plot
    p1 = ggplot2::ggplot()

    #image
    if (img) {
      p1 = p1 +
        ggplot2::geom_raster(ggplot2::aes(x, y, fill = colour),
                             alpha = imgAlpha,
                             data = imgdf) +
        ggplot2::scale_fill_identity() +
        ggnewscale::new_scale_fill()
    }

    #plot spots/circles
    if (circles) {
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
          ggplot2::aes(, , x0 = pxl_row_in_fullres * sf, y0 = pxl_col_in_fullres * sf,!!!aesmap),
          data = plotdf
        ) +
        ggplot2::update_geom_defaults(ggforce::GeomCircle, defaultmap)
    } else {
      p1 = p1 +
        ggplot2::geom_point(
          ggplot2::aes(x = pxl_row_in_fullres * sf, y = pxl_col_in_fullres * sf,!!!aesmap),
          data = plotdf
        ) +
        ggplot2::update_geom_defaults('point', defaultmap)
    }

    #SPEList?
    if (isSPEList) {
      p1 = p1 +
        ggplot2::facet_wrap(~ SPEListFacet, scales = 'free')
    }

    #----theme----
    p1 = p1 +
      vissE::bhuvad_theme(rl) +
      theme(axis.text = element_blank())

    return(p1)
  }

extractAnnotation <- function(spe) {
  if (is(spe, 'SpatialExperimentList')) {
    SummarizedExperiment::colData(spe, experimentData = TRUE)
  } else {
    SummarizedExperiment::colData(spe)
  }
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
  if (is(spe, 'ExperimentList')) {
    nimg = ExperimentList::nexp(spe)
  }

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

  #combine and add image index (to use for SPEList)
  imgdf = mapply(function(df, i) {
    df$SPEListFacet = i
    return(df)
  }, imgdf, seq(nimg), SIMPLIFY = FALSE)
  imgdf = do.call(rbind, imgdf)

  return(imgdf)
}
