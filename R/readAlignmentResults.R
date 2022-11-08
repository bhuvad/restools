#' Read results of a Subread analysis into a DGEList object
#'
#' Creates a DGEList object using files from a Subread analysis. The function
#' will look for a file called 'SraRunTable.txt' for sample annotations, a
#' '*_counts.txt' file(s) for expression data and 'counts.txt.summary' for
#' summary statistics (optional).
#'
#' @param fpath a character, specifying the path of the analysis (ideally using
#'   bhuva.d's scripts).
#' @param plot a logical, stating whether to plot summary statistics for the
#'   alignment or not.
#'
#' @return a DGEList object, containing transcriptomics data
#' @export
#'
#' @examples
readFeatureCountsRes <- function(fpath, plot = TRUE) {
  requirePkg('edgeR')

  count_file = list.files(fpath, 'fc_counts.txt$', recursive = TRUE, full.names = TRUE)
  stats_file = list.files(fpath, 'fc_counts.txt.summary$', recursive = TRUE, full.names = TRUE)
  annot_file = list.files(fpath, 'SraRunTable.txt$', recursive = TRUE, full.names = TRUE)

  #read counts
  if (length(count_file) == 1) {
    counts = read.table(count_file, header = TRUE, check.names = FALSE, row.names = 1)

    #read files
    gannot =  counts[, 1:7]
    counts = counts[, -(1:7)]
    colnames(counts) = gsub('\\.bam$', '', basename(colnames(counts)))
  } else if (length(count_file) > 1) {
    stop('Need to implement this bit!')
  } else {
    stop('Counts file not found')
  }

  #plot summary stats
  if (length(stats_file) == 1 & plot) {
    colmap = c('#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99', '#E31A1C', '#FDBF6F', '#FF7F00', '#CAB2D6', '#6A3D9A', '#FFFF99', '#B15928', '#666666', '#000000')

    stats = read.table(stats_file, header = TRUE, check.names = FALSE)
    colnames(stats)[-1] = gsub('\\.bam$', '', basename(colnames(stats)[-1]))
    p1 = stats |>
      tidyr::pivot_longer(!dplyr::contains('status'), names_to = 'Sample', values_to = 'NumReads') |>
      ggplot2::ggplot(ggplot2::aes(Sample, NumReads / 1e6, fill = Status)) +
        ggplot2::geom_bar(stat = 'identity') +
        ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(values = colmap) +
        ggplot2::labs(x = '', y = 'Million reads') +
        bhuvad_theme()
    print(p1)
  }

  #read sample annotation
  if (length(annot_file) == 1) {
    sannot = read.csv(annot_file, row.names = 1)
    if (!all(rownames(sannot) %in% colnames(counts))) {
      warning('some sample annotations are missing, not using sample annotations')
      sannot = sannot[rownames(sannot) %in% colnames(counts), ]
    }
  } else {
    sannot = NULL
  }

  #create DGEList
  dge = edgeR::DGEList(counts, genes = gannot, samples = sannot)

  return(dge)
}
