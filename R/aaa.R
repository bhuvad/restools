requirePkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)){
    stop(sprintf("'%s' needed for this function to work. Please install it.", pkg), call. = FALSE)
  }
}

bhuvad_theme <- function (rl = 1.1) {
  stopifnot(rl > 0)
  ggplot2::theme_minimal() + ggplot2::theme(
    panel.border = element_rect(colour = "black",
                                fill = NA),
    panel.grid = element_blank(),
    axis.title = element_text(size = rel(rl) *
                                1.1),
    axis.text = element_text(size = rel(rl)),
    plot.title = element_text(size = rel(rl) *
                                1.2),
    strip.background = element_rect(fill = NA, colour = "black"),
    strip.text = element_text(size = rel(rl)),
    legend.text = element_text(size = rel(rl)),
    legend.title = element_text(size = rel(rl), face = "italic")
  )
}
