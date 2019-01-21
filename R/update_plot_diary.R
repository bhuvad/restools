#' Update the plot_diary file for the project

update_plot_diary <- function(dsproj = '.') {
  #read README file
  rmdfile = file.path(dsproj, 'figures/plot_diary.Rmd')
  pltdiary = readLines(rmdfile)

  #update main sections
  readme = updateReadmeSection(dsproj, readme, 'code', 'raw_code')
  readme = updateReadmeSection(dsproj, readme, 'code', 'final_code')

  writeLines(readme, con = rmdfile)
  rmarkdown::render(rmdfile, quiet = TRUE)
}
