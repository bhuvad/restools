#' Create directory structure for a data science project
#'
#' Initialise a data science project. This function creates a directory
#' structure to help you organise your project. The structure is inspired from
#' JEFFERY LEEKS TUTORIAL. Along with the directory structure, it creates a
#' README file in Rmarkdown which can then be used to produce a GitHub
#' compatible markdown document. It also creates a plot diary to help you
#' organise plots generated during the project.
#'
#' @param pname a character, stating the name of the project
#' @param path a character, specifying the location of the project
#'
#' @return NULL
#'
#' @details The README.Rmd file is created in the project directory and the
#' plot_diary.Rmd created in the auto generated figures folder
#'
#' @export
#'
#' @seealso \code{\link{updateReadme}}
#'
#' @examples
#' \dontrun{
#' # create project named "ds_project" in home folder
#' createDSProject('ds_project', '~')
#' }
#'
createDSProject <- function(pname = NULL, path = '.') {
  stopifnot(!is.null(pname))

  #create project directory
  dpath = file.path(path, pname)
  created = dir.create(dpath)
  if (created) {
    message(paste('project directory created:', dpath))
  } else{
    warning('directory exists, proceeding with setup')
  }

  #create R project with defaults
  init_rproj(pname, path)

  #create subsequent directories
  dir.create(file.path(dpath, '.dsproject'))                    #.dsproject
  dir.create(file.path(dpath, 'literature'))                    #literature
  dir.create(file.path(dpath, 'data'))                          #data
  dir.create(file.path(dpath, 'code'))                          #code
  dir.create(file.path(dpath, 'code', 'raw_code'))              #code/raw_code
  dir.create(file.path(dpath, 'code', 'final_code'))            #code/final_code
  dir.create(file.path(dpath, 'figures'))                       #figures
  dir.create(file.path(dpath, 'figures', 'exploratory'))        #figures/exploratory
  dir.create(file.path(dpath, 'figures', 'explanatory'))        #figures/explanatory
  dir.create(file.path(dpath, 'products'))                      #products

  #create Rmds for README
  readmef = file.path(dpath, 'README.Rmd')
  rmarkdown::draft(readmef, template = 'readme', package = 'restools', edit = FALSE)
  rmarkdown::render(readmef, quiet = TRUE)

  #create Rmds for lit_review
  cat('', file = file.path(dpath, 'literature', 'bibliography.bib'))
  litreviewf = file.path(dpath, 'literature', 'lit_review.Rmd')
  rmarkdown::draft(litreviewf, template = 'lit_review', package = 'restools', edit = FALSE)
  rmarkdown::render(litreviewf, quiet = TRUE)
}

init_git <- function() {
  #initialise a git repo for the code folder
}

init_config <- function() {

}

init_rproj <- function(pname, path) {
  proj_text = 'Version: 1.0

  RestoreWorkspace: Default
  SaveWorkspace: Default
  AlwaysSaveHistory: Default

  EnableCodeIndexing: Yes
  UseSpacesForTab: Yes
  NumSpacesForTab: 2
  Encoding: UTF-8

  RnwWeave: Sweave
  LaTeX: pdfLaTeX'
  proj_text = strsplit(proj_text, '\n')[[1]]
  proj_text = gsub("^\\s+|\\s+$", "", proj_text)

  #write project string to project directory
  writeLines(proj_text, con = file.path(path, pname, paste0(pname, '.Rproj')))
}
