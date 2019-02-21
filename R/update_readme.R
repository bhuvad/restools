#' Update the README file for the project
#'
#' Scans the project directory provided and creates entries for files in the
#' README.Rmd file. Descriptions of files can then be modified manually.
#'
#' @param dsproj a character, providing the path of the project where the file
#' needs to be updated.
#'
#' @return NULL
#'
#' @details Only the description of files/directories should be modified in the
#' README.Rmd file. Modifying the comment tags or adding/removing text from the
#' file names will result in their removal in the next update.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # for "ds_project" in home folder
#' update_readme('~/ds_project')
#' }
update_readme <- function(dsproj = '.') {
  #read README file
  rmdfile = file.path(dsproj, 'README.Rmd')
  readme = readLines(rmdfile)

  #update main sections
  readme = updateReadmeSection(dsproj, readme, 'data')
  readme = updateReadmeSection(dsproj, readme, 'code', 'raw_code')
  readme = updateReadmeSection(dsproj, readme, 'code', 'final_code')
  readme = updateReadmeSection(dsproj, readme, 'figures', 'exploratory')
  readme = updateReadmeSection(dsproj, readme, 'figures', 'explanatory')
  readme = updateReadmeSection(dsproj, readme, 'products')

  writeLines(readme, con = rmdfile)
  rmarkdown::render(rmdfile, quiet = TRUE)
}

#locate a specific section of the readme file
locateReadmeSection <- function(dsproj = '.', readme, ...) {
  #create record templates
  dirpath = file.path(...)
  files = getFiles(dsproj, dirpath)
  frecords = paste0('* ** ', files, ' ** - description here')
  names(frecords) = files

  #create section tags for start and end and locate section using them
  btag = paste0('^\\h*<!--@', dirpath, ':begin-->') #begin
  etag = paste0('^\\h*<!--@', dirpath, ':end-->') #end
  bpos = which(grepl(btag, readme))
  epos = which(grepl(etag, readme))

  #error if either position is missing
  if (length(bpos) != 1 || length(epos) != 1) {
    stop("Tags may have been modified, make sure only descriptions have been changed")
  }

  #return null if no records
  if (epos - bpos == 1)
    return(NULL)

  return(list('begin' = bpos + 1, 'end' = epos - 1))
}

#extract a specific section of the readme file
extractReadmeSection <- function(readme, loc) {
  readme_section = readme[seq(loc$begin, loc$end)]
  return(readme_section)
}

#extract records from the readme section and create data.frame
getReadmeRecords <- function(content) {
  #create regular expressions for extraction
  regex_file = '`.+`' #file info in a single record
  regex_bullet = '^\\h*\\*\\h+'
  regex_dash = '\\h+-\\h+'
  regex_record = paste0(regex_bullet, regex_file, regex_dash)

  #locate records
  #add empty records to ease cutting
  record_loc = which(grepl(regex_record, content, perl = TRUE))
  if (length(record_loc) == 0)
    return(NULL)

  #split records
  cut_breaks = c(record_loc - 1, length(content))
  cut_factors = cut(seq(1, length(content)), breaks = cut_breaks)
  readme_records = split(content, cut_factors)
  readme_records = plyr::ldply(readme_records, function(rec) {
    #flatten record - no newlines
    rec = paste(rec, collapse = '\n')
    rec = stringr::str_replace_all(rec, '\\h*\n+\\h*', ' ')

    #extract file name
    fname = stringr::str_extract_all(rec, regex_file)[[1]]
    fname = stringr::str_replace_all(fname, '`', '')
    fname = stringr::str_trim(fname, 'both')

    #extract description
    desc = stringr::str_replace_all(rec, regex_record, '')
    desc = stringr::str_trim(desc, 'both')

    return(c('file' = fname, 'description' = desc))
  })[, -1]

  return(readme_records)
}

#convert each file to a markdown record
updateReadmeSection2 <- function(dsproj = '.', readme,  ...) {
  #locate section
  loc = locateReadmeSection(dsproj = '.', readme, ...)
  #return unmodified readme file  if no records exist for section
  if (is.null(loc))
    return(readme)

  #extract section content
  section_content = extractReadmeSection(readme, loc)
  record_df = getReadmeRecords(section_content)
  #return unmodified readme file  if no records exist for section
  if (is.null(record_df))
    return(readme)
}

#convert each file to a markdown record
updateReadmeSection <- function(dsproj = '.', readme,  ...) {
  #create record templates
  dirpath = file.path(...)
  files = getFiles(dsproj, dirpath)
  frecords = paste0('* ** ', files, ' ** - description here')
  names(frecords) = files

  #create section tags for start and end and locate section using them
  btag = paste0('<!--@', dirpath, ':begin-->') #begin
  etag = paste0('<!--@', dirpath, ':end-->') #end
  bpos = which(grepl(btag, readme))
  epos = which(grepl(etag, readme))

  #error if either position is missing
  if (length(bpos) != 1 || length(epos) != 1) {
    stop("Tags may have been modified, make sure only descriptions have been changed")
  }

  #extract relevant section and locate records
  dirrecord = readme[seq(bpos, epos)]
  dirrecord = readme[bpos + which(grepl('^*.*\\*\\*.*\\*\\*.*-', dirrecord)) - 1]
  dirrecord = plyr::ldply(stringr::str_split(dirrecord, ' - '))
  frecords =  plyr::ldply(stringr::str_split(frecords, ' - '))

  #if no files, remove the record from the README
  if (grepl('^*.*\\*\\*\\*\\*', frecords[1, 1])) {
    #replace original README section
    readme = c(
      readme[1:bpos],
      "",
      readme[epos:length(readme)]
    )

    return(readme)
  }

  #select information from old README where available
  if (nrow(dirrecord) > 0){
    colnames(frecords) = colnames(dirrecord) = c('file', 'desc')
    rownames(dirrecord) = dirrecord$file
    rownames(frecords) = frecords$file
    commonrecs = intersect(dirrecord$file, frecords$file)
    frecords[commonrecs, 'desc'] = dirrecord[commonrecs, 'desc']
  }
  frecords = apply(frecords, 1, paste, collapse = ' - ')
  names(frecords) = NULL

  #replace original README section
  readme = c(
    readme[1:bpos],
    "",
    frecords,
    "",
    readme[epos:length(readme)]
  )

  return(readme)
}

#get files within the folder of interest
getFiles <- function(dsproj = '.',  ...) {
  checkDSproject(dsproj)

  qrypath = file.path(dsproj, ...)
  files = list.files(qrypath, recursive = TRUE, include.dirs = TRUE, no.. = TRUE)

  #return if no files
  if (length(files) == 0)
    return(character())

  files = paste0('./', file.path(..., files))
  files = sort(files)

  return(files)
}

#get dirs within the folder of interest
getDirs <- function(dsproj = '.',  ...) {
  checkDSproject(dsproj)

  qrypath = file.path(dsproj, ...)
  files = list.dirs(qrypath, full.names = FALSE, recursive = TRUE)

  #return if no files
  if (length(files) == 0)
    return(character())

  files = paste0('./', file.path(..., files))
  files = sort(files)
  #remove trailing slashes
  pat = paste0(.Platform$file.sep, '$')
  files = stringr::str_remove(files, pat)

  return(files)
}

#stop if not a DS project
checkDSproject <- function(dsproj = '.') {
  stopifnot(isDSproject(dsproj))
}

#check if its a valid DS project
isDSproject <- function(dsproj = '.') {
  return(dir.exists(file.path(dsproj, '.dsproject')))
}
