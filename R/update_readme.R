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
#' updateReadme('~/ds_project')
#' }
#'
updateReadme <- function(dsproj = '.') {
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

#extract a specific section of the readme file
replaceReadmeSection <- function(readme, loc, new_records) {
  new_readme = c(
    readme[seq(1, loc$begin - 1)],
    '',
    new_records,
    '',
    readme[seq(loc$end + 1, length(readme))]
  )

  return(new_readme)
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
    fname = stringr::str_replace_all(fname, '/$', '')

    #extract description
    desc = stringr::str_replace_all(rec, regex_record, '')
    desc = stringr::str_trim(desc, 'both')

    return(c('file' = fname, 'description' = desc))
  })[, -1]

  return(readme_records)
}

#create records for the README file from a records dataframe
createReadmeRecords <- function(records) {
  records = apply(records, 1, function(x) {
    paste0('* `', x['file'], '` - ', x['description'])
  })
  names(records) = NULL

  return(records)
}

#convert each file to a markdown record
updateReadmeSection <- function(dsproj = '.', readme,  ...) {
  #----locate section----
  loc = locateReadmeSection(dsproj = '.', readme, ...)
  #return unmodified readme file  if no records exist for section
  if (is.null(loc))
    return(readme)

  #----extract section content----
  section_content = extractReadmeSection(readme, loc)
  record_df = getReadmeRecords(section_content)

  #----get list of existing files----
  dirpath = file.path(...)
  files = getFiles(dsproj, dirpath)

  #----filter out files whose dir record has a truncation mark----
  #e.g. ./dir1/dir2/* results in no annotations for files/dirs in dir2
  #i.e. effectively stop recursion

  if (!is.null(record_df)) {
    #identify truncation records
    trunc_dirs = record_df$file[grepl('\\*$', record_df$file)]
    #identify files within truncation directories
    trunc_dirs_names = stringr::str_sub(trunc_dirs, start = 3, end = -3)
    rm_files = unlist(lapply(trunc_dirs_names, function(d) getFiles(dsproj, d)))
    #add truncated dirs to rm_files for removal
    rm_files = c(rm_files, paste0('./', trunc_dirs_names))
    #remove these files from list of files whose records will be created
    files = setdiff(files, rm_files)
    files = union(files, trunc_dirs)
  } else{
    record_df = data.frame('file' = character(), 'description' = character(), stringsAsFactors = FALSE)
  }

  #----compute new record df----
  updated_record_df = merge(data.frame('file' = files), record_df, all.x = TRUE, sort = TRUE)
  updated_record_df$description[is.na(updated_record_df$description)] = 'description here'

  #----add figure thumbnails for figure sections----
  if (grepl('^figures/', dirpath, perl = TRUE) & nrow(updated_record_df) != 0)
    updated_record_df = updateFigureDescriptions(updated_record_df)

  #----replace section----
  new_records = createReadmeRecords(updated_record_df)
  readme = replaceReadmeSection(readme, loc, new_records)
}

#for figures, update the description to include a link to the figure
updateFigureDescriptions <- function(records, max_size_mb = 10) {
  regex_img_link = '\\h*\\!\\[.*\\]\\(.*\\)\\h*'

  max_size_mb = max_size_mb * 4 #knitr by default reduces height and width by 50%
  allowed_ftypes = c('png', 'jpg', 'jpeg', 'tiff', 'svg', 'bmp')

  records$file = as.character(records$file)
  records$description = as.character(records$description)

  #allow pdf display for explanatory figures
  if (all(grepl('explanatory', records$file))) {
    allowed_ftypes = c(allowed_ftypes, 'pdf')
  }

  #remove old links
  records$description = sapply(records$description, function(x) {
    x = stringr::str_replace_all(x, regex_img_link, '')
    return(x)
  })

  #add new links/thumbnails
  figext = tools::file_ext(records$file)
  figsizes = file.size(records$file)
  thumbs = figext %in% allowed_ftypes

  #remove large files such that the total size of all files is max_size_mb
  thumbs = thumbs & figsizes < getSizeThresh(figsizes[thumbs], max_size_mb)

  #if all items are too large, don't add thumbnails
  if (all(!thumbs)) {
    return(records)
  }

  records$description[thumbs] = mapply(paste0,
                                       records$description[thumbs],
                                       ' ![fig_thumbnail](',
                                       records$file[thumbs],
                                       ')')
  return(records)
}

#get a threshold for image file sizes that results in a total of N MB
getSizeThresh <- function(sizes, nMB = 10) {
  sizes = sort(sizes)
  thresh = max(sizes[cumsum(sizes / 2^20) < nMB], 0)
  return(thresh + 1)
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
