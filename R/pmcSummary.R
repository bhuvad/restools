getPubmedSummary <- function(ids, email) {
  stopifnot(all(grepl('[0-9]+', ids)))
  checkInstalled('xml2')

  #build & run query
  res = postNCBISummaryQuery('pubmed', ids, email)
  res = xml2::read_xml(res)

  #extract linksets
  docs = xml2::xml_children(res)
  resdf = lapply(docs, function(x) {
    id = xml2::xml_text(xml2::xml_find_all(x, 'Id'))
    x = xml2::xml_find_all(x, 'Item')

    #create data.frame
    df = data.frame(
      'PMID' = id,
      'Year' = gsub('([0-9][0-9][0-9][0-9]).*', '\\1', getItem(x, 'PubDate')),
      'Journal' = getItem(x, 'FullJournalName'),
      'Title' = getItem(x, 'Title'),
      'Authors' = getItem(x, 'AuthorList')
    )

    return(df)
  })

  #process results
  resdf = do.call(rbind, resdf)

  return(resdf)
}

getItem <- function(doc, field) {
  allfields = sapply(doc, xml2::xml_attr, 'Name')

  #return empty if field doesn't exist
  if (!field %in% allfields)
    return('')

  #retrieve record
  value = doc[allfields %in% field]

  #split lists
  if (xml2::xml_length(value) > 1) {
    value = xml2::xml_children(value)
  }

  #extract text
  value = paste(sapply(value, xml2::xml_text), collapse = ', ')

  return(value)
}

postNCBISummaryQuery <- function(db, ids, email) {
  checkInstalled('httr')

  #create query
  baseurl = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
  params = paste0(
    'db=', db,
    '&id=', paste(ids, collapse = ','),
    '&tool=restools',
    '&email=', email
  )

  #run POST request
  res = httr::POST(baseurl, body = params)

  return(res)
}
