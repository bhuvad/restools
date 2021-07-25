diffusePubmedNetwork <- function(pmgraph, ids, weights = c('equal', 'cites'), r = 0.5, correct.for.hubs = TRUE) {
  checkInstalled('diffusr')
  checkInstalled('igraph')
  weights = match.arg(weights)

  #set initial probabilities (equal for each query)
  p0 = rep(0, length(igraph::V(pmgraph)))
  names(p0) = igraph::V(pmgraph)$name
  if (weights == 'equal') {
    p0[ids] = 1 / length(ids)
  } else {
    cites = igraph::degree(pmgraph, mode = 'in')
    cites = cites[ids]
    p0[ids] = cites / sum(cites)
  }

  #perform RWR
  amat = igraph::as_adj(igraph::as.undirected(pmgraph), sparse = FALSE)
  pt = diffusr::random.walk(p0, amat, r = r, correct.for.hubs)
  pt = pt$p.inf[, 1]
  names(pt) = names(p0)

  return(pt)
}

createPubmedNetwork <- function(pmdf, email, filter.degree = 1) {
  checkInstalled('igraph')

  #get neighbourhood and create a graph
  pmgraph = igraph::graph_from_data_frame(pmdf)

  return(pmgraph)
}

annotatePubmedNetwork <- function(pmgraph, ids = TRUE, email) {
  checkInstalled('igraph')

  #get ids
  vids = igraph::V(pmgraph)[ids]

  #add node attributes
  nodedf = getPubmedSummary(vids$name, email)
  rownames(nodedf) = nodedf$PMID
  nodedf = nodedf[igraph::V(pmgraph)$name, ]

  #add attributes to graph
  for (i in colnames(nodedf)[-1]) {
    igraph::vertex_attr(pmgraph, i, vids) = nodedf[vids$name, i]
  }

  return(pmgraph)
}

getPubmedNeighbourhood <- function(ids, order = 1, email) {
  completedids = c()
  nbdf = data.frame(from = character(), to = character(), order = numeric())

  for(i in 1:order) {
    #avoid rerunning already executed queries
    ids = setdiff(ids, completedids)
    if (length(ids) == 0)
      break

    #get upstream nodes i.e. A cites B
    updf = getPubmedUpstream(ids, 1, email)
    #get downstream nodes i.e. A cited by B
    dndf = getPubmedDownstream(ids, 1, email)
    tmpdf = rbind(updf, dndf)
    #add order
    if (nrow(tmpdf) ==  0)
      tmpdf = cbind(tmpdf, data.frame('order' = numeric()))
    else
      tmpdf$order = i

    #prepare for next run
    nbdf = rbind(nbdf, tmpdf)
    completedids = union(completedids, ids)
    ids = union(nbdf$from, nbdf$to)
  }

  return(nbdf)
}

getPubmedUpstream <- function(ids, order = 1, email) {
  resdf = getPubmedRecursive(ids, 'pubmed_pubmed_refs', email, order)
  return(resdf)
}

getPubmedDownstream <- function(ids, order = 1, email) {
  resdf = getPubmedRecursive(ids, 'pubmed_pubmed_citedin', email, order)
  resdf = resdf[, c(2:1, 3)]
  colnames(resdf)[1:2] = c('from', 'to')
  return(resdf)
}

getPubmedRecursive <- function(ids, linkname, email, order = 1) {
  completedids = c()
  resdf = data.frame(from = character(), to = character(), order = numeric())

  for(i in 1:order) {
    #avoid rerunning already executed queries
    ids = setdiff(ids, completedids)
    if (length(ids) == 0)
      break

    #get upstream nodes i.e. A cited by / cites B
    tmpdf = getPubmedLinks(ids, 'pubmed', linkname, email)

    #add order
    if (nrow(tmpdf) ==  0)
      tmpdf = cbind(tmpdf, data.frame('order' = numeric()))
    else
      tmpdf$order = i

    #prepare for next run
    resdf = rbind(resdf, tmpdf)
    completedids = union(completedids, ids)
    ids = unique(resdf$to)
  }

  return(resdf)
}

getPubmedLinks <- function(ids, dbname, linkname, email) {
  stopifnot(all(grepl('[0-9]+', ids)))
  checkInstalled('xml2')

  #build & run query
  res = postNCBIRefQuery(dbname, linkname, ids, email)
  res = xml2::read_xml(res)

  #extract linksets
  linksets = xml2::xml_children(res)
  refdf = lapply(linksets, function(x) {
    #get ID
    id_from = xml2::xml_text(xml2::xml_child(x, 'IdList'))

    #get ref IDs
    linkset = xml2::xml_child(x, 'LinkSetDb')
    #filter for links only
    linkset = xml2::xml_find_all(linkset, 'Link')
    id_to = xml2::xml_text(xml2::xml_children(linkset))

    #check for links
    if (length(id_to) == 0)
      df = data.frame(from = character(), to = character())
    else
      df = data.frame(from = id_from, to = id_to)

    return(df)
  })

  #process results
  refdf = do.call(rbind, refdf)
  return(refdf)
}

postNCBIRefQuery <- function(dbfrom, linkname, ids, email) {
  checkInstalled('httr')

  #create query
  baseurl = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi'
  params = paste0(
    'dbfrom=', dbfrom,
    '&linkname=', linkname,
    '&', createIdString(ids),
    '&tool=restools',
    '&email=', email
  )

  #run POST request
  res = httr::POST(baseurl, body = params)

  return(res)
}

createIdString <- function(pmcids) {
  paste(paste0('id=', pmcids), collapse = '&')
}

checkInstalled <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop('This function requires the "', pkg, '" R package.')
}
