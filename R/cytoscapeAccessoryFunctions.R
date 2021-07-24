#' Get layout for a cytoscape network
#'
#' This function retrieves the layout of an active cytoscape network using the
#' `RCy3` package. The layout is retrieved as node coordinates.
#'
#' @param normalise a logical, indicating whether coordinates should be
#'   normalised using the z-transform or not.
#'
#' @return a data.frame with the x and y coordinates of nodes. Node names are
#'   stored as rownames of the data.frame.
#' @export
#'
#' @examples
#' # getCytoscapeLayout()
#'
getCytoscapeLayout <- function(normalise = TRUE) {
  #check whether cythoscape is running or not
  RCy3::cytoscapePing()

  #check for existing networks
  stopifnot(length(RCy3::getNetworkList()) > 0)

  #retrieve layout
  message('Retrieving layout for: ', RCy3::getNetworkName())

  #get layout
  lyt = RCy3::getNodePosition()

  #transform to numerics
  lyt$x_location = as.numeric(lyt$x_location)
  lyt$y_location = as.numeric(lyt$y_location)

  #switch y (cytoscape y-axis is inverted)
  lyt$y_location = -lyt$y_location

  #name axes
  colnames(lyt) = c('x', 'y')

  #normalise
  if (normalise){
    lyt$x = (lyt$x - mean(lyt$x)) / sd(lyt$x)
    lyt$y = (lyt$y - mean(lyt$y)) / sd(lyt$y)
  }

  return(lyt)
}
