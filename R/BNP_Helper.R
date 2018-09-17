####################
####   HELPER   ####
####################

#' Convert integer to binary state vector with node names.
#' '
#' @param x input integer representing the state
#' @param node.names network node names
#' 
#' @return Numeric binary vector of the same length as the number of nodes. Each 
#'     position corresponds to a node of the network. The values of each element 
#'     are 0 or 1. The name of each element corresponds to the name of the node in
#'     that position in the network.
#'     
#' @seealso \code{\link{intToBits}} which this function wraps
#' 
#' @examples
#' data(cellcycle)
#' int2binState(162,cellcycle$genes)
#' 
#' @keywords internal
#' @export
int2binState <- function(x, node.names){ 
    state <- as.integer( intToBits(x)[1:length(node.names)] )  
    names(state) <- node.names
    state
}

#' Convert binart state vector to integer.
#' '
#' @param x input vector representing the state
#' 
#' @return Corresponding int.
#' 
#' @seealso \code{\link{intToBits}} which this function wraps
#' 
#' @examples
#' bin2intState( c(0,1,0,0,0,1,0,1,0,0) )
#' 
#' @keywords internal
#' @export
bin2intState <- function(x){ 
  x <- rev(x)
  sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))
}

#' Validate and convert a number to a valid state.
#' 
#' @param x input integer representing the state, can be a vector of 0/1 or an int
#' @param node.names network node names
#' 
#' @return Numeric binary vector of the same length as the number of nodes. Each 
#'     position corresponds to a node of the network. The values of each element 
#'     are 0 or 1. The name of each element corresponds to the name of the node in
#'     that position in the network.
#'     
#' @seealso \code{\link{intToBits}} which this function wraps
#' 
#' @examples
#' validateState(162, cellcycle$genes)
#' validateState('162', cellcycle$genes)
#' validateState(c(0,1,0,0,0,1,0,1,0,0), cellcycle$genes)
#' 
#' @keywords internal
#' @export
validateState <- function(x, node.names) {
  if (is.character(x)) { #convert valid strings integers
    tryCatch({ x <- as.integer(x) }, 
             warning = function(w) stop(paste("Invalid x: ", x))
    ) }
  if (!is.numeric(x)) stop(paste("Invalid x: ", x))
  if (length(x) == 1) { #if integer pass to bool
    x <- as.integer( intToBits(x)[1:length(node.names)] )  
  }
  if ( length(node.names)!=length(x) ) stop("x and node.names should have the same length")
  names(x) <- node.names
  return(x)
}

#' Determine if a node is an input of the network.
#'
#' @param gene node to asses
#' @param net BoolNet network 
#' 
#' @return Boolean is the node an input in the network.
#' 
#' @examples
#' isGeneInput("CycD", cellcycle)
#' 
#' @keywords internal
#' @export
isGeneInput <- function(gene, net) {
    if ( 
        all((which( net$genes == gene )) == net$interactions[[gene]]$input)
        && 
            all(c(0,1) == net$interactions[[gene]]$func)
    ) { return(TRUE) } else return(FALSE)
}
