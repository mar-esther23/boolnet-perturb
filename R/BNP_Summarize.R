
#' Summarize attractors of a dataframe ignoring certain nodes.
#' 
#' ...
#' 
#' @param attr BoolNet attractor object
#' @param genes	 list of gene names to ignore, by default it detects inputs. To perturb multiple nodes at the same time use a vector inside the list.
#' @param returnBasin default FALSE, return a column with the combined basin size
#' @return dataframe, each column corresponds to the number of attractor, state, or node. If Boolean=FALSE return dataframe, each column corresponds to a property of the attractor.
#' 
#' @examples
#' data(netTh17Treg)
#' attr <- getAttractors(netTh17Treg)
#' summarizeAttractors(attr)
#' 
#' @export
summarizeAttractors <- function(attr, genes, returnBasin=F) {
  
}