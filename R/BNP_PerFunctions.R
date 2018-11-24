#' Perturb network functions knock-out (0) or over-expression (1).
#' 
#' Simulates fixed function perturbations (knock-out and over-expression) in a network. Takes a network, a list of gene names and corresponding values to perturb, fixes those values in each network and returns the corresponding attractors for all perturbed networks. By default it  returns all the single node knock-out and over-expression of a network.
#'
#' @param net network to perturb
#' @param genes list of gene names to perturb. To perturb multiple nodes at the same time use a vector inside the list.
#' @param values list of values of the perturbed genes. Knock-out is 0, over-expression is 1.
#' @param names names of the perturbation
#' @param returnDataFrame if returnDataFrame='occurrence' returns a df where each column corresponds to a network and the rows to the attractor/label or labels. The values indicate the frequency of the attractor/label
#'   if returnDataFrame='basinSize' returns a df where the values indicate the basin size of the attractor/label
#'   if returnDataFrame='attrList' returns a list of AttractorInfo objects
#' @param label.rules if label.rules=NULL (default) returns a df with the attractors as rows
#'   if label.rules is a dataframe with labels (1st col) and rules (2nd col), if more than one rule is true all labels are appendedl the node names present in the rules must be in node.names
#' @param sep string to join states in cyclic attractors, default "/"
#' @param verbose report the node being perturbed, useful for larger networks
#' @param ... Further parameters to getAttractors.
#' @return dataframe or list of attractors of the perturbed networks
#' @seealso \code{\link{fixGenes}} \code{\link{attractorListToDataframe}} 
#' @export
#' @examples
#' data(cellcycle)
#' # All single gene knock-out and over-expression of cellcycle network
#' perturbNetworkFixedNodes(cellcycle)
#' 
#' # Cellcycle network CycD over-expression
#' perturbNetworkFixedNodes( cellcycle, list("CycD"), list(1) )
#' 
#' # Cellcycle network CycD knock-out with Rb over-expression
#' perturbNetworkFixedNodes(  cellcycle, list(c("CycD", "Rb")), list(c(0,1)) , list("CycD-RB+") )
#' 
#' # All CycD double gene knock-outs
#' double.KO <- lapply(cellcycle$genes[-1], function(gene) c("CycD", gene))
#' number.double.KO <- length(double.KO)
#' perturbNetworkFixedNodes(  net=cellcycle, genes=double.KO, 
#'                         value=as.list(rep(0, number.double.KO))  ) 
perturbNetworkFixedNodes <- function(net, genes, values=0, names, 
                                     returnDataFrame=c('occurrence','basinSize','attrList'),
                                     label.rules=NULL, sep='/', verbose=FALSE, ...) {
  if (!is(net, "BooleanNetwork")) { stop("Error: non-valid network") }
  returnDataFrame <- match.arg(returnDataFrame)
  # generate genes to evaluate
  if (missing(genes)) { #Default, evaluate all single KOver
    genes <- c(NA, net$genes, net$genes)
    values <- c(NA, rep(0, length(net$genes)), rep(1, length(net$genes)))
  }
  # check lengths
  if (length(values)==1) values <- rep(values, length(genes))
  if (length(genes)!=length(values)) stop("genes and values must have the same number of elements, or values must have 1 element!")
  if (missing(names)) { #create names
    names <- c()
    for (i in seq(length(genes))) {
      n <- paste(paste(genes[[i]],collapse=''), paste(values[[i]],collapse=''),sep='_')
      names <- c(names,n)
    }
    names[match('NA_NA', names)] <- 'WT'
  }
  if (length(values)==1) values <- rep(values, length(genes))
  if (length(genes)!=length(names)) stop("genes and names must have the same number of elements!")
  
  # Calculate mutants
  mutants <- list()
  for (i in 1:length(genes)) {
    if (verbose) print( unlist(c(i,names[i],genes[i],values[i])) )
    not.WT <- ! (is.null(genes[i]) || is.na(genes[i]))
    #         print(not.WT)
    if ( not.WT ) { net <- fixGenes(net, unlist(genes[i]), unlist(values[i])) }
    mutants[[i]] <- getAttractors(net, ...)
    if ( not.WT ) net <- fixGenes(net, unlist(genes[i]), -1)
  }
  names(mutants) <- names
  # convert attractors to dataframe
  if (returnDataFrame!='attrList') {
    mutants <- attractorListToDataframe(mutants, sep=sep, returnDataFrame=returnDataFrame)
    if (!is.null(label.rules)) mutants <- aggregateByLabel(mutants, net$genes, label.rules, sep=sep)
  }
  mutants
}