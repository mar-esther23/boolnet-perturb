#####################
####   PERTURB   ####
#####################

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
#' @param ... Further parameters to getAttractors.
#' @return dataframe or list of attractors of the perturbed networks
#' @seealso \code{\link{fixGenes}} \code{\link{attractorsToDataframe}} 
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
                                     label.rules=NULL, sep='/', ...) {
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
    #print( unlist(c(i,names[i],genes[i],values[i])) )
    not.WT <- ! (is.null(genes[i]) || is.na(genes[i]))
    #         print(not.WT)
    if ( not.WT ) { net <- fixGenes(net, unlist(genes[i]), unlist(values[i])) }
   mutants[[i]] <- getAttractors(net, ...)
   if ( not.WT ) net <- fixGenes(net, unlist(genes[i]), -1)
  }
  names(mutants) <- names
  # convert attractors to dataframe
  if (returnDataFrame!='attrList') {
    mutants <- attractorsToDataframe(mutants, sep=sep, returnDataFrame=returnDataFrame)
    if (!is.null(label.rules)) mutants <- aggregateByLabel(mutants, net$genes, label.rules, sep=sep)
  }
  mutants
}

                           
#' Compare atractors with synchronous and asynchronous updating.
#' 
#' Takes a set of BoolNet attractor or states and verifies if they can be recovered with asynchronous updating. See \code{\link[BoolNet]{getAttractors}} with options type='async' and method="chosen".
#' 
#' @param net BoolNet network to perturb
#' @param attr BoolNet AttractorInfo object
#' @param label.rules dataframe with labels (1st col) and rules (2nd col), if more than one rule is true all labels are appendedl the node names present in the rules must be in node.names
#' @param ... Further parameters to  \code{\link{getAttractors}}.
#' @return dataframe comparing synchronous and asynchronous updating
#' @seealso \code{\link[BoolNet]{getAttractors}}
#' @export
#' @examples
#' data(cellcycle)
#' attr <- getAttractors(cellcycle)
#' verifyAttractorSyncronousVsAsyncronous(cellcycle, attr)
verifySyncronousVsAsyncronous <- function(net, attr, label.rules=NULL, ...) {
    if (!is(net, "BooleanNetwork")) { stop("Error: non-valid network") }
    if (missing(attr)) attr <- getAttractors(net,type = 'sync')
    if (!is(attr, "AttractorInfo")) { stop("Error: non-valid attractor") }
    # convert attr to list of vectors for getAttractors
    l <- attractorToDataframe(attr, Boolean=T)
    l["attractor"] <- NULL
    l["state"] <- NULL
    l <- simplify2array(l)
    l <- split(l, seq(nrow(l)) )
    attr.async <- getAttractors(net,type='async', 
                                method="chosen",startStates=l, ...)
    if (is.null(label.rules)) {
        attr <- attractorToDataframe(attr)$involvedStates
        attr.async <- attractorToDataframe(attr.async)$involvedStates
    } else {
        attr <- labelAttractors(attr, label.rules)
        attr.async <- labelAttractors(attr.async, label.rules)
    }
    attr <- table(attr)
    attr.async <- table(attr.async)
    #print(attr)
    #print(attr.async)
    res <- merge(data.frame(sync=attr,row.names=1),
            data.frame(async=attr.async,row.names=1),
            by="row.names", all=TRUE) 
    res <- rbind(res, data.frame(Row.names="Total",t(colSums(res[,-1], na.rm=T))))
    rownames(res) <- res$Row.names
    res$Row.names <- NULL
    res
}

      
      
#' Perturb network state or trajectory.
#' 
#' Modifies the values of a state and returns the resulting attractor or trajectory. If time=NULL the perturbation will be fixed permanently, if time=n the perturbation will last n time steps and then the rules will return to their original values.
#'
#' @param net BoolNet network to perturb
#' @param state state to perturb, default is a random state.
#' @param genes genes to perturb, default is a random gene.
#' @param values value of perturbed genes, default is bitflip of target genes.
#' @param time time of perturbation, default is 1. If time=NULL the perturbation will be fixed permanently, if time=n the perturbation will last n time steps and then the rules will return to their original values.
#' @param result select 'attractor','nextState','trajectory'
#' @param all.data returns list with initial state, perturb node, and result. states will be simplified as integer or string with sep='/' for cycles
#' @param int return attractor as int with sep="/"
#' @param ... Further parameters to \code{\link{stateTransition}} and \code{\link{getAttractors}}.
#' @return If result='attractor' returns the attractor the state reached after the perturbation
#'         if result=''nextState' returns the next state reached after the perturbation
#'         if result='trajectory' return the trajectory after perturbation (only synchronous update)
#' @seealso \code{\link{stateTransition}}, \code{\link{getAttractors}}, \code{\link{perturbTrajectories}}
#' @export
#' @examples
#' # Perturb two nodes of a state for one time step, return attractor object
#' > perturbState(cellcycle, 162, c('CycD', 'Rb'))
#'  |--<---------|
#'  V            |
#'  0100010100   |
#'  V            |
#'  |-->---------|
#' 
#' # Perturb two nodes of a state for one time step, return nextState
#' > perturbState(cellcycle, 162, c('CycD', 'Rb'), result='nextState')
#'   CycD     Rb    E2F   CycE   CycA    p27  Cdc20   Cdh1 UbcH10   CycB 
#'      0      0      1      0      0      0      0      1      0      0 
#' 
#' # Perturb an state by permanetly fixing the value of nodes, return the trajectory
#' > perturbState(cellcycle, 162, c('CycD', 'Rb'), result = 'trajectory')
#'   CycD Rb E2F CycE CycA p27 Cdc20 Cdh1 UbcH10 CycB
#' 1    0  1   0    0    0   1     0    1      0    0
#' 2    1  0   0    0    0   1     0    1      0    0
#' 3    0  0   1    0    0   0     0    1      0    0
#' 4    0  1   1    1    1   1     0    1      0    0
#' 5    0  1   0    0    0   0     0    1      0    0
#' 6    0  1   0    0    0   1     0    1      0    0
#' 
#' # Perturb random node for 3 time steps
#' perturbState(cellcycle, time=3)
#' 
perturbState <- function(net, state, genes, time=1, result=c('attractor','nextState','trajectory'), all.data=F, int=F, ...) {    
  #' Modify state, set new value to target nodes.
  setStateValues <- function(state, new.nodes, new.values) {
    for (i in 1:length(new.nodes)) { state[new.nodes[i]] <- new.values[i] }
    state
  }
  # validate and select random if default
  if (!is(net, "BooleanNetwork")) { stop("Error: non-valid network") }
  if (missing(state)) state <- sample.int(2^length(net$genes),1) #random initial state
  state <- validateState(state, net$genes) #validate state
  if (missing(genes)) genes <- sample(net$genes,1) # random gene
  result <- match.arg(result)
  # initialize things
  values <- as.integer(!state[genes]) #bitflip target genes
  initial.state <- state 
  if (result=='trajectory') path.perturbed <- list(initial.state) #save initial state in path
  
  if (!is.null(time)) { # if transient perturbation
    # determine original value of inputs
    inputs <- unlist(sapply(net$genes, function(g) {
      if (isGeneInput(g,net)) return(T) else F
    }))
    # We need to save the inputs, because the default input function is input=input
    # So if we don't explicitly return them to the original value 
    # they will stay in the modified value !
    inputs.values <- state[inputs]
    inputs <- names(inputs.values)
    # apply perturbation
    net <- fixGenes(net, genes, values) 
    initial.state <- setStateValues(state, genes, values)
    for (t in seq(time)) { #iterate n times
      if (result=='trajectory') path.perturbed <- append(path.perturbed, list(initial.state))
      initial.state <- stateTransition(net, initial.state, ...)
    }
    # recover original network and inputs
    net <- fixGenes(net, genes, -1) 
    initial.state <- setStateValues(initial.state, inputs, inputs.values) 
  } else { # if fixed perturbation
    # apply perturbation
    net <- fixGenes(net, genes, values) 
    # change the values according to perturbation
    initial.state <- setStateValues(state, genes, values)
    #if (result) path.perturbed <- append(path.perturbed, list(initial.state))
  }
  if (result=='nextState') {
      pert.state <- stateTransition(net,initial.state, ...)
      if (! all.data) return(pert.state)
      else if (int) {
          state <- bin2intState(state)
          pert.state <- bin2intState(pert.state)
      }
      return(list(
          initial=state,final=pert.state,
          genes=genes,values=values  ))
  }
  
  if (result=='attractor') { #just return the attractor
      attr <- getAttractors(net, startStates = list(initial.state), returnTable=F, ...)
      if (! all.data) return(attr)
      else if (int) {
          state <- bin2intState(state)
          attr <- attr$attractors[[1]]$involvedStates
          if(length(attr)>1) attr <- paste(attr, collapse="/")
      }
      return(list(
          initial=state,final=attr,
          genes=genes,values=values  ))
  }
  
  # calculate the rest of the trajectory
  path <- getPathToAttractor(net, initial.state)
  
  #if (is.null(time)) return(path)  # if fixed perturbation and trajectory
  # join both paths
  path.perturbed <- sapply(path.perturbed, function(f) f) # simplify path.perturbed
  path <- t(path)
  path.res <- t(cbind(path.perturbed, path)) #join lists
  rownames(path.res) <- c(1:(dim(path.res)[1]))
  if (! all.data) return(path.res)
  return(list(
          path=path.res,
          genes=genes,values=values  ))
}

#' Get cell fate map of a network.
#' 
#' Runs all perturbations of time N for all states and all genes, returns a dataframe with initial states, final attractor, pertubed gene and value. See \code{\link{perturbState}}
#'
#' @param net BoolNet network to perturb
#' @param states state to perturb, default is all attractors.
#' @param genes genes to perturb, default is all genes.
#' @param values value of perturbed genes, default is bitflip of target genes.
#' @param time time of perturbation, default is 1. If time=NULL the perturbation will be fixed permanently, if time=n the perturbation will last n time steps and then the rules will return to their original values.
#' @param label.rules dataframe with labels (1st col) and rules (2nd col), if more than one rule is true all labels are appendedl the node names present in the rules must be in node.names
#' @return dataframe with initial state, final attractor, pertubed gene and value. If label.rules is given it will label all states according to the rules.
#' @seealso \code{\link{perturbState}}
#' @export
#' @examples
#' data(cellcycle)
#' cellFateMap(cellcycle, time=NULL)
#' data(cellcycle)
#' cellFateMap(cellcycle, time=NULL)
#'       
cellFateMap <- function(net,states,genes,time=1, label.rules, ...) {    
    # validate and select random if default
    if (!is(net, "BooleanNetwork")) { stop("Error: non-valid network") }
    if (missing(states)) { #run for all attractors
        states <- getAttractors(net, method="sat.exhaustive")
        states <- attractorToDataframe(states, Boolean=T)
        states["attractor"] <-NULL
        states["state"] <- NULL
    }
    if (missing(genes)) {genes = net$genes}
    res <- apply(states, 1, function(s) {
        r <- sapply(genes, function(g) {
            perturbState(net,s,g,time,all.data=T,int=T) 
        })
        t(data.frame(r))
    })
    res <- do.call("rbind", res)
    row.names(res) <- NULL
    res <- data.frame(res, row.names=NULL)
    
    if (!missing(label.rules)) { #label
        res["initial"] <- suppressWarnings( 
            labelList(res["initial"], label.rules,net$genes) )
        res["final"] <- suppressWarnings( 
            labelList(res["final"], label.rules,net$genes) )
    }
    
    return(res)
}


derrida <- function(net, repetitions=1000) {
  #states <- getAttractors(net, method="sat.exhaustive")
  repetitions = round(repetitions/length(net$genes))
  res = 1:length(net$genes)
  for (nodes.to.perturb in res) { 
    for (n in 1:repetitions) {
      state <- sample.int(2^length(net$genes),1) #random initial state
      state <- validateState(state, net$genes)
      genes <- sample(net$genes,nodes.to.perturb) #random nodes to perturb
      
      attr.ori <- getAttractors(net,startStates = list(state) )
      attr.ori <- attr.ori$attractors[[1]]$involvedStates
      attr.per <- perturbState(net,state,genes)
      attr.per <- attr.per$attractors[[1]]$involvedStates
      
      hamming <- bitwXor(attr.ori,attr.per)
      hamming <- sum(as.integer(intToBits(hamming)))
      hamming
    }
    
  }
}
