

      
      
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
    if (length(inputs)>0) {
      initial.state <- setStateValues(initial.state, inputs, inputs.values)
    }
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
#' @param dataframe of states state to perturb, default is all attractors.
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


