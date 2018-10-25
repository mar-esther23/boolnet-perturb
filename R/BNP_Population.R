#################################
####   POPULATION ANALYSIS   ####
#################################

#' Analyze and graph the average of multiple asynchronous simulations for an specific number of time steps
#' 
#' @param net Boolnet network
#' @param tStates Number of time steps running each simulation, default value set to 1000.
#' @param simulations Number of simulations equivalent to "population size", default value set to 1000.
#' @param new.nodes Nodes to fix with an specific value at the initial state.
#' @param new.values Values for nodes in new.nodes.
#' @param gene Gene to perturb at an specific time-step (pertTime), default value set NULL.
#' @param duration Time of perturbation, default value set to 1. 
#' @param value Value of perturbed gene, default value set to 1.
#' @param pertTime Specific time-step at which the gene will be perturbed, default value set to 400.
#' @param onePert If this is TRUE the perturbation will only occur one time during the simulation. If FALSE, de perturbation wil occur periodically, default value set to TRUE.
#' 
#' @return Data frame with the average of the simulations for each node in each time step
#' 
#' @export
#' 
#' @examples 
#' data(netTh17Treg)
#' pop_Th17 <- popSimulation(netTh17Treg, new.nodes = c("IL2e","IL21e","TGFBe"), new.values = c(1,1,1), gene ="TGFBe", value = 0)
#' 

populationSimulation <- function(net, tStates=1000, simulations=1000, new.nodes=list(), new.values=list(), gene=NULL, duration = 1, value = 1, pertTime=400, onePert = TRUE) {
  nodesNames <- net$genes
  noNodes <- length(nodesNames)
  for(fn in new.nodes){
    if(is.na(match(fn,nodesNames))==TRUE){
      stop('Node in new.nodes do not match with any gene in network')
    }
  }
  if(length(gene)!=0){
    if(is.na(match(gene,nodesNames))==TRUE){
      stop('Gene name do not match with any gene in network')
    }
  }
  if(duration > pertTime){
    stop('Duration must be smaller or equal than pertTime')
  }
  pulse <- c() # Make a list with the time steps where the perturbation must occur
  pt = tStates/pertTime
  if(onePert==TRUE){
    pt = 1
  }
  for(pT in 1:pt){ 
    for(b in 1:duration){
      pulse <- c(pulse, ((pT*pertTime)+b))
    }
  }
  itP = 0 # Number of iteration
  s <- match(gene, nodesNames) #aquire the index of the node for stimulation
  final <- matrix(0, nrow=tStates+1, ncol=noNodes) #create an empty matrix for the output
  Porct1 = 0 #Tracker for the perecentage of simulation progress
  cat('START (% progress) ')
  for(i in 1:simulations) {
    x <- sample(0:1, noNodes, replace=T) #Create a random initial state of the network
    if(length(new.nodes) != 0){ #Set the values to new.nodes
      x[match(new.nodes,nodesNames)] <- new.values
    }
    finalTemp <- matrix(0, nrow=tStates+1, ncol=noNodes) #Temporal matrix to obtain the average activation values
    finalTemp[1,] <- x
    for(j in 2:(tStates+1)) {
      if(TRUE%in%(j==pulse)){ #Evaluate if the time-step correspond to some value in pulse (time-steps for perturbation)
        x[s] <- value #Change value of gene 
      }
      x <- as.vector(stateTransition(net, state = x, type = 'asynchronous'))
      finalTemp[j,] <- x
    }
    final = final + finalTemp
    itP = itP + 1
    Porct2 = round((itP/simulations)*100)
    if (Porct2 != Porct1) {
      cat('... ')
      cat(Porct2)
    }
    Porct1 = Porct2
  }
  finalGraph = final/itP
  rownames(finalGraph) <- c(0:tStates)
  colnames(finalGraph) <- nodesNames
  cat(' END')
  listReturn <- list('totalData'=round(finalGraph,5), 'perturbedGene'=c(gene,' = ', value))
  return(listReturn)
}

##############################################################################################
#' @param net Boolnet network
#' @param popSim_output Output object from popSimulation
#' @param genes List with genes to plot (to many genes will be difficult to visualize). Default, graph all.
#' 
#' @export
#' 
#' @examples 
#' data(netTh17Treg)
#' pop_Th17 <- popSimulation(netTh17Treg, new.nodes = c("IL2e","IL21e","TGFBe"), new.values = c(1,1,1), gene ="TGFBe", value = 0)
#' plot_popSimulation(netTh17Treg, pop_Th17)
#' 

plotPopulationSimulation <- function(net, popSim_output, genes=NULL){
  #if(  !("reshape2" %in% (.packages()))  ) warning("reshape2 is not attached")
  if(  !("ggplot2" %in% (.packages()))  ) warning("ggplot2 is not attached")
  
  AnalyzeNetTable <- popSim_output$totalData
  if(length(genes)==0){
    genes = net$genes
  }
  for(g in 1:length(genes)){
    if(is.na(match(genes[g],net$genes))==TRUE){
      stop('Gene name do not match with any gene in network')
    }
  }
  for(a in 1:length(genes)){
    b <- match(genes, net$genes)
    TimeSteps <- as.numeric(rownames(AnalyzeNetTable))
    dF <- data.frame(TimeSteps,AnalyzeNetTable[,b])
    dF2 <- reshape2::melt(data = dF, id.vars = "TimeSteps")
    plot <- ggplot(data = dF2, aes(x = TimeSteps, y = value, colour = variable)) + geom_line() + ggtitle(paste(genes, collapse = ", ")) + ylab("Average activation value")  + 
      geom_line(size = 1)
    return(plot)
  }
}
