#' Derrida analysis.
#' 
#' Derrida curves are used to determine the dynamic behaviour of Boolean networks. We take two states with Hamming distance h in the time t, and determine the hamming distance in time t+1. This process is repeated multiple times. If the network is chaotic, the curve tends to be over the diagonal. If the network is ordered it will be below the diagonal.
#' 
#' @param net BoolNet network to perturb
#' @param repetitions number of repetitions
#' 
#' @return vector with average hamming distance in t+1, names are hamming distance in t
#' 
#' @examples
#' data(cellcycle)
#' derrida(cellcycle)
#' 
#' @export
derrida <- function(net, repetitions=1000) {
  #states <- getAttractors(net, method="sat.exhaustive")
  repetitions = round(repetitions/length(net$genes))
  res = 1:length(net$genes)
  for (i in res) { 
    hamming <- sapply(1:repetitions, function(n) {
      state <- sample.int(2^length(net$genes),1) #random initial state
      state <- validateState(state, net$genes)
      genes <- sample(net$genes,i) #random nodes to perturb
      ori <- stateTransition(net,state)
      per <- perturbState(net, state, genes, result = "nextState")
      sum(ori!=per) #hamming between two states
    })
    res[[i]] <- mean(hamming)
  }
  res <- c(0, res)
  names(res) <- 0:length(net$genes)
  res
}
