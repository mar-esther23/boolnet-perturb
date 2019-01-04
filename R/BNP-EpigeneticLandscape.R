#library(dplyr)
#library(BoolNet)
#library(reshape2)


#' Calculate the epigenetic landscape of a network
#' 
#' ...
#'  
#' @param net BoolNet network to perturb
#' @param p probability of noise in the simulation
#' @param returnDataFrame returns a df where ...
#' @return The transition probabilities between states.
#' @export
#' @examples
#' data(netTh17Treg)
#' epigeneticLandscape(netTh17Treg, p = 0.01, returnDataFrame = T)
epigeneticLandscape = function(net, p = 0.01, returnDataFrame = F){
  HammingDistance = function(x, len = length(net$genes)){
    z = bitwXor(x[1],x[2])
    z = BoolNet:::dec2bin(z, len = len)
    z = sum(z)
    return(z)
  }
  
  Normalize = function(x){
    x = x/sum(x)
    return(x)
  }
  
  attr = getAttractors(net)
  
  # Obtain the number of S(t) that lead to each S(t + 1) states.
  kin = table(attr$stateInfo$table)
  # Obtain the S(t + 1) states.
  t.plus = as.numeric(names(kin))
  t = seq(0, 2**length(net$genes) - 1)
  
  kin = rep(kin, each = length(t))
  t.plus = rep(t.plus, each = length(t))
  t = rep(t, length.out = length(t.plus))
  
  states = cbind(t = t, tplus = t.plus)
  # Obtain the Hamming distances against all S(t+1) states and all the S(t)
  HDistance = apply(states, 1, HammingDistance)
  numgenes = length(net$genes)
  # Convert the Hamming distances into transtion probabilities.
  prob = p**(HDistance)*(1-p)**(numgenes - HDistance)
  
  # Associate the S(t) states to its respective attractor basin.
  t.basin = attr$stateInfo$attractorAssignment
  t.basin = rep(t.basin, length.out = length(t.plus))
  basin = attr$stateInfo$attractorAssignment
  # Associate the S(t+1) states to its respective attractor basin.
  t.plus.basin = basin[t.plus + 1]
  
  data = data.frame(t.plus = t.plus, t = t, prob = kin*prob, BasinOUT = t.plus.basin, BasinIN = t.basin)
  subset(data, data$BasinOUT == 1 & data$BasinIN == 1)
  
  # Summarize transition probabilities among individual states into 
  # transition probabilities among attractor basins.
  data = data %>% group_by(BasinOUT, BasinIN)
  data = summarize(data, prob = sum(prob))
  data = data %>% group_by(BasinOUT)
  
  data$prob = data$prob/sum(data$prob)
  data = mutate(data, probability = Normalize(prob))
  data = data[,-3]
  
  if(!returnDataFrame){
    data = dcast(data, BasinOUT~BasinIN, value.var = "probability")
    data = data[,-1]
  }
  return(data)
}
