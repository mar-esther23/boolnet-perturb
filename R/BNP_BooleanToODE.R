# References:
# Alvarez-Buylla et al (2008). Floral morphogenesis: stochastic explorations of a gene network epigenetic landscape. Plos one, 3(11), e3626.
# Sanchez-Corrales et al (2010). The Arabidopsis thaliana flower organ specification gene regulatory network determines a robust differentiation process. Journal of theoretical biology, 264(3), 971-983.
# Villarreal et al (2012). General theory of genotype to phenotype mapping: derivation of epigenetic landscapes from N-node complex gene regulatory networks. Physical review letters, 109(11), 118102.
# Davila-Velderrain (2015). Reshaping the epigenetic landscape during early flower development: induction of attractor transitions by relative differences in gene decay rates. BMC systems biology, 9(1), 1.

#' Translate from Boolean Network to Continuous ODE Network.
#' 
#'  The function translate the Boolean Network to a Determinist Continuous Network 
#'  composed by two principal parts. The first one containes the translation from 
#'  classic to fuzzy logic using one of two methods. The second one contains an 
#'  ODE-based description of the rate of changes for each of the nodes. 
#' IMPORTANT: This function call the perl script "ConvertTocontinuos.pl". 
#'
#' @param net BoolNet network
#' @param logic is a character string specificating either if Zadeh or Probabilistic 
#' logic is employed during the translation. Default: Zadeh.
#' @param eq is a character string specificating either if the equation proposed 
#' by Sanchez-Corrales et al., 2010 or Villarreal et al., 2012 is employed by 
#' describe the rate of change for the nodes in the network. 
#' @param sep between targets and factors in file, default ","
#' @param mutants list of genes to mutate, the derivate will be 0
#' @param keep.input if TRUE set the derivate of the inputs to 0
#' 
#' 
#' @return Returns a list with func (ode function), parameters and example state and time.
#'         The parameters are a vector where h=10, w=0.5 and decay rate = 1.
#' 
#' @examples
#' library(deSolve)
#' data(cellcycle)
#' 
#' net.ode <- booleanToODE(cellcycle)
#' ode(func = net.ode$func, 
#'     parms = net.ode$parameters, 
#'     y = net.ode$state, 
#'     times = seq(0, 5, 0.1))
#' 
#' @export
booleanToODE <- function(net, logic = "Zadeh", eq = "SQUAD", sep = ",", 
                         mutants = c(), keep.input = FALSE){
  
  saveNetwork(net, "./net.csv") #write and
  net <- "./net.csv"
  
  if(logic != "Zadeh" && logic != "Probabilistic"){
    stop(paste("Invalid logic value: logic must be either 'Fuzzy' or 'Probabilistic' "))
  }
  if(eq != "SQUAD" && eq != "Villarreal" && eq != "HillCube"){
    stop(paste("Invalid eq value: eq must be either 'LuisMendoza' or 'CarlosVillarreal'  or 'HillCube' "))
  }
  
  if(eq == "HillCube"){
    net = loadNetwork(net)
    
    numGenes = length(net$genes)
    
    sink("HillCubeInfo.txt")
    for(i in 1:numGenes){
      gene     = net$genes[i]
      regs     = net$genes[net$interactions[[i]]$input]
      regs     = paste(regs, sep ='', collapse = ",")
      pos.true = which(net$interactions[[i]]$func == 1) - 1
      pos.true = paste(pos.true, sep ='', collapse = ",")
      info = paste(gene, regs, pos.true, sep = '\t', collapse = "")
      cat(info)
      cat("\n")
    }
    sink()
    
    net = "HillCubeInfo.txt"
    
  }
  
  
  if(length(mutants) != 0){
    mutants = paste(mutants, sep = ",", collapse = ",")
  } else {
    mutants = "."
  }
  if(keep.input){
    input = "1"
  } else{
    input = "0"
  }
  perl.scrip.path = paste0(system.file(package = "BoolNetPerturb"),"/perl/ConvertToContinuos.pl")
  cmd = paste("perl", perl.scrip.path, net, logic, eq, sep, mutants, input, "> cont.r")
  #print(cmd)
  system(cmd)
  source("./cont.r") #load
  system("rm ./cont.r") #and delete
  system("rm ./net.csv") #and delete
  return(BoolODE)
}
