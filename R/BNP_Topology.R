##################
#### TOPOLOGY ####
##################

#' Returns the topology with signs of a network
#' 
#' Takes a BoolNet network
#' Returns a DataFrame of the form: Source, Target, Interaction
#' Interactions types can be:
#'     '+'  : Mandatory, positive, unambiguous
#'     '-'  : Mandatory, negatice, unambiguous
#'     'MA' : Mandatory, ambiguous
#'     'NR' : Non-functional regulation
#' 
#' @param net BoolNet function
#' @return df.graph dataframe with source and target nodes and the type of interactions
#' 
#' @examples
#' data(cellcycle)
#' getNetTopology(cellcycle)
#' 
#' @export
getNetTopology <- function(net) {
  if (!is(net, "BooleanNetwork")) { stop("Error: non-valid network") }
  df <- data.frame(matrix(vector(), 0, 3,
                   dimnames=list(c(), c("Source","Target","Interaction"))),
                   stringsAsFactors=F)
  for (i in net$genes) {
    for (j in net$genes) {
      sign = getInteractionSign(j, i, net)
      if (!is.null(sign)) {
        de = data.frame(j,i,sign, stringsAsFactors=F)
        names(de) <- c("Source","Target","Interaction")
        df = rbind(df, de)
      }
    }
  }
  return(df)
}



#' Determine the sign of the interaction between two elements of a network
#' 
#' Returns a dataframe with the following interaction types:
#'     '+'  : Mandatory, positive, unambiguous
#'     '-'  : Mandatory, negatice, unambiguous
#'     'MA' : Mandatory, ambiguous
#'     'NR' : No regulation
#'     NULL : No interaction
#'     
#' @param source source node name
#' @param target target node name
#' @param net BoolNet function  
#' @return interaction str interaction type between source and target


#' Get sign of interaction
#' 
#' @keywords internal
#' @export
getInteractionSign <-  function(source, target, net) {
  if (! source %in% net$genes) stop(paste(source, "is not in network"))
  if (! target %in% net$genes) stop(paste(target, "is not in network"))
  # get target interaction
  inter = net$interactions[[target]]
  # check if interaction exists
  index <- match(source, net$genes)
  if (! index %in% inter$input) return(NULL)
  index = match(index, inter$input)
  jump  = 2^(length(inter$input)-index)
  #print(c(source, index, jump))
  #print(inter$func)
  
  check = seq(length(inter$func))
  sign = c()
  #print(check)
  while (length(check)>0) {
    i = check[1]
    #print(c('index',i, i+jump,' ',inter$func[i], inter$func[i+jump]))
    if (inter$func[i] < inter$func[i+jump]) sign = c(sign,'+')
    if (inter$func[i] > inter$func[i+jump]) sign = c(sign,'-')
    check = check[! check %in% c(i, i+jump)]
    #print(check)
  }
  #print(sign)
  if (length(sign)==0)  return('NA')
  if (all(sign=='+'))   return('+')
  if (all(sign=='-'))   return('-')
  return('MA')
}
