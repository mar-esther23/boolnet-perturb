######################
####   LABELING   ####
######################

#' Label a state.
#' 
#' Labels a binary state using a set of labelling rules.
#' The elements of the state correspond the network node names. If a rule is satisfied the corresponding label is appended.
#'
#' @param state binary state to label
#' @param node.names node names of the state, the length must be the same that the state's
#' @param label.rules dataframe with labels (1st col) and rules (2nd col), if more than one rule is true all labels are appendedl the node names present in the rules must be in node.names
#' @param sep string to separate the labels when more than one can be applied to the state.
#' @return String corresponding to the label of the state.
#' 
#' @examples
#' labelState(state, net$genes, label.rules, sep='')
#' 
#' @export
labelState <- function(state, node.names, label.rules, sep='') {
    labels <- label.rules[[1]]
    rules  <- label.rules[[2]]
    state <- validateState(state, node.names)
    label = c()
    for (j in 1:length(rules)) { #evaluate rules
        # create string with function
        f <- paste('function(', paste(node.names, collapse=','), ') { if (',
                   rules[j], ') \'', labels[j], '\' }' , sep='')
        f <- eval(parse(text=f)) # evaluate string to create function
        label <- append(label, do.call(f, as.list(state))) # apply function
    }
    # format label
    if (is.null(label)) { label <-c('NoLabel')}
    label <- paste(label, collapse=sep)
}

#' Label attractor.
#' 
#' Returns a list of labels for a the involved states in a BoolNet attractor object using a set of labelling rules, returns a list of labels.
#'
#' @param attr Boolnet attractor object
#' @param label.rules dataframe with labels (1st col) and rules (2nd col), if more than one rule is true all labels are appendedl the node names present in the rules must be in node.names
#' @param sep string to separate the labels when more than one can be applied to the state, 
#'        if NULL it will return the label of each state in a cycle separately
#'        if string it will paste the states of a cycle with sep
#' @return List of strings corresponding to the label of the attractor, if an attractor has multiple states it will return a list of strings for that state.
#' @seealso \code{\link{labelState}} which this function wraps
#' 
#' @examples
#' attr <- getAttractors(net$genes)
#' labelAttractors(attr, net$genes, labels, rules, sep='')
#' 
#' @export
labelAttractors <- function(attr, label.rules, node.names=NULL, sep="/") {
  # takes an attractors object created by BoolNet
  # returns a list of the labels for each attractor in order.
  # If an attractor has multiple states it will return a label for each state.
  if (!is(attr, "AttractorInfo")) { stop("Error: non-valid attractor") }
  node.names <- attr$stateInfo$genes
  res <- list()
  for (i in 1:length(attr$attractors)) {
      label <- sapply(attr$attractors[[i]]$involvedStates, function(state) {
      state <- int2binState(state, node.names) #state to binary
      l <- labelState(state, node.names, label.rules, sep='') #label
    })
    if (!is.null(sep)) { label <- paste(label, collapse=sep) }
    res <- append(res, list(label))
  }
  unlist(res)
}


#' Label state list from dataframe column 
#' 
#' @keywords internal
#' @export
labelList <- function(states, label.rules, node.names=NULL, sep='/') {
    labels <- sapply(unlist(states), function(state) {
        state <- strsplit(as.character(state), sep)
        #print(c(state))
        state <- as.numeric(unlist(state))
        label <- lapply(state, function(s) labelState(s, node.names, label.rules)) 
        #print(label)
        label <- paste(label, collapse=sep)
        #print(label)
        label
    })
    labels
}



