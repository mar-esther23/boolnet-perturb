###################
#### DATAFRAME ####
###################

#' Convert a BoolNet attractor to dataframe.
#'
#' If Booleans converts a BoolNet attractor to data frame with nodes displayed in Boolean format. First column is the attractor number, second is the number of state inside the attractor, the rest of the columns correspond to each node.
#' If not Boolean it converts a BoolNet attractor to dataframe with properties as columns. The rownames correspond to the int value of each attractor, in the case of cycles the state are joined by sep. Each property of attr$attractors corresponds to a dataframe column. If the property has elements with length > 1 it converts them to a string and joins them with sep.
#' 
#' @param attr BoolNet attractor object
#' @param node.names node names, by default taken from attractor object 
#' @param sep string to join elements with length > 1, default "/"
#' @param Boolean return attractor in Boolean or integer format, default FALSE
#' @return If Boolean=TRUE return dataframe, each column corresponds to the numebr of attractor, state, or node. If Boolean=FALSE return dataframe, each column corresponds to a property of the attractor
#'
#' @examples
#' attr <- getAttractors(cellcycle)
#' attractorToDataframe(attr)
#' #              involvedStates basinSize
#' #1                        162       512
#' #2 25/785/849/449/389/141/157       512
#' 
#' attractorToDataframe(attr, Boolean=TRUE)
#' #   attractor state CycD Rb E2F CycE CycA p27 Cdc20 Cdh1 UbcH10 CycB
#' #1         1     1    0  1   0    0    0   1     0    1      0    0
#' #2         2     1    1  0   0    1    1   0     0    0      0    0
#' #3         2     2    1  0   0    0    1   0     0    0      1    1
#' #4         2     3    1  0   0    0    1   0     1    0      1    1
#' #5         2     4    1  0   0    0    0   0     1    1      1    0
#' #6         2     5    1  0   1    0    0   0     0    1      1    0
#' #7         2     6    1  0   1    1    0   0     0    1      0    0
#' #8         2     7    1  0   1    1    1   0     0    1      0    0
#' 
#' @export
attractorToDataframe <- function(attr, sep="/", node.names=NULL, Boolean=FALSE) {
    if (Boolean) {
          if (is(attr, "AttractorInfo")) {
        if (is.null(node.names)) node.names <- attr$stateInfo$genes
        attr <- sapply(attr$attractors, function(a) paste(a$involvedStates, collapse=sep) )
      }
      if (is.null(node.names)) stop("Invalid node.names")
      if (is.list(attr)) { attr <- unlist(attr) }

      df <- data.frame(matrix(ncol=length(node.names)+2, nrow=0))
      for (i in seq(length(attr))) {
        s <- attr[i]
        if(is.character(s)) s<-unlist(strsplit(s,sep))
        for (j in seq(length(s))) {
          df <- rbind(df, c(attractor=i,state=j,int2binState(s[j],node.names)))
        }
      }
      colnames(df)<-c('attractor','state',node.names)
      return(df)
    }
    
    else {
        # check if valid attr object
        if (!is(attr, "AttractorInfo")) { stop("Error: non-valid attractor") }
        attr <- attr$attractors
        # create properties list, if labeled we will have more
        attr.properties <- vector("list", length(attr[[(1)]]))
        names(attr.properties) <- names(attr[[(1)]])
        #print(attr.properties)

        for (n in names(attr.properties) ) { #create list for each property
            attr.properties[[n]] <- lapply(attr, function(a) a[[n]]) 
            #verify number of elements inside list
            ncol <- max(sapply(attr.properties[[n]], length))
            if ( ncol > 1) { #collapse
                  attr.properties[[n]] <- sapply(attr.properties[[n]], function(a) {
                  paste(as.character(a), collapse=sep)
                  })}
            attr.properties[[n]] <- unlist(attr.properties[[n]])
            #print(attr.properties[[n]])
        }
        return(data.frame(attr.properties, stringsAsFactors=F))
    }
}



#' Convert a list of attractors to a data frame.
#' 
#' Convert a list of BoolNet attractor objects to a data frame. Each property of each attr$attractors corresponds to a dataframe column. Columns are named attrName.propertyName, if the list has no names numbers will be used. If the property has elements with length > 1 it converts them to a string and joins them with sep.
#'
#' @param attr.list list of BoolNet attractor objects
#' @param sep string to join elements with length > 1, default "/"
#' @param returnDataFrame if returnDataFrame='occurrence' returns a df where each column corresponds to a network and the rows to the attractor/label or labels. The values indicate the frequency of the attractor/label
#' if returnDataFrame='basinSize' returns a df where the values indicate the basin size of the attractor/label
#' if returnDataFrame='attrList' returns a list of AttractorInfo objects
#' 
#' @return Dataframe, each column corresponds to a property of the attractor
#' 
#' @examples
#' data(cellcycle)
#' attrs <- list(getAttractors(cellcycle))
#' attractorListToDataframe(attrs)
#' 
#' @keywords internal
attractorListToDataframe <- function(attr.list, sep='/', returnDataFrame=c('occurrence','basinSize'), ...) {
  # Receives a list of BoolNet attractors and return a dataframe
  # Each column is named attrName.propertyName
  returnDataFrame <- match.arg(returnDataFrame)
  # Transform from attr to dataframe
  attr.list <- lapply(attr.list, attractorToDataframe, sep)
  # Verify names exist
  if ( is.null(names(attr.list)) ) names(attr.list) <- 1:length(attr.list)
  # set involvedStates as rowname, delete and rename df columns 
  for (n in names(attr.list)) {
    rownames(attr.list[[n]]) <- attr.list[[n]]$involvedStates #set involvedStates as rowname
    if (returnDataFrame=='occurrence') attr.list[[n]] <- replace(attr.list[[n]],!is.na(attr.list[[n]]),1)
    attr.list[[n]]$involvedStates <- NULL #delete
    names(attr.list[[n]]) <- n # rename df columns
    #print(attr.list[[n]])
  } 
  
  #merge and reduce by rownames
  attr.df <- Reduce(function(x, y){
    df <- merge(x, y, by= "row.names", all=TRUE)
    rownames(df) <- df$Row.names
    df$Row.names <- NULL
    return(df)
  }, attr.list)
  attr.df
}

#' Convert a data frame with nodes displayed in Boolean format to a BoolNet attractor.
#' 
#' Convert a data frame with nodes displayed in Boolean format to a BoolNet attractor. First column is the attractor number, second is the number of state inside the attractor, the rest of the columns correspond to each node.
#'
#' @param Dataframe, see\code{\link{attractorToDataframe}} each column corresponds to the number of attractor, state, or node
#' @param fixedGenes fixedGenes of network
#' 
#' @return attr BoolNet attractor object
#' 
#' @examples
#' > data("cellcycle")
#' > attr <- getAttractors(cellcycle)
#' > attr.df <- attractorToDataframe(attr)
#' > print(dataframeToAttractor(attr.df))
#' 
#' @keywords internal
#' @export
dataframeToAttractor <- function(df, fixedGenes) {
  bin2intState <- function(x){ 
    x <- rev(x)
    sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))
  }
  attractors = vector("list", length = max(df$attractor))
  #df <- df[seq(dim(df)[1],1),]
  for(i in 1:nrow(df)) {
    row <- df[i,]
    n <- row[[1]]
    state <- bin2intState(row[c(-1,-2)])
    attractors[[ n ]] <- c(attractors[[ n ]], state)
  }
  for (i in seq(length(attractors))) {
    l = length(attractors[[i]])
    attractors[[i]] <- list(
      involvedStates = array(attractors[[i]], dim=c(1,l)),
      basinSize = NA
    )  
  }
  node.names <- colnames(df)[c(-1,-2)]
  if (missing(fixedGenes)) {
    fixedGenes <- rep(-1, length(attr$stateInfo$genes))
    names(fixedGenes) <- node.names
  }
  stateInfo = list( genes = node.names, fixedGenes = fixedGenes )
  
  result <- list( stateInfo = stateInfo,attractors = attractors )
  class(result) <- "AttractorInfo"
  result
}
      
      
      

#' Label a list of int attractors and agregate by label
#' 
#' Takes a dataframe where the rownames are states in integer format (cycles are strings joined with sep). The function labels the states and agregates the dataframe using label. Agregate splits the data into subsets by label, computes summary statistics for each, and returns the result in a convenient form. See   \code{\link[stats]{aggregate}}
#'
#' @param df dataframe with states as rownames
#' @param node.names node names of the state, the length must be the same that the state's
#' @param label.rules dataframe with labels (1st col) and rules (2nd col), if more than one rule is true all labels are appendedl the node names present in the rules must be in node.names
#' @param sep string to join elements with length > 1, default "/"
#' 
#' @return Dataframe, each row corresponds to a label and each column corresponds to a property of the original dataframe 
#' 
#' @keywords internal
#' @export
aggregateByLabel <- function(df, node.names, label.rules, sep='/') {
  labels <- lapply(rownames(df), function(state) {
    state <- as.numeric(unlist(strsplit(state, sep)))
    label <- lapply(state, function(s) {
      l <- labelState(s, node.names, label.rules)
    })  
    label <- paste(label, collapse=sep)
  })
  labels <- unlist(labels)
  df[is.na(df)] <- 0
  df <- aggregate(df, list(labels), sum)
  rownames(df) <- df$Group.1
  df$Group.1 <- NULL
  as.data.frame(df)    
}

      
#' Count won or lost values in comparison with those present in a reference column.
#' 
#' Counts how many values were won or lost in comparison between dataframe. It measures ocurrences and then counts the values with sum.
#'
#' @param df dataframe with numeric values
#' @param reference name of column to use as reference
#' 
#' @return axis, axis to obtain new/lost values, see \code{\link[base]{sum}}
#' 
#' @examples
#' df <- data.frame(WT=c(1,2,0,0), a=c(1,2,0,0), 
#'                  b=c(0,2,1,0), c=c(1,0,3,1))
#'countChangesDataframe(df)
#' #     	WT	    a	b	c
#' #     	new	    0	0	1	2
#' #     	lost	0	0	1	1
#'countChangesDataframe(df, axis=1)
#' #     	WT	    a	b	c
#' #     	new	    0	0	2	1
#' #     	lost	1	1	0	0
#' 
#' @keywords internal
#' @export
countChangesDataframe <- function(df, reference='WT', axis=2) {
    df <- df/df
    df[is.na(df)] <- 0
    df <- df-df[[reference]]
    new <- apply(df, axis, function(x) sum(x==1))
    lost <- apply(df, axis, function(x) sum(x==-1))
    rbind(new, lost)
}
