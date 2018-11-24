#' Summarize attractors of a dataframe ignoring certain nodes.
#' 
#' Summarizes the attractors as a dataframe, ignoring nodes that are not relevant for the analysis (input nodes by default). Returns a dataframe in binary format. It can return the combined basin size of the attractors. The function does not summarize cycles to avoid mistakes.
#' 
#' @param attr BoolNet attractor object
#' @param genes	 list of gene names to ignore, by default it detects inputs. To perturb multiple nodes at the same time use a vector inside the list.
#' @param net net to detect inputs from
#' @param returnBasin default FALSE, return a column with the combined basin size
#' @return dataframe, each column corresponds to the number of attractor, state, or node. If Boolean=FALSE return dataframe, each column corresponds to a property of the attractor.
#' 
#' @examples
#' data(netTh17Treg)
#' attr <- getAttractors(netTh17Treg)
#' summarizeAttractors(attr, net=netTh17Treg, returnBasin=T)
#' #    attractor state IL2 RORGT STAT3 FOXP3 TGFB basinSize
#' # 1          1     1   0     0     0     0    0         7
#' # 2          2     1   1     0     0     0    0        21
#' # 3          3     1   0     0     1     0    0        60
#' # 4          4     1   0     0     0     0    1        27
#' # 5          7     1   1     0     0     1    1        41
#' # 6         14     1   0     1     1     0    1        62
#' # 7         20     1   0     0     0     0    0        16
#' # 8         20     2   1     0     1     0    0        NA
#' # 9         21     1   0     1     0     0    1         6
#' # 10        21     2   1     0     1     0    1        NA
#' # 11        22     1   0     1     0     0    1        16
#' # 12        22     2   1     0     1     0    1        NA
#' 
#' summarizeAttractors(attr, c("IL2","STAT3","IL2e","IL21e"))
#' 
#' @export
summarizeAttractors <- function(attr, genes=NULL, net=NULL, returnBasin=F) {
    # Default, detect inputs
    if (is.null(genes)) {
        if (is.null(net) ) stop("Error: function requieres a gene list or network")
        genes <- sapply(net$genes, isGeneInput, net)
        genes <- net$genes[genes]
    }
    # Create dataframe
    df <- attractorToDataframe(attr, Boolean = T)
    if (returnBasin) {
      remove <- c("attractor","state","basinSize",genes) 
      } else remove <- c("attractor","state",genes) 
    group.genes <- colnames(df)[!colnames(df) %in% remove]
    # Obtain basins
    if (returnBasin) {
        df[,"basinSize"] <- NA
        basin <- attractorToDataframe(attr)$basinSize
        df$basinSize[df$state==1] <- basin
    } 
    # Remove columns 
    df[,c(genes)]<- NULL
    # Separate cycles
    cycles <- df$attractor[duplicated(df$attractor)]
    steady <- df[!df$attractor %in% cycles,]
    cycles <- df[df$attractor %in% cycles,]
    # Group steady
    if (returnBasin) {
        steady <- steady %>% 
                  group_by_(.dots=group.genes) %>% 
                  summarize(attractor=min(attractor),state=min(state),
                            basinSize=sum(basinSize)) %>%
                  ungroup %>% as.data.frame()
        steady <- steady[,c("attractor","state",group.genes,"basinSize")]
    } else {
        steady <- steady %>% 
                  group_by_(.dots=group.genes) %>% 
                  summarize(attractor=min(attractor),state=min(state)) %>%
                  ungroup %>% as.data.frame()
        steady <- steady[,c("attractor","state",group.genes)]
    }
    
    steady <- steady[order(steady$attractor),]
    # join steady and cycles
    df <- merge(steady, cycles, all=T)
    return(df)
}