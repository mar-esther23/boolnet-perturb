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